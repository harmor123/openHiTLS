/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */

/*
 * P-384 specialized point multiplication for openHiTLS.
 *
 * This is a local integration of the project's fixed 6x64-bit CIOS /
 * width-4 table / width-5 wNAF implementation.  The arithmetic core is the
 * template from algorithms/ecc/src/curve_impl.h; this file only converts
 * between openHiTLS BN/ECC_Point objects and big-endian byte arrays.
 */

#include "hitls_build.h"
#if defined(HITLS_CRYPTO_CURVE_NISTP384) && defined(HITLS_CRYPTO_NIST_USE_ACCEL)

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include "bsl_err_internal.h"
#include "bsl_sal.h"
#include "crypt_utils.h"
#include "crypt_errno.h"
#include "crypt_bn.h"
#include "crypt_ecc.h"
#include "ecc_local.h"
#include "ecc_utils.h"

#ifndef HITLS_PLATFORM_INT128
#error "This nistp384 implementation requires 128-bit integer support."
#endif

#define CURVE_PREFIX p384
#define CURVE_LIMBS 6
#define CURVE_BYTES 48
#define CURVE_BITS 384
#define CURVE_WINDOW_BITS 4
#define CURVE_WINDOWS 96

static const uint64_t P[6] = {
    0x00000000ffffffffull, 0xffffffff00000000ull,
    0xfffffffffffffffeull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
};
static const uint64_t N[6] = {
    0xecec196accc52973ull, 0x581a0db248b0a77aull,
    0xc7634d81f4372ddfull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
};
static const uint64_t B[6] = {
    0x2a85c8edd3ec2aefull, 0xc656398d8a2ed19dull,
    0x0314088f5013875aull, 0x181d9c6efe814112ull,
    0x988e056be3f82d19ull, 0xb3312fa7e23ee7e4ull,
};
static const uint64_t GX[6] = {
    0x3a545e3872760ab7ull, 0x5502f25dbf55296cull,
    0x59f741e082542a38ull, 0x6e1d3b628ba79b98ull,
    0x8eb1c71ef320ad74ull, 0xaa87ca22be8b0537ull,
};
static const uint64_t GY[6] = {
    0x7a431d7c90ea0e5full, 0x0a60b1ce1d7e819dull,
    0xe9da3113b5f0b8c0ull, 0xf8f41dbd289a147cull,
    0x5d9e98bf9292dc29ull, 0x3617de4a96262c6full,
};
static const uint64_t P_MINUS_2[6] = {
    0x00000000fffffffdull, 0xffffffff00000000ull,
    0xfffffffffffffffeull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
};
static const uint64_t N_MINUS_2[6] = {
    0xecec196accc52971ull, 0x581a0db248b0a77aull,
    0xc7634d81f4372ddfull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
};
static const uint64_t P_R[6] = {
    0xffffffff00000001ull, 0x00000000ffffffffull,
    0x0000000000000001ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
};
static const uint64_t P_R2[6] = {
    0xfffffffe00000001ull, 0x0000000200000000ull,
    0xfffffffe00000000ull, 0x0000000200000000ull,
    0x0000000000000001ull, 0x0000000000000000ull,
};
static const uint64_t N_R[6] = {
    0x1313e695333ad68dull, 0xa7e5f24db74f5885ull,
    0x389cb27e0bc8d220ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
};
static const uint64_t N_R2[6] = {
    0x2d319b2419b409a9ull, 0xff3d81e5df1aa419ull,
    0xbc3e483afcb82947ull, 0xd40d49174aab1cc5ull,
    0x3fb05b7a28266895ull, 0x0c84ee012b39bf21ull,
};

#define P_N0 0x0000000100000001ull
#define N_N0 0x6ed46089e88fdc45ull

#include "curve_p384_impl.h"

static int32_t P384PointMulBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
    const uint8_t k[CURVE_BYTES], const uint8_t *px, const uint8_t *py)
{
    sc kv;
    jac_point r;
    fe x;
    fe y;

    u_from_bytes(&kv, k);
    init_constants();
    if (px == NULL) {
        scalar_mul_table(&r, g_table, kv.v);
    } else {
        jac_point q;
        jac_point table[CURVE_TABLE_SIZE];

        if (parse_point(&q, px, py) != 0) {
            return CRYPT_ECC_POINT_NOT_ON_CURVE;
        }
        build_table(table, &q);
        scalar_mul_table(&r, table, kv.v);
    }
    if (point_is_infinity(&r)) {
        return CRYPT_ECC_POINT_AT_INFINITY;
    }
    point_to_affine(&x, &y, &r);
    fe_from_mont(&x, &x);
    fe_from_mont(&y, &y);
    u_to_bytes(qx, &x);
    u_to_bytes(qy, &y);
    return CRYPT_SUCCESS;
}

static int32_t P384PointMulAddBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
    const uint8_t k1[CURVE_BYTES], const uint8_t k2[CURVE_BYTES],
    const uint8_t px[CURVE_BYTES], const uint8_t py[CURVE_BYTES])
{
    sc k1v;
    sc k2v;
    jac_point q;
    jac_point xp;
    affine_point qtab[8];
    int8_t naf1[CURVE_BITS + 2];
    int8_t naf2[CURVE_BITS + 2];
    int len1;
    int len2;
    int len;
    int i;
    fe x;
    fe y;

    u_from_bytes(&k1v, k1);
    u_from_bytes(&k2v, k2);
    init_constants();
    if (parse_point(&q, px, py) != 0) {
        return CRYPT_ECC_POINT_NOT_ON_CURVE;
    }
    len1 = scalar_wnaf(naf1, k1v.v, 5);
    len2 = scalar_wnaf(naf2, k2v.v, 5);
    build_naf_table(qtab, &q);
    len = len1 > len2 ? len1 : len2;
    point_set_infinity(&xp);
    for (i = len - 1; i >= 0; i--) {
        if (!point_is_infinity(&xp)) {
            point_double(&xp, &xp);
        }
        if (i < len1 && naf1[i] != 0) {
            int d = naf1[i];
            int idx = ((d < 0 ? -d : d) - 1) / 2;

            point_add_signed(&xp, &g_naf[idx], d < 0);
        }
        if (i < len2 && naf2[i] != 0) {
            int d = naf2[i];
            int idx = ((d < 0 ? -d : d) - 1) / 2;

            point_add_signed(&xp, &qtab[idx], d < 0);
        }
    }
    if (point_is_infinity(&xp)) {
        return CRYPT_ECC_POINT_AT_INFINITY;
    }
    point_to_affine(&x, &y, &xp);
    fe_from_mont(&x, &x);
    fe_from_mont(&y, &y);
    u_to_bytes(qx, &x);
    u_to_bytes(qy, &y);
    return CRYPT_SUCCESS;
}

static int32_t P384Point2AffineBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
    const uint8_t px[CURVE_BYTES], const uint8_t py[CURVE_BYTES],
    const uint8_t pz[CURVE_BYTES])
{
    fe x;
    fe y;
    fe z;
    fe ox;
    fe oy;
    jac_point p;

    if (fe_parse(&x, px) != 0 || fe_parse(&y, py) != 0 || fe_parse(&z, pz) != 0) {
        return CRYPT_ECC_POINT_NOT_ON_CURVE;
    }
    fe_to_mont(&p.x, &x);
    fe_to_mont(&p.y, &y);
    fe_to_mont(&p.z, &z);
    point_to_affine(&ox, &oy, &p);
    fe_from_mont(&ox, &ox);
    fe_from_mont(&oy, &oy);
    u_to_bytes(qx, &ox);
    u_to_bytes(qy, &oy);
    return CRYPT_SUCCESS;
}

static int32_t PointToBytes(const ECC_Point *pt, uint8_t x[CURVE_BYTES],
    uint8_t y[CURVE_BYTES], uint8_t z[CURVE_BYTES])
{
    int32_t ret;

    ret = BN_Bn2BinFixZero(&pt->x, x, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = BN_Bn2BinFixZero(&pt->y, y, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    return BN_Bn2BinFixZero(&pt->z, z, CURVE_BYTES);
}

static int32_t BytesToPoint(ECC_Point *pt, const uint8_t x[CURVE_BYTES],
    const uint8_t y[CURVE_BYTES], const uint8_t z[CURVE_BYTES])
{
    int32_t ret;

    ret = BN_Bin2Bn(&pt->x, x, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = BN_Bin2Bn(&pt->y, y, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    return BN_Bin2Bn(&pt->z, z, CURVE_BYTES);
}

static int32_t PointToAffineBytes(const ECC_Point *pt, uint8_t x[CURVE_BYTES],
    uint8_t y[CURVE_BYTES])
{
    uint8_t px[CURVE_BYTES];
    uint8_t py[CURVE_BYTES];
    uint8_t pz[CURVE_BYTES];
    uint8_t ax[CURVE_BYTES];
    uint8_t ay[CURVE_BYTES];
    int32_t ret;
    int32_t isOne = 1;
    int32_t i;

    ret = PointToBytes(pt, px, py, pz);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    for (i = 0; i < CURVE_BYTES - 1; i++) {
        isOne &= (pz[i] == 0);
    }
    isOne &= (pz[CURVE_BYTES - 1] == 1);
    if (isOne) {
        memcpy(x, px, CURVE_BYTES);
        memcpy(y, py, CURVE_BYTES);
        return CRYPT_SUCCESS;
    }
    ret = P384Point2AffineBytes(ax, ay, px, py, pz);
    if (ret == CRYPT_SUCCESS) {
        memcpy(x, ax, CURVE_BYTES);
        memcpy(y, ay, CURVE_BYTES);
    }
    return ret;
}

int32_t ECP384_PointMul(ECC_Para *para, ECC_Point *r, const BN_BigNum *k,
    const ECC_Point *pt)
{
    (void)para;
    uint8_t kb[CURVE_BYTES];
    uint8_t px[CURVE_BYTES];
    uint8_t py[CURVE_BYTES];
    uint8_t qx[CURVE_BYTES];
    uint8_t qy[CURVE_BYTES];
    uint8_t qz[CURVE_BYTES] = {0};
    int32_t ret;

    ret = CheckBnValid(k, CURVE_BITS);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = BN_Bn2BinFixZero(k, kb, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    if (pt == NULL) {
        ret = P384PointMulBytes(qx, qy, kb, NULL, NULL);
    } else {
        ret = PointToAffineBytes(pt, px, py);
        if (ret == CRYPT_SUCCESS) {
            ret = P384PointMulBytes(qx, qy, kb, px, py);
        }
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    BSL_SAL_CleanseData(kb, sizeof(kb));
    return ret;
}

int32_t ECP384_PointMulAdd(ECC_Para *para, ECC_Point *r, const BN_BigNum *k1,
    const BN_BigNum *k2, const ECC_Point *pt)
{
    (void)para;
    uint8_t k1b[CURVE_BYTES];
    uint8_t k2b[CURVE_BYTES];
    uint8_t px[CURVE_BYTES];
    uint8_t py[CURVE_BYTES];
    uint8_t qx[CURVE_BYTES];
    uint8_t qy[CURVE_BYTES];
    uint8_t qz[CURVE_BYTES] = {0};
    int32_t ret;

    ret = CheckBnValid(k1, CURVE_BITS);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = CheckBnValid(k2, CURVE_BITS);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = BN_Bn2BinFixZero(k1, k1b, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = BN_Bn2BinFixZero(k2, k2b, CURVE_BYTES);
    if (ret != CRYPT_SUCCESS) {
        return ret;
    }
    ret = PointToAffineBytes(pt, px, py);
    if (ret == CRYPT_SUCCESS) {
        ret = P384PointMulAddBytes(qx, qy, k1b, k2b, px, py);
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    BSL_SAL_CleanseData(k1b, sizeof(k1b));
    BSL_SAL_CleanseData(k2b, sizeof(k2b));
    return ret;
}

int32_t ECP384_Point2Affine(const ECC_Para *para, ECC_Point *r, const ECC_Point *pt)
{
    (void)para;
    uint8_t px[CURVE_BYTES];
    uint8_t py[CURVE_BYTES];
    uint8_t pz[CURVE_BYTES];
    uint8_t qx[CURVE_BYTES];
    uint8_t qy[CURVE_BYTES];
    uint8_t qz[CURVE_BYTES] = {0};
    int32_t ret;

    ret = PointToBytes(pt, px, py, pz);
    if (ret == CRYPT_SUCCESS) {
        ret = P384Point2AffineBytes(qx, qy, px, py, pz);
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    return ret;
}

#endif /* defined(HITLS_CRYPTO_CURVE_NISTP384) && defined(HITLS_CRYPTO_NIST_USE_ACCEL) */
