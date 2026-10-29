/*
 * P-521 specialized point multiplication for openHiTLS.
 *
 * This is a local integration of the project's fixed 9x64-bit CIOS /
 * width-4 table / width-5 wNAF implementation.  The arithmetic core is the
 * template from algorithms/ecc/src/curve_impl.h; this file only converts
 * between openHiTLS BN/ECC_Point objects and big-endian byte arrays.
 */

#include "hitls_build.h"
#if defined(HITLS_CRYPTO_CURVE_NISTP521) && defined(HITLS_CRYPTO_NIST_USE_ACCEL)

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
#error "This nistp521 implementation requires 128-bit integer support."
#endif

#define CURVE_PREFIX p521
#define CURVE_LIMBS 9
#define CURVE_BYTES 66
#define CURVE_BITS 521
#define CURVE_WINDOW_BITS 4
#define CURVE_WINDOWS 131

static const uint64_t P[9] = {
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0x00000000000001ffull,
};
static const uint64_t N[9] = {
    0xbb6fb71e91386409ull, 0x3bb5c9b8899c47aeull,
    0x7fcc0148f709a5d0ull, 0x51868783bf2f966bull,
    0xfffffffffffffffaull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0x00000000000001ffull,
};
static const uint64_t B[9] = {
    0xef451fd46b503f00ull, 0x3573df883d2c34f1ull,
    0x1652c0bd3bb1bf07ull, 0x56193951ec7e937bull,
    0xb8b489918ef109e1ull, 0xa2da725b99b315f3ull,
    0x929a21a0b68540eeull, 0x953eb9618e1c9a1full,
    0x0000000000000051ull,
};
static const uint64_t GX[9] = {
    0xf97e7e31c2e5bd66ull, 0x3348b3c1856a429bull,
    0xfe1dc127a2ffa8deull, 0xa14b5e77efe75928ull,
    0xf828af606b4d3dbaull, 0x9c648139053fb521ull,
    0x9e3ecb662395b442ull, 0x858e06b70404e9cdull,
    0x00000000000000c6ull,
};
static const uint64_t GY[9] = {
    0x88be94769fd16650ull, 0x353c7086a272c240ull,
    0xc550b9013fad0761ull, 0x97ee72995ef42640ull,
    0x17afbd17273e662cull, 0x98f54449579b4468ull,
    0x5c8a5fb42c7d1bd9ull, 0x39296a789a3bc004ull,
    0x0000000000000118ull,
};
static const uint64_t P_MINUS_2[9] = {
    0xfffffffffffffffdull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0x00000000000001ffull,
};
static const uint64_t N_MINUS_2[9] = {
    0xbb6fb71e91386407ull, 0x3bb5c9b8899c47aeull,
    0x7fcc0148f709a5d0ull, 0x51868783bf2f966bull,
    0xfffffffffffffffaull, 0xffffffffffffffffull,
    0xffffffffffffffffull, 0xffffffffffffffffull,
    0x00000000000001ffull,
};
static const uint64_t P_R[9] = {
    0x0080000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull,
};
static const uint64_t P_R2[9] = {
    0x0000000000000000ull, 0x0000400000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull,
};
static const uint64_t N_R[9] = {
    0xfb80000000000000ull, 0x28a2482470b763cdull,
    0x17e2251b23bb31dcull, 0xca4019ff5b847b2dull,
    0x02d73cbc3e206834ull, 0x0000000000000000ull,
    0x0000000000000000ull, 0x0000000000000000ull,
    0x0000000000000000ull,
};
static const uint64_t N_R2[9] = {
    0x137cd04dcf15dd04ull, 0xf707badce5547ea3ull,
    0x12a78d38794573ffull, 0xd3721ef557f75e06ull,
    0xdd6e23d82e49c7dbull, 0xcff3d142b7756e3eull,
    0x5bcc6d61a8e567bcull, 0x2d8e03d1492d0d45ull,
    0x000000000000003dull,
};

#define P_N0 0x0000000000000001ull
#define N_N0 0x1d2f5ccd79a995c7ull

#include "curve_p521_proj_impl.h"

static int32_t P521PointMulBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
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

static int32_t P521PointMulAddBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
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

static int32_t P521Point2AffineBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
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
    ret = P521Point2AffineBytes(ax, ay, px, py, pz);
    if (ret == CRYPT_SUCCESS) {
        memcpy(x, ax, CURVE_BYTES);
        memcpy(y, ay, CURVE_BYTES);
    }
    return ret;
}

int32_t ECP521_PROJ_PointMul(ECC_Para *para, ECC_Point *r, const BN_BigNum *k,
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
        ret = P521PointMulBytes(qx, qy, kb, NULL, NULL);
    } else {
        ret = PointToAffineBytes(pt, px, py);
        if (ret == CRYPT_SUCCESS) {
            ret = P521PointMulBytes(qx, qy, kb, px, py);
        }
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    BSL_SAL_CleanseData(kb, sizeof(kb));
    return ret;
}

int32_t ECP521_PROJ_PointMulAdd(ECC_Para *para, ECC_Point *r, const BN_BigNum *k1,
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
        ret = P521PointMulAddBytes(qx, qy, k1b, k2b, px, py);
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    BSL_SAL_CleanseData(k1b, sizeof(k1b));
    BSL_SAL_CleanseData(k2b, sizeof(k2b));
    return ret;
}

int32_t ECP521_PROJ_Point2Affine(const ECC_Para *para, ECC_Point *r, const ECC_Point *pt)
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
        ret = P521Point2AffineBytes(qx, qy, px, py, pz);
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    return ret;
}

#endif /* defined(HITLS_CRYPTO_CURVE_NISTP521) && defined(HITLS_CRYPTO_NIST_USE_ACCEL) */
