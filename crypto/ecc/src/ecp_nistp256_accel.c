/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */

/*
 * RISC-V P-256 acceleration adapter for openHiTLS.
 *
 * The arithmetic and point implementation is the project's P-256 source;
 * P256_USE_ASM selects p256_mont_mul_rv64.S for the field multiplication.
 */

#include "hitls_build.h"
#if defined(HITLS_CRYPTO_CURVE_NISTP256) && defined(HITLS_CRYPTO_NIST_USE_ACCEL) && \
    defined(HITLS_ASM_RISCV64)

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
#error "This nistp256 RISC-V implementation requires 128-bit integer support."
#endif

#ifndef P256_NO_ASM
#define P256_USE_ASM 1
#endif
#include "curve_init_once.h"
#include "p256_accel_impl.h"

#define CURVE_BYTES 32
#define CURVE_BITS 256

static int32_t P256PointMulBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
    const uint8_t k[CURVE_BYTES], const uint8_t *px, const uint8_t *py)
{
    sc kv;
    jac_point r;
    fe x;
    fe y;

    u256_from_bytes(&kv, k);
    init_constants();
    if (px == NULL) {
        scalar_mul_table(&r, g_table, kv.v);
    } else {
        jac_point q;
        jac_point table[16];

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
    u256_to_bytes(qx, &x);
    u256_to_bytes(qy, &y);
    return CRYPT_SUCCESS;
}

static int32_t P256PointMulAddBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
    const uint8_t k1[CURVE_BYTES], const uint8_t k2[CURVE_BYTES],
    const uint8_t px[CURVE_BYTES], const uint8_t py[CURVE_BYTES])
{
    sc k1v;
    sc k2v;
    jac_point q;
    jac_point qt[16];
    jac_point r;
    fe x;
    fe y;

    u256_from_bytes(&k1v, k1);
    u256_from_bytes(&k2v, k2);
    init_constants();
    if (parse_point(&q, px, py) != 0) {
        return CRYPT_ECC_POINT_NOT_ON_CURVE;
    }
    build_table(qt, &q);
    scalar_mul_joint(&r, g_table, k1v.v, qt, k2v.v);
    if (point_is_infinity(&r)) {
        return CRYPT_ECC_POINT_AT_INFINITY;
    }
    point_to_affine(&x, &y, &r);
    fe_from_mont(&x, &x);
    fe_from_mont(&y, &y);
    u256_to_bytes(qx, &x);
    u256_to_bytes(qy, &y);
    return CRYPT_SUCCESS;
}

static int32_t P256Point2AffineBytes(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
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
    u256_to_bytes(qx, &ox);
    u256_to_bytes(qy, &oy);
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
    ret = P256Point2AffineBytes(ax, ay, px, py, pz);
    if (ret == CRYPT_SUCCESS) {
        memcpy(x, ax, CURVE_BYTES);
        memcpy(y, ay, CURVE_BYTES);
    }
    return ret;
}

int32_t ECP256_ACCEL_PointMul(ECC_Para *para, ECC_Point *r, const BN_BigNum *k,
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
        ret = P256PointMulBytes(qx, qy, kb, NULL, NULL);
    } else {
        ret = PointToAffineBytes(pt, px, py);
        if (ret == CRYPT_SUCCESS) {
            ret = P256PointMulBytes(qx, qy, kb, px, py);
        }
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    BSL_SAL_CleanseData(kb, sizeof(kb));
    return ret;
}

int32_t ECP256_ACCEL_PointMulAdd(ECC_Para *para, ECC_Point *r,
    const BN_BigNum *k1, const BN_BigNum *k2, const ECC_Point *pt)
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
        ret = P256PointMulAddBytes(qx, qy, k1b, k2b, px, py);
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    BSL_SAL_CleanseData(k1b, sizeof(k1b));
    BSL_SAL_CleanseData(k2b, sizeof(k2b));
    return ret;
}

int32_t ECP256_ACCEL_Point2Affine(const ECC_Para *para, ECC_Point *r,
    const ECC_Point *pt)
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
        ret = P256Point2AffineBytes(qx, qy, px, py, pz);
    }
    if (ret == CRYPT_SUCCESS) {
        qz[CURVE_BYTES - 1] = 1;
        ret = BytesToPoint(r, qx, qy, qz);
    }
    return ret;
}

#endif /* curve/accel condition */
