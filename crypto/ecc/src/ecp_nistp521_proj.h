/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */
#ifndef ECP_NISTP521_PROJ_H
#define ECP_NISTP521_PROJ_H

#include "hitls_build.h"
#if defined(HITLS_CRYPTO_CURVE_NISTP521) && defined(HITLS_CRYPTO_NIST_USE_ACCEL)

#include "ecc_local.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t ECP521_PROJ_Point2Affine(const ECC_Para *para, ECC_Point *r, const ECC_Point *pt);

int32_t ECP521_PROJ_PointMulAdd(ECC_Para *para, ECC_Point *r, const BN_BigNum *k1,
    const BN_BigNum *k2, const ECC_Point *pt);

int32_t ECP521_PROJ_PointMul(ECC_Para *para, ECC_Point *r, const BN_BigNum *k,
    const ECC_Point *pt);

#ifdef __cplusplus
}
#endif

#endif /* defined(HITLS_CRYPTO_CURVE_NISTP521) && defined(HITLS_CRYPTO_NIST_USE_ACCEL) */

#endif /* ECP_NISTP521_PROJ_H */
