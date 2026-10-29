/*
 * RISC-V P-256 acceleration adapter for openHiTLS.
 */

#ifndef ECP_NISTP256_ACCEL_H
#define ECP_NISTP256_ACCEL_H

#include "hitls_build.h"
#if defined(HITLS_CRYPTO_CURVE_NISTP256) && defined(HITLS_CRYPTO_NIST_USE_ACCEL) && \
    defined(HITLS_ASM_RISCV64)

#include "ecc_local.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t ECP256_ACCEL_Point2Affine(const ECC_Para *para, ECC_Point *r,
    const ECC_Point *pt);

int32_t ECP256_ACCEL_PointMulAdd(ECC_Para *para, ECC_Point *r,
    const BN_BigNum *k1, const BN_BigNum *k2, const ECC_Point *pt);

int32_t ECP256_ACCEL_PointMul(ECC_Para *para, ECC_Point *r,
    const BN_BigNum *k, const ECC_Point *pt);

#ifdef __cplusplus
}
#endif

#endif /* curve/accel condition */

#endif /* ECP_NISTP256_ACCEL_H */
