/*
 * This file is part of the openHiTLS project.
 *
 * openHiTLS is licensed under the Mulan PSL v2.
 * You can use this software according to the terms and conditions of the Mulan PSL v2.
 * You may obtain a copy of Mulan PSL v2 at:
 *
 *     http://license.coscl.org.cn/MulanPSL2
 *
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
 * EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
 * MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 * See the Mulan PSL v2 for more details.
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
