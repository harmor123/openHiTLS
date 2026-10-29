/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */

#ifndef P256_H
#define P256_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * P-256 (secp256r1) primitives.  All byte strings are big-endian and fixed
 * width: scalars/coordinates are 32 bytes.
 *
 * Return conventions:
 *   0  success
 *  -1  invalid input or operation failure (e.g. point at infinity)
 *   1  signature/public-key valid (verify/validate only)
 *   0  signature/public-key invalid (verify/validate only)
 */

int p256_public_key(uint8_t qx[32], uint8_t qy[32], const uint8_t d[32]);

int p256_validate_public_key(const uint8_t qx[32], const uint8_t qy[32]);

/* z is the (already truncated) message digest as a big-endian integer. */
int p256_ecdsa_sign(uint8_t r[32], uint8_t s[32], const uint8_t d[32],
                    const uint8_t z[32], const uint8_t k[32]);

int p256_ecdsa_verify(const uint8_t qx[32], const uint8_t qy[32],
                      const uint8_t z[32], const uint8_t r[32],
                      const uint8_t s[32]);

int p256_ecdh(uint8_t out[32], const uint8_t d[32], const uint8_t qx[32],
              const uint8_t qy[32]);

#ifdef __cplusplus
}
#endif

#endif /* P256_H */
