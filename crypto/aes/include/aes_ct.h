/*
 * Copyright (c) 2026 harmor.
 *
 * Constant-time AES interface used by the openHiTLS AES/CT build.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */
#ifndef AES_CT_H
#define AES_CT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AES_MAX_ROUND_KEYS 60

typedef struct {
    uint32_t rk[AES_MAX_ROUND_KEYS];   /* big-endian AES words, C path */
    uint32_t rk_z[AES_MAX_ROUND_KEYS]; /* byte-swapped for Zvkned vle32 */
    uint32_t rk_d[AES_MAX_ROUND_KEYS]; /* equivalent inverse cipher keys */
    uint8_t rk_b[AES_MAX_ROUND_KEYS * 4]
        __attribute__((aligned(16)));  /* round-key bytes in state order */
    int rounds;
} aes_key;

int aes_set_encrypt_key(const uint8_t *key, int bits, aes_key *ctx);

/* Runtime C/Zvkned selection.  Returns 0 on success, -1 when the requested
 * implementation was not compiled in.  Default is Zvkned when available,
 * otherwise the constant-time C fallback. */
int aes_set_impl(int use_asm);
int aes_get_impl(void);

/* Optional RVV 1.0 block implementation (VLEN >= 256).  Returns 0 on
 * success, -1 when RVV is unavailable or too narrow. */
int aes_set_rvv(int use_rvv);
int aes_get_rvv(void);
int aes_rvv_available(void);

void aes_encrypt_block_c(uint8_t out[16], const uint8_t in[16],
                         const aes_key *ctx);
void aes_encrypt_block(uint8_t out[16], const uint8_t in[16],
                       const aes_key *ctx);

void aes_decrypt_block_c(uint8_t out[16], const uint8_t in[16],
                         const aes_key *ctx);
void aes_decrypt_block(uint8_t out[16], const uint8_t in[16],
                       const aes_key *ctx);

void aes_ecb_encrypt(uint8_t *out, const uint8_t *in, size_t len,
                     const aes_key *ctx);
void aes_cbc_encrypt(uint8_t *out, const uint8_t *in, size_t len,
                     uint8_t iv[16], const aes_key *ctx);
void aes_ecb_decrypt(uint8_t *out, const uint8_t *in, size_t len,
                     const aes_key *ctx);
void aes_cbc_decrypt(uint8_t *out, const uint8_t *in, size_t len,
                     uint8_t iv[16], const aes_key *ctx);
void aes_cfb128_encrypt(uint8_t *out, const uint8_t *in, size_t len,
                        uint8_t iv[16], const aes_key *ctx);
void aes_cfb128_decrypt(uint8_t *out, const uint8_t *in, size_t len,
                        uint8_t iv[16], const aes_key *ctx);
void aes_cfb8_crypt(uint8_t *out, const uint8_t *in, size_t len,
                    uint8_t iv[16], const aes_key *ctx, int encrypt);
void aes_cfb1_crypt(uint8_t *out_bits, const uint8_t *in_bits, size_t nbits,
                    uint8_t iv[16], const aes_key *ctx, int encrypt);
void aes_ofb_crypt(uint8_t *out, const uint8_t *in, size_t len,
                   uint8_t iv[16], const aes_key *ctx);

#if defined(__riscv) && defined(AES_USE_RVV)
void aes_encrypt_block_rvv(uint8_t out[16], const uint8_t in[16],
                           const uint8_t *rk_b, int rounds);
void aes_decrypt_block_rvv(uint8_t out[16], const uint8_t in[16],
                           const uint8_t *rk_b, int rounds);
#endif

#if defined(__riscv) && defined(AES_USE_ZVKNED)
void aes_encrypt_block_zvkned(uint8_t out[16], const uint8_t in[16],
                              const uint32_t *rk, int rounds);
void aes_decrypt_block_zvkned(uint8_t out[16], const uint8_t in[16],
                              const uint32_t *rk, int rounds);
void aes_ecb_encrypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                                const uint32_t *rk);
void aes_cbc_encrypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                                uint8_t iv[16], const uint32_t *rk);
void aes_ecb_encrypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                                const uint32_t *rk);
void aes_cbc_encrypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                                uint8_t iv[16], const uint32_t *rk);
void aes_ecb_encrypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                                const uint32_t *rk);
void aes_cbc_encrypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                                uint8_t iv[16], const uint32_t *rk);
void aes_ecb_decrypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                                const uint32_t *rk);
void aes_cbc_decrypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                                uint8_t iv[16], const uint32_t *rk);
void aes_ecb_decrypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                                const uint32_t *rk);
void aes_cbc_decrypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                                uint8_t iv[16], const uint32_t *rk);
void aes_ecb_decrypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                                const uint32_t *rk);
void aes_cbc_decrypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                                uint8_t iv[16], const uint32_t *rk);
void aes_cfb128_encrypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb128_encrypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb128_encrypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb128_decrypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb128_decrypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb128_decrypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_ofb_crypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                              uint8_t iv[16], const uint32_t *rk);
void aes_ofb_crypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                              uint8_t iv[16], const uint32_t *rk);
void aes_ofb_crypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                              uint8_t iv[16], const uint32_t *rk);
void aes_cfb8_crypt_zvkned_128(uint8_t *out, const uint8_t *in, size_t len,
                               uint8_t iv[16], const uint32_t *rk);
void aes_cfb8_crypt_zvkned_128_dec(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb8_crypt_zvkned_192(uint8_t *out, const uint8_t *in, size_t len,
                               uint8_t iv[16], const uint32_t *rk);
void aes_cfb8_crypt_zvkned_192_dec(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
void aes_cfb8_crypt_zvkned_256(uint8_t *out, const uint8_t *in, size_t len,
                               uint8_t iv[16], const uint32_t *rk);
void aes_cfb8_crypt_zvkned_256_dec(uint8_t *out, const uint8_t *in, size_t len,
                                   uint8_t iv[16], const uint32_t *rk);
#endif

#ifdef __cplusplus
}
#endif

#endif /* AES_CT_H */
