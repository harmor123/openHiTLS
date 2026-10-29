/*
 * Copyright (c) 2026 harmor.
 *
 * Constant-time AES adapter for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 *
 * The block core is the project's bitsliced constant-time AES implementation.
 * This adapter exposes the openHiTLS CRYPT_AES API for a CT-only build.
 */
#include "hitls_build.h"
#ifdef HITLS_CRYPTO_AES

#include <stdint.h>
#include <string.h>
#include "crypt_aes.h"
#include "crypt_errno.h"
#include "aes_ct.h"

static int key_bits(uint32_t len)
{
    if (len == 16) return 128;
    if (len == 24) return 192;
    if (len == 32) return 256;
    return 0;
}

static int32_t set_key(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{
    int bits = key_bits(len);
    if (ctx == NULL || key == NULL || bits == 0) return CRYPT_AES_ERR_KEYLEN;
    if (aes_set_encrypt_key(key, bits, &ctx->ct) != 0) return CRYPT_AES_ERR_KEYLEN;
    ctx->rounds = (uint32_t)ctx->ct.rounds;
    return CRYPT_SUCCESS;
}

int32_t CRYPT_AES_SetEncryptKey128(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{ return set_key(ctx, key, len); }
int32_t CRYPT_AES_SetEncryptKey192(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{ return set_key(ctx, key, len); }
int32_t CRYPT_AES_SetEncryptKey256(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{ return set_key(ctx, key, len); }
int32_t CRYPT_AES_SetDecryptKey128(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{ return set_key(ctx, key, len); }
int32_t CRYPT_AES_SetDecryptKey192(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{ return set_key(ctx, key, len); }
int32_t CRYPT_AES_SetDecryptKey256(CRYPT_AES_Key *ctx, const uint8_t *key, uint32_t len)
{ return set_key(ctx, key, len); }

int32_t CRYPT_AES_Encrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out, uint32_t len)
{
    if (ctx == NULL || in == NULL || out == NULL || (len & 15u) != 0) return CRYPT_AES_ERR_KEYLEN;
    aes_ecb_encrypt(out, in, len, &ctx->ct);
    return CRYPT_SUCCESS;
}

int32_t CRYPT_AES_Decrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out, uint32_t len)
{
    if (ctx == NULL || in == NULL || out == NULL || (len & 15u) != 0) return CRYPT_AES_ERR_KEYLEN;
    aes_ecb_decrypt(out, in, len, &ctx->ct);
    return CRYPT_SUCCESS;
}

int32_t CRYPT_AES_ECB_Encrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out, uint32_t len)
{ return CRYPT_AES_Encrypt(ctx, in, out, len); }
int32_t CRYPT_AES_ECB_Decrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out, uint32_t len)
{ return CRYPT_AES_Decrypt(ctx, in, out, len); }

int32_t CRYPT_AES_CBC_Encrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out,
                              uint32_t len, uint8_t *iv)
{
    if (ctx == NULL || in == NULL || out == NULL || iv == NULL) return CRYPT_AES_ERR_KEYLEN;
    aes_cbc_encrypt(out, in, len, iv, &ctx->ct);
    return CRYPT_SUCCESS;
}

int32_t CRYPT_AES_CBC_Decrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out,
                              uint32_t len, uint8_t *iv)
{
    if (ctx == NULL || in == NULL || out == NULL || iv == NULL) return CRYPT_AES_ERR_KEYLEN;
    aes_cbc_decrypt(out, in, len, iv, &ctx->ct);
    return CRYPT_SUCCESS;
}

int32_t CRYPT_AES_CFB_Decrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out,
                              uint32_t len, uint8_t *iv)
{
    if (ctx == NULL || in == NULL || out == NULL || iv == NULL) return CRYPT_AES_ERR_KEYLEN;
    aes_cfb128_decrypt(out, in, len, iv, &ctx->ct);
    return CRYPT_SUCCESS;
}

int32_t CRYPT_AES_CTR_Encrypt(const CRYPT_AES_Key *ctx, const uint8_t *in, uint8_t *out,
                              uint32_t len, uint8_t *iv)
{
    uint8_t block[16];
    uint32_t off = 0;
    int i;

    if (ctx == NULL || in == NULL || out == NULL || iv == NULL) return CRYPT_AES_ERR_KEYLEN;
    while (off < len) {
        uint32_t n = len - off > 16 ? 16 : len - off;
        aes_encrypt_block(block, iv, &ctx->ct);
        for (i = 0; i < (int)n; i++) out[off + i] = in[off + i] ^ block[i];
        for (i = 15; i >= 0; i--) {
            iv[i]++;
            if (iv[i] != 0) break;
        }
        off += n;
    }
    return CRYPT_SUCCESS;
}

void CRYPT_AES_Clean(CRYPT_AES_Key *ctx)
{
    if (ctx != NULL) memset(ctx, 0, sizeof(*ctx));
}

#endif /* HITLS_CRYPTO_AES */
