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

#include "hitls_build.h"
#if defined(HITLS_CRYPTO_AES) && defined(HITLS_CRYPTO_CTR)

#include "crypt_modes_ctr.h"
#include "modes_local.h"
#if defined(HITLS_CRYPTO_AES_RISCV64) || defined(HITLS_CRYPTO_AES_CT)
#include "crypt_aes.h"
#include "crypt_errno.h"
#endif

int32_t AES_CTR_Update(MODES_CipherCtx *modeCtx, const uint8_t *in, uint32_t inLen, uint8_t *out, uint32_t *outLen)
{
#if defined(HITLS_CRYPTO_AES_RISCV64) || defined(HITLS_CRYPTO_AES_CT)
    if (modeCtx != NULL && modeCtx->commonCtx.ciphCtx != NULL && inLen != 0 &&
        modeCtx->commonCtx.offset == 0 && (inLen & 15u) == 0) {
        int32_t ret = CRYPT_AES_CTR_Encrypt((const CRYPT_AES_Key *)modeCtx->commonCtx.ciphCtx,
                                            in, out, inLen, modeCtx->commonCtx.iv);
        if (ret == CRYPT_SUCCESS) {
            *outLen = inLen;
            return CRYPT_SUCCESS;
        }
    }
#endif
    return MODES_CTR_Update(modeCtx, in, inLen, out, outLen);
}
#endif
