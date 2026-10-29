/*
 * ChaCha20 multi-block update for RISC-V RVV 1.0 (VLEN >= 256).
 *
 * Selected instead of chacha20block.c when HITLS_CRYPTO_CHACHA20_RISCV64
 * is enabled.  Full 8-block groups use the vector kernel; large inputs use
 * the nine-block vector+scalar paired kernel.  Blocks not covered by the
 * vector kernels fall back to CHACHA20_Block from chacha20.c.
 */

#include "hitls_build.h"
#ifdef HITLS_CRYPTO_CHACHA20

#include <stdint.h>
#include "crypt_utils.h"
#include "chacha20_local.h"

void chacha20_rvv_blocks(uint8_t *out, const uint8_t *in, size_t groups,
                         const uint32_t key[8], uint32_t counter,
                         const uint8_t nonce[12]);
void chacha20_rvv_blocks9(uint8_t *out, const uint8_t *in, size_t groups,
                          const uint32_t key[8], uint32_t counter,
                          const uint8_t nonce[12]);

static int ChaCha20HasV256(void)
{
    unsigned long vlenb;

    if (!IsSupportV()) {
        return 0;
    }
    __asm__ volatile ("csrr %0, vlenb" : "=r"(vlenb));
    return vlenb >= 32;
}

void CHACHA20_Update(CRYPT_CHACHA20_Ctx *ctx, const uint8_t *in,
                     uint8_t *out, uint32_t len)
{
    uint32_t counter = ctx->state[12];
    uint32_t off = 0;

    /* vlsseg8e32 requires naturally aligned element accesses.  Callers may
     * hand us arbitrary byte pointers; use the C path for those. */
    if (ChaCha20HasV256() &&
        (((uintptr_t)in | (uintptr_t)out) & (sizeof(uint32_t) - 1u)) == 0) {
        const uint32_t *key = ctx->state + 4;
        const uint8_t *nonce = (const uint8_t *)(ctx->state + 13);

        if (len >= 9 * 64 + 8 * 64) {
            size_t groups = (len - 8 * 64) / (9 * 64);
            uint64_t room = ((uint64_t)0xffffffffu - counter + 1) / 9;

            if ((uint64_t)groups > room) {
                groups = (size_t)room;
            }
            if (groups != 0) {
                chacha20_rvv_blocks9(out, in, groups, key, counter, nonce);
                off = (uint32_t)(groups * 9 * 64);
                counter += (uint32_t)(groups * 9);
            }
        }
        while (len - off >= 8 * 64) {
            size_t groups = (len - off) / (8 * 64);
            uint64_t room = ((uint64_t)0xffffffffu - counter + 1) / 8;

            if ((uint64_t)groups > room) {
                groups = (size_t)room;
            }
            if (groups == 0) {
                break;
            }
            chacha20_rvv_blocks(out + off, in + off, groups, key, counter,
                                nonce);
            off += (uint32_t)(groups * 8 * 64);
            counter += (uint32_t)(groups * 8);
        }
        ctx->state[12] = counter;
    }

    while (off + 64 <= len) {
        CHACHA20_Block(ctx);
        DATA64_XOR(ctx->last.u, in + off, out + off, 64);
        off += 64;
    }
}

#endif /* HITLS_CRYPTO_CHACHA20 */
