/*
 * AES-128/192/256 block functions using RVV 1.0.
 *
 * The target is the MUSE Pi Pro M1 (VLEN=256, no Zvkned).  The S-box is
 * looked up from eight 32-byte vectors held in registers; the secret byte
 * only feeds vrgather indices, never a memory address.  ShiftRows and the
 * MixColumns rotations use fixed permutation vectors, and xtime is
 * branchless.  VLEN must be at least 256; narrower implementations fall
 * back to the constant-time C path.
 */

#include "aes_ct.h"

#if defined(__riscv) && defined(AES_USE_RVV)

#include <riscv_vector.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define RVV_VL 16
#define RVV_TABLE_BYTES 32

static const uint8_t rvv_sbox[256] __attribute__((aligned(32))) = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
    0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
    0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
    0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
    0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
    0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
    0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
    0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
    0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
    0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
    0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
    0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16,
};

static const uint8_t rvv_inv_sbox[256] __attribute__((aligned(32))) = {
    0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38,
    0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
    0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87,
    0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
    0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d,
    0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
    0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2,
    0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
    0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16,
    0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
    0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda,
    0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
    0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a,
    0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
    0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02,
    0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
    0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea,
    0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
    0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85,
    0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
    0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89,
    0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
    0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20,
    0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
    0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31,
    0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
    0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d,
    0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
    0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0,
    0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26,
    0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d,
};

static const uint8_t rvv_shift_idx[16] __attribute__((aligned(16))) = {
    0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11,
};

static const uint8_t rvv_inv_shift_idx[16] __attribute__((aligned(16))) = {
    0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3,
};

static const uint8_t rvv_rot1_idx[16] __attribute__((aligned(16))) = {
    1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12,
};

static const uint8_t rvv_rot2_idx[16] __attribute__((aligned(16))) = {
    2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13,
};

static const uint8_t rvv_rot3_idx[16] __attribute__((aligned(16))) = {
    3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14,
};

static const uint8_t rvv_inv_shift_idx2[32] __attribute__((aligned(32))) = {
    0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3,
    16, 29, 26, 23, 20, 17, 30, 27, 24, 21, 18, 31, 28, 25, 22, 19,
};
static const uint8_t rvv_rot1_idx2[32] __attribute__((aligned(32))) = {
    1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12,
    17, 18, 19, 16, 21, 22, 23, 20, 25, 26, 27, 24, 29, 30, 31, 28,
};
static const uint8_t rvv_rot2_idx2[32] __attribute__((aligned(32))) = {
    2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13,
    18, 19, 16, 17, 22, 23, 20, 21, 26, 27, 24, 25, 30, 31, 28, 29,
};
static const uint8_t rvv_rot3_idx2[32] __attribute__((aligned(32))) = {
    3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14,
    19, 16, 17, 18, 23, 20, 21, 22, 27, 24, 25, 26, 31, 28, 29, 30,
};

static const uint8_t rvv_shift_idx2[32] __attribute__((aligned(32))) = {
    0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11,
    16, 21, 26, 31, 20, 25, 30, 19, 24, 29, 18, 23, 28, 17, 22, 27,
};

int aes_rvv_available(void)
{
    unsigned long vlenb;

    __asm__ volatile ("csrr %0, vlenb" : "=r"(vlenb));
    return vlenb >= RVV_TABLE_BYTES;
}

static inline __attribute__((always_inline))
vuint8m1_t rvv_sbox_lookup(vuint8m1_t x, size_t vl,
                           vuint8m1_t t0, vuint8m1_t t1,
                           vuint8m1_t t2, vuint8m1_t t3,
                           vuint8m1_t t4, vuint8m1_t t5,
                           vuint8m1_t t6, vuint8m1_t t7)
{
    vuint8m1_t idx = __riscv_vand_vx_u8m1(x, 31, vl);
    vuint8m1_t hi = __riscv_vsrl_vx_u8m1(x, 5, vl);
    vuint8m1_t r = __riscv_vrgather_vv_u8m1(t0, idx, vl);

#define RVV_SBOX_MERGE(tab, g)                                            \
    r = __riscv_vmerge_vvm_u8m1(                                          \
        r, __riscv_vrgather_vv_u8m1((tab), idx, vl),                      \
        __riscv_vmseq_vx_u8m1_b8(hi, (g), vl), vl)

    RVV_SBOX_MERGE(t1, 1);
    RVV_SBOX_MERGE(t2, 2);
    RVV_SBOX_MERGE(t3, 3);
    RVV_SBOX_MERGE(t4, 4);
    RVV_SBOX_MERGE(t5, 5);
    RVV_SBOX_MERGE(t6, 6);
    RVV_SBOX_MERGE(t7, 7);

#undef RVV_SBOX_MERGE

    return r;
}

static inline __attribute__((always_inline))
vuint8m1_t rvv_xtime(vuint8m1_t x, size_t vl)
{
    vuint8m1_t doubled = __riscv_vadd_vv_u8m1(x, x, vl);
    vbool8_t hi = __riscv_vmsgtu_vx_u8m1_b8(x, 0x7f, vl);
    vuint8m1_t red = __riscv_vmerge_vxm_u8m1(
        __riscv_vmv_v_x_u8m1(0, vl), 0x1b, hi, vl);

    return __riscv_vxor_vv_u8m1(doubled, red, vl);
}

static inline __attribute__((always_inline))
vuint8m1_t rvv_mix(vuint8m1_t x, size_t vl,
                    vuint8m1_t i1, vuint8m1_t i2, vuint8m1_t i3)
{
    vuint8m1_t r1 = __riscv_vrgather_vv_u8m1(x, i1, vl);
    vuint8m1_t r2 = __riscv_vrgather_vv_u8m1(x, i2, vl);
    vuint8m1_t r3 = __riscv_vrgather_vv_u8m1(x, i3, vl);
    vuint8m1_t t = rvv_xtime(__riscv_vxor_vv_u8m1(x, r1, vl), vl);

    return __riscv_vxor_vv_u8m1(
        t, __riscv_vxor_vv_u8m1(r1, __riscv_vxor_vv_u8m1(r2, r3, vl), vl),
        vl);
}

static inline __attribute__((always_inline))
vuint8m1_t rvv_inv_mix(vuint8m1_t x, size_t vl,
                        vuint8m1_t i1, vuint8m1_t i2, vuint8m1_t i3)
{
    vuint8m1_t r1 = __riscv_vrgather_vv_u8m1(x, i1, vl);
    vuint8m1_t r2 = __riscv_vrgather_vv_u8m1(x, i2, vl);
    vuint8m1_t r3 = __riscv_vrgather_vv_u8m1(x, i3, vl);
    vuint8m1_t sum = __riscv_vxor_vv_u8m1(
        x, __riscv_vxor_vv_u8m1(r1, __riscv_vxor_vv_u8m1(r2, r3, vl), vl), vl);
    vuint8m1_t t8 = rvv_xtime(rvv_xtime(rvv_xtime(sum, vl), vl), vl);
    vuint8m1_t t4 = rvv_xtime(
        rvv_xtime(__riscv_vxor_vv_u8m1(x, r2, vl), vl), vl);
    vuint8m1_t t2 = rvv_xtime(__riscv_vxor_vv_u8m1(x, r1, vl), vl);

    return __riscv_vxor_vv_u8m1(
        t8, __riscv_vxor_vv_u8m1(t4, __riscv_vxor_vv_u8m1(
            t2, __riscv_vxor_vv_u8m1(sum, x, vl), vl), vl), vl);
}

static inline __attribute__((always_inline))
vuint8m1_t rvv_load_table(const uint8_t *p)
{
    return __riscv_vle8_v_u8m1(p, RVV_TABLE_BYTES);
}

void aes_encrypt_block_rvv(uint8_t out[16], const uint8_t in[16],
                           const uint8_t *rk_b, int rounds)
{
    uint8_t buf[16] __attribute__((aligned(16)));
    const size_t vl = RVV_VL;
    vuint8m1_t t0, t1, t2, t3, t4, t5, t6, t7;
    vuint8m1_t shift, i1, i2, i3;
    vuint8m1_t s;
    int round;

    memcpy(buf, in, sizeof(buf));

    t0 = rvv_load_table(rvv_sbox + 0);
    t1 = rvv_load_table(rvv_sbox + 32);
    t2 = rvv_load_table(rvv_sbox + 64);
    t3 = rvv_load_table(rvv_sbox + 96);
    t4 = rvv_load_table(rvv_sbox + 128);
    t5 = rvv_load_table(rvv_sbox + 160);
    t6 = rvv_load_table(rvv_sbox + 192);
    t7 = rvv_load_table(rvv_sbox + 224);
    shift = __riscv_vle8_v_u8m1(rvv_shift_idx, vl);
    i1 = __riscv_vle8_v_u8m1(rvv_rot1_idx, vl);
    i2 = __riscv_vle8_v_u8m1(rvv_rot2_idx, vl);
    i3 = __riscv_vle8_v_u8m1(rvv_rot3_idx, vl);

    s = __riscv_vle8_v_u8m1(buf, vl);
    s = __riscv_vxor_vv_u8m1(s, __riscv_vle8_v_u8m1(rk_b, vl), vl);
    for (round = 1; round < rounds; round++) {
        s = rvv_sbox_lookup(s, vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s = __riscv_vrgather_vv_u8m1(s, shift, vl);
        s = rvv_mix(s, vl, i1, i2, i3);
        s = __riscv_vxor_vv_u8m1(
            s, __riscv_vle8_v_u8m1(rk_b + 16 * round, vl), vl);
    }
    s = rvv_sbox_lookup(s, vl, t0, t1, t2, t3, t4, t5, t6, t7);
    s = __riscv_vrgather_vv_u8m1(s, shift, vl);
    s = __riscv_vxor_vv_u8m1(
        s, __riscv_vle8_v_u8m1(rk_b + 16 * rounds, vl), vl);

    __riscv_vse8_v_u8m1(buf, s, vl);
    memcpy(out, buf, sizeof(buf));
}

void aes_decrypt_block_rvv(uint8_t out[16], const uint8_t in[16],
                           const uint8_t *rk_b, int rounds)
{
    uint8_t buf[16] __attribute__((aligned(16)));
    const size_t vl = RVV_VL;
    vuint8m1_t t0, t1, t2, t3, t4, t5, t6, t7;
    vuint8m1_t inv_shift, i1, i2, i3;
    vuint8m1_t s;
    int round;

    memcpy(buf, in, sizeof(buf));

    t0 = rvv_load_table(rvv_inv_sbox + 0);
    t1 = rvv_load_table(rvv_inv_sbox + 32);
    t2 = rvv_load_table(rvv_inv_sbox + 64);
    t3 = rvv_load_table(rvv_inv_sbox + 96);
    t4 = rvv_load_table(rvv_inv_sbox + 128);
    t5 = rvv_load_table(rvv_inv_sbox + 160);
    t6 = rvv_load_table(rvv_inv_sbox + 192);
    t7 = rvv_load_table(rvv_inv_sbox + 224);
    inv_shift = __riscv_vle8_v_u8m1(rvv_inv_shift_idx, vl);
    i1 = __riscv_vle8_v_u8m1(rvv_rot1_idx, vl);
    i2 = __riscv_vle8_v_u8m1(rvv_rot2_idx, vl);
    i3 = __riscv_vle8_v_u8m1(rvv_rot3_idx, vl);

    s = __riscv_vle8_v_u8m1(buf, vl);
    s = __riscv_vxor_vv_u8m1(
        s, __riscv_vle8_v_u8m1(rk_b + 16 * rounds, vl), vl);
    for (round = rounds - 1; round > 0; round--) {
        s = __riscv_vrgather_vv_u8m1(s, inv_shift, vl);
        s = rvv_sbox_lookup(s, vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s = __riscv_vxor_vv_u8m1(
            s, __riscv_vle8_v_u8m1(rk_b + 16 * round, vl), vl);
        s = rvv_inv_mix(s, vl, i1, i2, i3);
    }
    s = __riscv_vrgather_vv_u8m1(s, inv_shift, vl);
    s = rvv_sbox_lookup(s, vl, t0, t1, t2, t3, t4, t5, t6, t7);
    s = __riscv_vxor_vv_u8m1(s, __riscv_vle8_v_u8m1(rk_b, vl), vl);

    __riscv_vse8_v_u8m1(buf, s, vl);
    memcpy(out, buf, sizeof(buf));
}


void aes_encrypt_blocks_rvv(uint8_t *out, const uint8_t *in, size_t blocks,
                            const uint8_t *rk_b, int rounds)
{
    const size_t vl = RVV_VL;
    const size_t vl2 = 32;
    vuint8m1_t t0, t1, t2, t3, t4, t5, t6, t7;
    vuint8m1_t shift, i1, i2, i3;
    vuint8m1_t shift2, j1, j2, j3;
    uint8_t rkb2[15][32];
    size_t b = 0;
    int round;

    for (round = 0; round <= rounds; round++) {
        memcpy(rkb2[round], rk_b + 16 * round, 16);
        memcpy(rkb2[round] + 16, rk_b + 16 * round, 16);
    }
    t0 = rvv_load_table(rvv_sbox + 0);
    t1 = rvv_load_table(rvv_sbox + 32);
    t2 = rvv_load_table(rvv_sbox + 64);
    t3 = rvv_load_table(rvv_sbox + 96);
    t4 = rvv_load_table(rvv_sbox + 128);
    t5 = rvv_load_table(rvv_sbox + 160);
    t6 = rvv_load_table(rvv_sbox + 192);
    t7 = rvv_load_table(rvv_sbox + 224);
    shift2 = __riscv_vle8_v_u8m1(rvv_shift_idx2, vl2);
    j1 = __riscv_vle8_v_u8m1(rvv_rot1_idx2, vl2);
    j2 = __riscv_vle8_v_u8m1(rvv_rot2_idx2, vl2);
    j3 = __riscv_vle8_v_u8m1(rvv_rot3_idx2, vl2);

    for (; b + 8 <= blocks; b += 8) {
        vuint8m1_t s0 = __riscv_vle8_v_u8m1(in + 16 * b, vl2);
        vuint8m1_t s1 = __riscv_vle8_v_u8m1(in + 16 * b + 32, vl2);
        vuint8m1_t s2 = __riscv_vle8_v_u8m1(in + 16 * b + 64, vl2);
        vuint8m1_t s3 = __riscv_vle8_v_u8m1(in + 16 * b + 96, vl2);
        vuint8m1_t k = __riscv_vle8_v_u8m1(rkb2[0], vl2);

        s0 = __riscv_vxor_vv_u8m1(s0, k, vl2);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl2);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl2);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl2);
        for (round = 1; round < rounds; round++) {
            k = __riscv_vle8_v_u8m1(rkb2[round], vl2);
            s0 = rvv_sbox_lookup(s0, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s1 = rvv_sbox_lookup(s1, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s2 = rvv_sbox_lookup(s2, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s3 = rvv_sbox_lookup(s3, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s0 = __riscv_vrgather_vv_u8m1(s0, shift2, vl2);
            s1 = __riscv_vrgather_vv_u8m1(s1, shift2, vl2);
            s2 = __riscv_vrgather_vv_u8m1(s2, shift2, vl2);
            s3 = __riscv_vrgather_vv_u8m1(s3, shift2, vl2);
            s0 = rvv_mix(s0, vl2, j1, j2, j3);
            s1 = rvv_mix(s1, vl2, j1, j2, j3);
            s2 = rvv_mix(s2, vl2, j1, j2, j3);
            s3 = rvv_mix(s3, vl2, j1, j2, j3);
            s0 = __riscv_vxor_vv_u8m1(s0, k, vl2);
            s1 = __riscv_vxor_vv_u8m1(s1, k, vl2);
            s2 = __riscv_vxor_vv_u8m1(s2, k, vl2);
            s3 = __riscv_vxor_vv_u8m1(s3, k, vl2);
        }
        k = __riscv_vle8_v_u8m1(rkb2[rounds], vl2);
        s0 = rvv_sbox_lookup(s0, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s1 = rvv_sbox_lookup(s1, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s2 = rvv_sbox_lookup(s2, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s3 = rvv_sbox_lookup(s3, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s0 = __riscv_vrgather_vv_u8m1(s0, shift2, vl2);
        s1 = __riscv_vrgather_vv_u8m1(s1, shift2, vl2);
        s2 = __riscv_vrgather_vv_u8m1(s2, shift2, vl2);
        s3 = __riscv_vrgather_vv_u8m1(s3, shift2, vl2);
        s0 = __riscv_vxor_vv_u8m1(s0, k, vl2);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl2);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl2);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b, s0, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b + 32, s1, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b + 64, s2, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b + 96, s3, vl2);
    }
    shift = __riscv_vle8_v_u8m1(rvv_shift_idx, vl);
    i1 = __riscv_vle8_v_u8m1(rvv_rot1_idx, vl);
    i2 = __riscv_vle8_v_u8m1(rvv_rot2_idx, vl);
    i3 = __riscv_vle8_v_u8m1(rvv_rot3_idx, vl);

    for (; b + 4 <= blocks; b += 4) {
        vuint8m1_t s0 = __riscv_vle8_v_u8m1(in + 16 * b, vl);
        vuint8m1_t s1 = __riscv_vle8_v_u8m1(in + 16 * b + 16, vl);
        vuint8m1_t s2 = __riscv_vle8_v_u8m1(in + 16 * b + 32, vl);
        vuint8m1_t s3 = __riscv_vle8_v_u8m1(in + 16 * b + 48, vl);
        vuint8m1_t k = __riscv_vle8_v_u8m1(rk_b, vl);

        s0 = __riscv_vxor_vv_u8m1(s0, k, vl);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl);
        for (round = 1; round < rounds; round++) {
            k = __riscv_vle8_v_u8m1(rk_b + 16 * round, vl);
            s0 = rvv_sbox_lookup(s0, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s1 = rvv_sbox_lookup(s1, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s2 = rvv_sbox_lookup(s2, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s3 = rvv_sbox_lookup(s3, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s0 = __riscv_vrgather_vv_u8m1(s0, shift, vl);
            s1 = __riscv_vrgather_vv_u8m1(s1, shift, vl);
            s2 = __riscv_vrgather_vv_u8m1(s2, shift, vl);
            s3 = __riscv_vrgather_vv_u8m1(s3, shift, vl);
            s0 = rvv_mix(s0, vl, i1, i2, i3);
            s1 = rvv_mix(s1, vl, i1, i2, i3);
            s2 = rvv_mix(s2, vl, i1, i2, i3);
            s3 = rvv_mix(s3, vl, i1, i2, i3);
            s0 = __riscv_vxor_vv_u8m1(s0, k, vl);
            s1 = __riscv_vxor_vv_u8m1(s1, k, vl);
            s2 = __riscv_vxor_vv_u8m1(s2, k, vl);
            s3 = __riscv_vxor_vv_u8m1(s3, k, vl);
        }
        k = __riscv_vle8_v_u8m1(rk_b + 16 * rounds, vl);
        s0 = rvv_sbox_lookup(s0, vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s1 = rvv_sbox_lookup(s1, vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s2 = rvv_sbox_lookup(s2, vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s3 = rvv_sbox_lookup(s3, vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s0 = __riscv_vrgather_vv_u8m1(s0, shift, vl);
        s1 = __riscv_vrgather_vv_u8m1(s1, shift, vl);
        s2 = __riscv_vrgather_vv_u8m1(s2, shift, vl);
        s3 = __riscv_vrgather_vv_u8m1(s3, shift, vl);
        s0 = __riscv_vxor_vv_u8m1(s0, k, vl);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl);
        __riscv_vse8_v_u8m1(out + 16 * b, s0, vl);
        __riscv_vse8_v_u8m1(out + 16 * b + 16, s1, vl);
        __riscv_vse8_v_u8m1(out + 16 * b + 32, s2, vl);
        __riscv_vse8_v_u8m1(out + 16 * b + 48, s3, vl);
    }
}

void aes_decrypt_blocks_rvv(uint8_t *out, const uint8_t *in, size_t blocks,
                            const uint8_t *rk_b, int rounds)
{
    const size_t vl = RVV_VL;
    const size_t vl2 = 32;
    vuint8m1_t t0, t1, t2, t3, t4, t5, t6, t7;
    vuint8m1_t inv_shift, i1, i2, i3;
    vuint8m1_t inv_shift2, j1, j2, j3;
    uint8_t rkb2[15][32];
    size_t b = 0;
    int round;

    for (round = 0; round <= rounds; round++) {
        memcpy(rkb2[round], rk_b + 16 * round, 16);
        memcpy(rkb2[round] + 16, rk_b + 16 * round, 16);
    }
    t0 = rvv_load_table(rvv_inv_sbox + 0);
    t1 = rvv_load_table(rvv_inv_sbox + 32);
    t2 = rvv_load_table(rvv_inv_sbox + 64);
    t3 = rvv_load_table(rvv_inv_sbox + 96);
    t4 = rvv_load_table(rvv_inv_sbox + 128);
    t5 = rvv_load_table(rvv_inv_sbox + 160);
    t6 = rvv_load_table(rvv_inv_sbox + 192);
    t7 = rvv_load_table(rvv_inv_sbox + 224);
    inv_shift2 = __riscv_vle8_v_u8m1(rvv_inv_shift_idx2, vl2);
    j1 = __riscv_vle8_v_u8m1(rvv_rot1_idx2, vl2);
    j2 = __riscv_vle8_v_u8m1(rvv_rot2_idx2, vl2);
    j3 = __riscv_vle8_v_u8m1(rvv_rot3_idx2, vl2);

    for (; b + 8 <= blocks; b += 8) {
        vuint8m1_t s0 = __riscv_vle8_v_u8m1(in + 16 * b, vl2);
        vuint8m1_t s1 = __riscv_vle8_v_u8m1(in + 16 * b + 32, vl2);
        vuint8m1_t s2 = __riscv_vle8_v_u8m1(in + 16 * b + 64, vl2);
        vuint8m1_t s3 = __riscv_vle8_v_u8m1(in + 16 * b + 96, vl2);
        vuint8m1_t k = __riscv_vle8_v_u8m1(rkb2[rounds], vl2);

        s0 = __riscv_vxor_vv_u8m1(s0, k, vl2);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl2);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl2);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl2);
        for (round = rounds - 1; round > 0; round--) {
            k = __riscv_vle8_v_u8m1(rkb2[round], vl2);
            s0 = __riscv_vrgather_vv_u8m1(s0, inv_shift2, vl2);
            s1 = __riscv_vrgather_vv_u8m1(s1, inv_shift2, vl2);
            s2 = __riscv_vrgather_vv_u8m1(s2, inv_shift2, vl2);
            s3 = __riscv_vrgather_vv_u8m1(s3, inv_shift2, vl2);
            s0 = rvv_sbox_lookup(s0, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s1 = rvv_sbox_lookup(s1, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s2 = rvv_sbox_lookup(s2, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s3 = rvv_sbox_lookup(s3, vl2, t0, t1, t2, t3, t4, t5, t6, t7);
            s0 = __riscv_vxor_vv_u8m1(s0, k, vl2);
            s1 = __riscv_vxor_vv_u8m1(s1, k, vl2);
            s2 = __riscv_vxor_vv_u8m1(s2, k, vl2);
            s3 = __riscv_vxor_vv_u8m1(s3, k, vl2);
            s0 = rvv_inv_mix(s0, vl2, j1, j2, j3);
            s1 = rvv_inv_mix(s1, vl2, j1, j2, j3);
            s2 = rvv_inv_mix(s2, vl2, j1, j2, j3);
            s3 = rvv_inv_mix(s3, vl2, j1, j2, j3);
        }
        k = __riscv_vle8_v_u8m1(rkb2[0], vl2);
        s0 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s0, inv_shift2, vl2),
                             vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s1 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s1, inv_shift2, vl2),
                             vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s2 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s2, inv_shift2, vl2),
                             vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s3 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s3, inv_shift2, vl2),
                             vl2, t0, t1, t2, t3, t4, t5, t6, t7);
        s0 = __riscv_vxor_vv_u8m1(s0, k, vl2);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl2);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl2);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b, s0, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b + 32, s1, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b + 64, s2, vl2);
        __riscv_vse8_v_u8m1(out + 16 * b + 96, s3, vl2);
    }
    inv_shift = __riscv_vle8_v_u8m1(rvv_inv_shift_idx, vl);
    i1 = __riscv_vle8_v_u8m1(rvv_rot1_idx, vl);
    i2 = __riscv_vle8_v_u8m1(rvv_rot2_idx, vl);
    i3 = __riscv_vle8_v_u8m1(rvv_rot3_idx, vl);

    for (; b + 4 <= blocks; b += 4) {
        vuint8m1_t s0 = __riscv_vle8_v_u8m1(in + 16 * b, vl);
        vuint8m1_t s1 = __riscv_vle8_v_u8m1(in + 16 * b + 16, vl);
        vuint8m1_t s2 = __riscv_vle8_v_u8m1(in + 16 * b + 32, vl);
        vuint8m1_t s3 = __riscv_vle8_v_u8m1(in + 16 * b + 48, vl);
        vuint8m1_t k = __riscv_vle8_v_u8m1(rk_b + 16 * rounds, vl);

        s0 = __riscv_vxor_vv_u8m1(s0, k, vl);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl);
        for (round = rounds - 1; round > 0; round--) {
            k = __riscv_vle8_v_u8m1(rk_b + 16 * round, vl);
            s0 = __riscv_vrgather_vv_u8m1(s0, inv_shift, vl);
            s1 = __riscv_vrgather_vv_u8m1(s1, inv_shift, vl);
            s2 = __riscv_vrgather_vv_u8m1(s2, inv_shift, vl);
            s3 = __riscv_vrgather_vv_u8m1(s3, inv_shift, vl);
            s0 = rvv_sbox_lookup(s0, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s1 = rvv_sbox_lookup(s1, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s2 = rvv_sbox_lookup(s2, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s3 = rvv_sbox_lookup(s3, vl, t0, t1, t2, t3, t4, t5, t6, t7);
            s0 = __riscv_vxor_vv_u8m1(s0, k, vl);
            s1 = __riscv_vxor_vv_u8m1(s1, k, vl);
            s2 = __riscv_vxor_vv_u8m1(s2, k, vl);
            s3 = __riscv_vxor_vv_u8m1(s3, k, vl);
            s0 = rvv_inv_mix(s0, vl, i1, i2, i3);
            s1 = rvv_inv_mix(s1, vl, i1, i2, i3);
            s2 = rvv_inv_mix(s2, vl, i1, i2, i3);
            s3 = rvv_inv_mix(s3, vl, i1, i2, i3);
        }
        k = __riscv_vle8_v_u8m1(rk_b, vl);
        s0 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s0, inv_shift, vl),
                             vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s1 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s1, inv_shift, vl),
                             vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s2 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s2, inv_shift, vl),
                             vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s3 = rvv_sbox_lookup(__riscv_vrgather_vv_u8m1(s3, inv_shift, vl),
                             vl, t0, t1, t2, t3, t4, t5, t6, t7);
        s0 = __riscv_vxor_vv_u8m1(s0, k, vl);
        s1 = __riscv_vxor_vv_u8m1(s1, k, vl);
        s2 = __riscv_vxor_vv_u8m1(s2, k, vl);
        s3 = __riscv_vxor_vv_u8m1(s3, k, vl);
        __riscv_vse8_v_u8m1(out + 16 * b, s0, vl);
        __riscv_vse8_v_u8m1(out + 16 * b + 16, s1, vl);
        __riscv_vse8_v_u8m1(out + 16 * b + 32, s2, vl);
        __riscv_vse8_v_u8m1(out + 16 * b + 48, s3, vl);
    }
}

#endif /* __riscv && AES_USE_RVV */
