/*
 * Copyright (c) 2026 harmor.
 *
 * Four-block bitsliced AES-128/192/256 encryption (Boyar-Peralta S-box).
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 *
 * The 64-bit planes carry four 16-bit block lanes.  The S-box and
 * MixColumns circuits use only bitwise/rotate operations, so no secret
 * value reaches a memory address.
 */

#include "aes_ct.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define M16 0xffffu

static uint64_t transpose8(uint64_t x)
{
    uint64_t t;
    t = (x ^ (x >> 7)) & 0x00aa00aa00aa00aaULL;
    x ^= t ^ (t << 7);
    t = (x ^ (x >> 14)) & 0x0000cccc0000ccccULL;
    x ^= t ^ (t << 14);
    t = (x ^ (x >> 28)) & 0x00000000f0f0f0f0ULL;
    x ^= t ^ (t << 28);
    return x;
}

static void block_to_planes(uint16_t p[8], const uint8_t s[16])
{
    uint64_t x0, x1, t0, t1;
    int k;
    memcpy(&x0, s, 8);
    memcpy(&x1, s + 8, 8);
    t0 = transpose8(x0);
    t1 = transpose8(x1);
    for (k = 0; k < 8; k++) {
        p[k] = (uint16_t)(((t0 >> (8 * k)) & 0xffu) |
                          (((t1 >> (8 * k)) & 0xffu) << 8));
    }
}

static void planes_to_block(uint8_t s[16], const uint16_t p[8])
{
    uint64_t t0 = 0, t1 = 0, x0, x1;
    int k;
    for (k = 0; k < 8; k++) {
        t0 |= (uint64_t)(p[k] & 0xffu) << (8 * k);
        t1 |= (uint64_t)((p[k] >> 8) & 0xffu) << (8 * k);
    }
    x0 = transpose8(t0);
    x1 = transpose8(t1);
    memcpy(s, &x0, 8);
    memcpy(s + 8, &x1, 8);
}

static void load4(uint64_t p[8], const uint8_t in[64])
{
    int b, k;
    for (k = 0; k < 8; k++) p[k] = 0;
    for (b = 0; b < 4; b++) {
        uint16_t q[8];
        block_to_planes(q, in + 16 * b);
        for (k = 0; k < 8; k++) p[k] |= (uint64_t)q[k] << (16 * b);
    }
}

static void store4(uint8_t out[64], const uint64_t p[8])
{
    int b, k;
    for (b = 0; b < 4; b++) {
        uint16_t q[8];
        for (k = 0; k < 8; k++) q[k] = (uint16_t)((p[k] >> (16 * b)) & M16);
        planes_to_block(out + 16 * b, q);
    }
}

static inline __attribute__((always_inline)) void bp_sbox(uint64_t a[8])
{
    uint64_t x0=a[7],x1=a[6],x2=a[5],x3=a[4],x4=a[3],x5=a[2],x6=a[1],x7=a[0];
    uint64_t y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18,y19,y20,y21;
    uint64_t z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17;
    uint64_t t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19;
    uint64_t t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39;
    uint64_t t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59;
    uint64_t t60,t61,t62,t63,t64,t65,t66,t67,s0,s1,s2,s3,s4,s5,s6,s7;

    y14=x3^x5; y13=x0^x6; y9=x0^x3; y8=x0^x5; t0=x1^x2; y1=t0^x7;
    y4=y1^x3; y12=y13^y14; y2=y1^x0; y5=y1^x6; y3=y5^y8; t1=x4^y12;
    y15=t1^x5; y20=t1^x1; y6=y15^x7; y10=y15^t0; y11=y20^y9; y7=x7^y11;
    y17=y10^y11; y19=y10^y8; y16=t0^y11; y21=y13^y16; y18=x0^y16;
    t2=y12&y15; t3=y3&y6; t4=t3^t2; t5=y4&x7; t6=t5^t2; t7=y13&y16;
    t8=y5&y1; t9=t8^t7; t10=y2&y7; t11=t10^t7; t12=y9&y11; t13=y14&y17;
    t14=t13^t12; t15=y8&y10; t16=t15^t12; t17=t4^t14; t18=t6^t16;
    t19=t9^t14; t20=t11^t16; t21=t17^y20; t22=t18^y19; t23=t19^y21;
    t24=t20^y18; t25=t21^t22; t26=t21&t23; t27=t24^t26; t28=t25&t27;
    t29=t28^t22; t30=t23^t24; t31=t22^t26; t32=t31&t30; t33=t32^t24;
    t34=t23^t33; t35=t27^t33; t36=t24&t35; t37=t36^t34; t38=t27^t36;
    t39=t29&t38; t40=t25^t39; t41=t40^t37; t42=t29^t33; t43=t29^t40;
    t44=t33^t37; t45=t42^t41; z0=t44&y15; z1=t37&y6; z2=t33&x7;
    z3=t43&y16; z4=t40&y1; z5=t29&y7; z6=t42&y11; z7=t45&y17;
    z8=t41&y10; z9=t44&y12; z10=t37&y3; z11=t33&y4; z12=t43&y13;
    z13=t40&y5; z14=t29&y2; z15=t42&y9; z16=t45&y14; z17=t41&y8;
    t46=z15^z16; t47=z10^z11; t48=z5^z13; t49=z9^z10; t50=z2^z12;
    t51=z2^z5; t52=z7^z8; t53=z0^z3; t54=z6^z7; t55=z16^z17;
    t56=z12^t48; t57=t50^t53; t58=z4^t46; t59=z3^t54; t60=t46^t57;
    t61=z14^t57; t62=t52^t58; t63=t49^t58; t64=z4^t59; t65=t61^t62;
    t66=z1^t63; s0=t59^t63; s6=t56^(~t62); s7=t48^(~t60); t67=t64^t65;
    s3=t53^t66; s4=t51^t66; s5=t47^t65; s1=t64^(~s3); s2=t55^(~t67);

    a[7]=s0; a[6]=s1; a[5]=s2; a[4]=s3;
    a[3]=s4; a[2]=s5; a[1]=s6; a[0]=s7;
}

static uint64_t shift_rows_word(uint64_t x)
{
    uint64_t y = x & ~(0x2222222222222222ULL | 0x4444444444444444ULL |
                       0x8888888888888888ULL);
    y |= (((x & 0x2220222022202220ULL) >> 4) |
          ((x & 0x0002000200020002ULL) << 12)) & 0x2222222222222222ULL;
    y |= (((x & 0x4400440044004400ULL) >> 8) |
          ((x & 0x0044004400440044ULL) << 8)) & 0x4444444444444444ULL;
    y |= (((x & 0x8000800080008000ULL) >> 12) |
          ((x & 0x0888088808880888ULL) << 4)) & 0x8888888888888888ULL;
    return y;
}

#define ROT1(x) ((((x) >> 1) & 0x7777777777777777ULL) | \
                 (((x) & 0x1111111111111111ULL) << 3))
#define ROT2(x) ((((x) >> 2) & 0x3333333333333333ULL) | \
                 (((x) & 0x3333333333333333ULL) << 2))
#define ROT3(x) ((((x) >> 3) & 0x1111111111111111ULL) | \
                 (((x) & 0x7777777777777777ULL) << 1))

static void mix_columns(uint64_t p[8])
{
    uint64_t a0=p[0],a1=p[1],a2=p[2],a3=p[3],a4=p[4],a5=p[5],a6=p[6],a7=p[7];
    uint64_t t0=a7, t1=a0^a7, t2=a1, t3=a2^a7;
    uint64_t t4=a3^a7, t5=a4, t6=a5, t7=a6;
    p[0]=t0^ROT1(t0)^ROT1(a0)^ROT2(a0)^ROT3(a0);
    p[1]=t1^ROT1(t1)^ROT1(a1)^ROT2(a1)^ROT3(a1);
    p[2]=t2^ROT1(t2)^ROT1(a2)^ROT2(a2)^ROT3(a2);
    p[3]=t3^ROT1(t3)^ROT1(a3)^ROT2(a3)^ROT3(a3);
    p[4]=t4^ROT1(t4)^ROT1(a4)^ROT2(a4)^ROT3(a4);
    p[5]=t5^ROT1(t5)^ROT1(a5)^ROT2(a5)^ROT3(a5);
    p[6]=t6^ROT1(t6)^ROT1(a6)^ROT2(a6)^ROT3(a6);
    p[7]=t7^ROT1(t7)^ROT1(a7)^ROT2(a7)^ROT3(a7);
}

static void enc4(uint8_t out[64], const uint8_t in[64],
                 const uint64_t rk[15][8], int rounds)
{
    uint64_t p[8];
    int r, k;
    load4(p, in);
    for (k = 0; k < 8; k++) p[k] ^= rk[0][k];
    for (r = 1; r < rounds; r++) {
        bp_sbox(p);
        for (k = 0; k < 8; k++) p[k] = shift_rows_word(p[k]);
        mix_columns(p);
        for (k = 0; k < 8; k++) p[k] ^= rk[r][k];
    }
    bp_sbox(p);
    for (k = 0; k < 8; k++) p[k] = shift_rows_word(p[k]);
    for (k = 0; k < 8; k++) p[k] ^= rk[rounds][k];
    store4(out, p);
}

void aes_bs4_encrypt_blocks(uint8_t *out, const uint8_t *in, size_t blocks,
                            const uint8_t *rk_b, int rounds)
{
    uint64_t rk[15][8];
    uint16_t q[8];
    size_t b;
    int i;
    int j;

    for (i = 0; i <= rounds; i++) {
        block_to_planes(q, rk_b + 16 * i);
        for (j = 0; j < 8; j++) {
            rk[i][j] = (uint64_t)q[j] * 0x0001000100010001ULL;
        }
    }
    for (b = 0; b + 4 <= blocks; b += 4) {
        enc4(out + 16 * b, in + 16 * b, rk, rounds);
    }
}
