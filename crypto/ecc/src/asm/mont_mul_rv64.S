/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */
// Generic multi-precision Montgomery multiplication for RV64.
//
// void hitls_mont_mul_rv64(uint64_t *r, const uint64_t *a,
//                          const uint64_t *b, const uint64_t *mod,
//                          uint64_t n0, uint64_t limbs);
//
// CIOS: r = a * b * R^-1 mod mod, with R = 2^(64*limbs).
// The loop counts are public (curve size); the final conditional subtraction
// is performed with a mask, not a secret-dependent branch.

    .text
    .p2align 3
    .globl hitls_mont_mul_rv64
    .type hitls_mont_mul_rv64, @function
hitls_mont_mul_rv64:
    addi    sp, sp, -352
    sd      s0, 256(sp)
    sd      s1, 264(sp)
    sd      s2, 272(sp)
    sd      s3, 280(sp)
    sd      s4, 288(sp)
    sd      s5, 296(sp)
    sd      s6, 304(sp)
    sd      s7, 312(sp)
    sd      s8, 320(sp)
    sd      s9, 328(sp)
    sd      s10, 336(sp)
    sd      s11, 344(sp)

    mv      s0, a0              // r
    mv      s1, a1              // a
    mv      s2, a2              // b
    mv      s3, a3              // modulus
    mv      s4, a4              // n0
    mv      s5, a5              // limbs

    // t[0..15] = 0 at sp[0..127]
    mv      t0, sp
    li      t1, 16
1:
    sd      zero, 0(t0)
    addi    t0, t0, 8
    addi    t1, t1, -1
    bnez    t1, 1b

    slli    t0, s5, 3
    add     s6, sp, t0          // end of t[0..limbs-1]
    add     s11, s1, t0         // end of a
    addi    s7, sp, 128         // d[0] scratch area

    li      s8, 0               // i = 0
outer_loop:
    li      s9, 0               // carry = 0
    slli    t0, s8, 3
    add     t0, s2, t0
    ld      a7, 0(t0)           // b[i]

    mv      t2, sp              // t[j]
    mv      t3, s1              // a[j]
mul_loop:
    ld      a6, 0(t3)
    mul     t4, a6, a7
    mulhu   t5, a6, a7
    ld      t6, 0(t2)
    add     a6, t4, t6
    sltu    t0, a6, t4
    add     a6, a6, s9
    sltu    t1, a6, s9
    add     t0, t0, t1
    add     t5, t5, t0
    sd      a6, 0(t2)
    mv      s9, t5
    addi    t2, t2, 8
    addi    t3, t3, 8
    bne     t2, s6, mul_loop

    ld      t0, 0(s6)
    add     t0, t0, s9
    sltu    t1, t0, s9
    sd      t0, 0(s6)
    sd      t1, 8(s6)

    // m = t[0] * n0 mod 2^64
    ld      t0, 0(sp)
    mul     s10, t0, s4

    li      s9, 0
    mv      t2, sp
    mv      t3, s3
red_loop:
    ld      a6, 0(t3)
    mul     t4, s10, a6
    mulhu   t5, s10, a6
    ld      t6, 0(t2)
    add     a6, t4, t6
    sltu    t0, a6, t4
    add     a6, a6, s9
    sltu    t1, a6, s9
    add     t0, t0, t1
    add     t5, t5, t0
    sd      a6, 0(t2)
    mv      s9, t5
    addi    t2, t2, 8
    addi    t3, t3, 8
    bne     t2, s6, red_loop

    ld      t0, 0(s6)
    add     t0, t0, s9
    sltu    t1, t0, s9
    sd      t0, 0(s6)
    ld      t2, 8(s6)
    add     t2, t2, t1
    sd      t2, 8(s6)

    // t >>= 64
    mv      t2, sp
    mv      t3, sp
    addi    t3, t3, 8
shift_loop:
    ld      t4, 0(t3)
    sd      t4, 0(t2)
    addi    t2, t2, 8
    addi    t3, t3, 8
    bne     t2, s6, shift_loop
    ld      t4, 8(s6)
    sd      t4, 0(s6)
    sd      zero, 8(s6)

    addi    s8, s8, 1
    bne     s8, s5, outer_loop

    // d = t - modulus, borrow in s9
    li      s9, 0
    mv      t2, sp
    mv      t3, s3
    mv      t4, s7
sub_loop:
    ld      t0, 0(t2)
    ld      t1, 0(t3)
    sub     t5, t0, t1
    sltu    a6, t0, t1
    mv      a7, t5
    sub     t5, t5, s9
    sltu    t0, a7, s9
    or      a6, a6, t0
    sd      t5, 0(t4)
    mv      s9, a6
    addi    t2, t2, 8
    addi    t3, t3, 8
    addi    t4, t4, 8
    bne     t2, s6, sub_loop

    // need = (t[limbs] != 0) || (borrow == 0)
    ld      t0, 0(s6)
    sltu    t0, zero, t0
    xori    t1, s9, 1
    or      t0, t0, t1
    sub     s10, zero, t0       // mask = -need

    mv      t2, sp
    mv      t3, s7
    mv      t4, s0
    not     a7, s10
select_loop:
    ld      t5, 0(t2)
    ld      a6, 0(t3)
    and     t5, t5, a7
    and     a6, a6, s10
    or      t5, t5, a6
    sd      t5, 0(t4)
    addi    t2, t2, 8
    addi    t3, t3, 8
    addi    t4, t4, 8
    bne     t2, s6, select_loop

    ld      s0, 256(sp)
    ld      s1, 264(sp)
    ld      s2, 272(sp)
    ld      s3, 280(sp)
    ld      s4, 288(sp)
    ld      s5, 296(sp)
    ld      s6, 304(sp)
    ld      s7, 312(sp)
    ld      s8, 320(sp)
    ld      s9, 328(sp)
    ld      s10, 336(sp)
    ld      s11, 344(sp)
    addi    sp, sp, 352
    ret
    .size hitls_mont_mul_rv64, .-hitls_mont_mul_rv64
