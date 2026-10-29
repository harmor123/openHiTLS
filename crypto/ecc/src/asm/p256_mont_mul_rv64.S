// P-256 Montgomery multiplication, 4x64-bit limbs.
//
// void p256_mont_mul_p_rv64(uint64_t r[4], const uint64_t a[4],
//                           const uint64_t b[4]);
//
// Computes a*b*R^-1 mod p with p = 2^256 - 2^224 + 2^192 + 2^96 - 1 and
// R = 2^256.  The least-significant limb of p is -1, so -p^-1 mod 2^64 = 1
// and the Montgomery multiplier is simply the current low limb.

    .text
    .p2align 3

// acc += a3 * b[j] for one 64x64 product; t is kept on the stack.
.macro MULADD toff, boff
    ld    a5, \toff(sp)
    ld    a6, \boff(s1)
    mul   a7, a3, a6
    mulhu t2, a3, a6
    add   a7, a7, a5
    sltu  t3, a7, a5
    add   a7, a7, a4
    sltu  t4, a7, a4
    add   t2, t2, t3
    add   t2, t2, t4
    sd    a7, \toff(sp)
    mv    a4, t2
.endm

// Same as MULADD, but the multiplier limb comes from a register.
.macro MULADD_P toff, preg
    ld    a5, \toff(sp)
    mul   a7, a3, \preg
    mulhu t2, a3, \preg
    add   a7, a7, a5
    sltu  t3, a7, a5
    add   a7, a7, a4
    sltu  t4, a7, a4
    add   t2, t2, t3
    add   t2, t2, t4
    sd    a7, \toff(sp)
    mv    a4, t2
.endm

.macro PROD_ROW i
    ld    a3, (\i * 8)(s0)
    li    a4, 0
    MULADD (\i + 0) * 8, 0
    MULADD (\i + 1) * 8, 8
    MULADD (\i + 2) * 8, 16
    MULADD (\i + 3) * 8, 24
    sd    a4, (\i + 4) * 8(sp)
.endm

.macro RED_ROW i
    ld    a3, (\i * 8)(sp)
    li    a4, 0
    MULADD_P (\i + 0) * 8, s5
    MULADD_P (\i + 1) * 8, s6
    MULADD_P (\i + 2) * 8, s7
    MULADD_P (\i + 3) * 8, s8
    addi  t0, sp, (\i + 4) * 8
    /* Fixed-length carry propagation: the old loop branched on a carry
     * derived from secret operands.  t[] has 12 limbs, so row i needs
     * 8 - i further words. */
.rept (8 - \i)
    ld    a5, 0(t0)
    add   a5, a5, a4
    sd    a5, 0(t0)
    sltu  a4, a5, a4
    addi  t0, t0, 8
.endr
.endm

    .globl p256_mont_mul_p_rv64
    .type p256_mont_mul_p_rv64, @function
p256_mont_mul_p_rv64:
    addi sp, sp, -192
    sd   s0, 96(sp)
    sd   s1, 104(sp)
    sd   s2, 112(sp)
    sd   s3, 120(sp)
    sd   s4, 128(sp)
    sd   s5, 136(sp)
    sd   s6, 144(sp)
    sd   s7, 152(sp)
    sd   s8, 160(sp)
    sd   s9, 168(sp)

    mv   s0, a1
    mv   s1, a2
    mv   s2, a0
    li   s5, -1
    li   s6, 0xffffffff
    li   s7, 0
    li   s8, 0xffffffff00000001

    sd   zero, 0(sp)
    sd   zero, 8(sp)
    sd   zero, 16(sp)
    sd   zero, 24(sp)
    sd   zero, 32(sp)
    sd   zero, 40(sp)
    sd   zero, 48(sp)
    sd   zero, 56(sp)
    sd   zero, 64(sp)
    sd   zero, 72(sp)
    sd   zero, 80(sp)
    sd   zero, 88(sp)

    PROD_ROW 0
    PROD_ROW 1
    PROD_ROW 2
    PROD_ROW 3

    RED_ROW 0
    RED_ROW 1
    RED_ROW 2
    RED_ROW 3

    // result = t[4..7], high = t[8]; conditionally subtract p.
    ld   t0, 32(sp)
    ld   t1, 40(sp)
    ld   t2, 48(sp)
    ld   t3, 56(sp)
    ld   t4, 64(sp)

    sub  a3, t0, s5
    sltu a4, t0, s5             // borrow0

    sub  a5, t1, s6
    sltu a6, t1, s6             // raw borrow1
    sltu t5, a5, a4             // underflow subtracting borrow0
    sub  a5, a5, a4
    or   a4, a6, t5             // borrow1

    sub  a6, t2, s7
    sltu a7, t2, s7             // raw borrow2
    sltu t5, a6, a4
    sub  a6, a6, a4
    or   a4, a7, t5             // borrow2

    sub  a7, t3, s8
    sltu t5, t3, s8             // raw borrow3
    sltu t6, a7, a4
    sub  a7, a7, a4
    or   a4, t5, t6             // borrow3

    sltu t6, zero, t4           // high != 0
    xori t5, a4, 1              // no borrow
    or   t6, t6, t5
    neg  t6, t6                 // mask

    xor  a3, t0, a3
    and  a3, a3, t6
    xor  a3, t0, a3
    sd   a3, 0(s2)

    xor  a5, t1, a5
    and  a5, a5, t6
    xor  a5, t1, a5
    sd   a5, 8(s2)

    xor  a6, t2, a6
    and  a6, a6, t6
    xor  a6, t2, a6
    sd   a6, 16(s2)

    xor  a7, t3, a7
    and  a7, a7, t6
    xor  a7, t3, a7
    sd   a7, 24(s2)

    ld   s0, 96(sp)
    ld   s1, 104(sp)
    ld   s2, 112(sp)
    ld   s3, 120(sp)
    ld   s4, 128(sp)
    ld   s5, 136(sp)
    ld   s6, 144(sp)
    ld   s7, 152(sp)
    ld   s8, 160(sp)
    ld   s9, 168(sp)
    addi sp, sp, 192
    ret

    .size p256_mont_mul_p_rv64, . - p256_mont_mul_p_rv64

    .section .note.GNU-stack,"",@progbits
