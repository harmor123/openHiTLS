// ChaCha20 nine-block RVV 1.0 + scalar hybrid kernel (paired issue).
//
// void chacha20_rvv_blocks9(uint8_t *out, const uint8_t *in, size_t groups,
//                           const uint32_t key[8], uint32_t counter,
//                           const uint8_t nonce[12]);
//
// Eight vector blocks (lanes 0..7) plus one scalar block.  Every vector
// instruction is paired with the corresponding scalar instruction so the
// integer and vector pipelines can overlap.
//
// Structure follows the same vector+scalar pairing used by OpenSSL's
// chacha-riscv64-v-zbb implementation.

    .text
    .p2align 2

// One quarter-round sub-step for four independent QRs.
//   A += B; C ^= A; C = rotl(C, shift)
.macro SUBROUND va0, va1, va2, va3, vb0, vb1, vb2, vb3, vc0, vc1, vc2, vc3, sa0, sa1, sa2, sa3, sb0, sb1, sb2, sb3, sc0, sc1, sc2, sc3, shift, vt0, vt1, vt2, vt3
    vadd.vv \va0, \va0, \vb0
    add  \sa0, \sa0, \sb0
    vadd.vv \va1, \va1, \vb1
    add  \sa1, \sa1, \sb1
    vadd.vv \va2, \va2, \vb2
    add  \sa2, \sa2, \sb2
    vadd.vv \va3, \va3, \vb3
    add  \sa3, \sa3, \sb3

    vxor.vv \vc0, \vc0, \va0
    xor  \sc0, \sc0, \sa0
    vxor.vv \vc1, \vc1, \va1
    xor  \sc1, \sc1, \sa1
    vxor.vv \vc2, \vc2, \va2
    xor  \sc2, \sc2, \sa2
    vxor.vv \vc3, \vc3, \va3
    xor  \sc3, \sc3, \sa3

    vsll.vi \vt0, \vc0, \shift
    vsll.vi \vt1, \vc1, \shift
    vsll.vi \vt2, \vc2, \shift
    vsll.vi \vt3, \vc3, \shift
    vsrl.vi \vc0, \vc0, 32-\shift
    vsrl.vi \vc1, \vc1, 32-\shift
    vsrl.vi \vc2, \vc2, 32-\shift
    vsrl.vi \vc3, \vc3, 32-\shift
    vor.vv  \vc0, \vc0, \vt0
    roriw \sc0, \sc0, 32-\shift
    vor.vv  \vc1, \vc1, \vt1
    roriw \sc1, \sc1, 32-\shift
    vor.vv  \vc2, \vc2, \vt2
    roriw \sc2, \sc2, 32-\shift
    vor.vv  \vc3, \vc3, \vt3
    roriw \sc3, \sc3, 32-\shift
.endm

// Four QRs using x0,x4,x8,x12 / x1,x5,x9,x13 / x2,x6,x10,x14 / x3,x7,x11,x15.
.macro QRGROUP_COLUMN
    SUBROUND v0, v1, v2, v3, v4, v5, v6, v7, v12, v13, v14, v15, a5, a6, a7, s0, s1, s2, s3, s4, s9, s10, s11, t0, 16, v24, v25, v26, v27
    SUBROUND v8, v9, v10, v11, v12, v13, v14, v15, v4, v5, v6, v7, s5, s6, s7, s8, s9, s10, s11, t0, s1, s2, s3, s4, 12, v24, v25, v26, v27
    SUBROUND v0, v1, v2, v3, v4, v5, v6, v7, v12, v13, v14, v15, a5, a6, a7, s0, s1, s2, s3, s4, s9, s10, s11, t0, 8, v24, v25, v26, v27
    SUBROUND v8, v9, v10, v11, v12, v13, v14, v15, v4, v5, v6, v7, s5, s6, s7, s8, s9, s10, s11, t0, s1, s2, s3, s4, 7, v24, v25, v26, v27
.endm

// Diagonal QRs:
//   (x0,x5,x10,x15) (x1,x6,x11,x12) (x2,x7,x8,x13) (x3,x4,x9,x14)
.macro QRGROUP_DIAG
    SUBROUND v0, v1, v2, v3, v5, v6, v7, v4, v15, v12, v13, v14, a5, a6, a7, s0, s2, s3, s4, s1, t0, s9, s10, s11, 16, v24, v25, v26, v27
    SUBROUND v10, v11, v8, v9, v15, v12, v13, v14, v5, v6, v7, v4, s7, s8, s5, s6, t0, s9, s10, s11, s2, s3, s4, s1, 12, v24, v25, v26, v27
    SUBROUND v0, v1, v2, v3, v5, v6, v7, v4, v15, v12, v13, v14, a5, a6, a7, s0, s2, s3, s4, s1, t0, s9, s10, s11, 8, v24, v25, v26, v27
    SUBROUND v10, v11, v8, v9, v15, v12, v13, v14, v5, v6, v7, v4, s7, s8, s5, s6, t0, s9, s10, s11, s2, s3, s4, s1, 7, v24, v25, v26, v27
.endm

    .globl chacha20_rvv_blocks9
    .type chacha20_rvv_blocks9, @function
    .p2align 4
chacha20_rvv_blocks9:
    beqz a2, .Ldone9

    addi sp, sp, -736
    sd   s0, 576(sp)
    sd   s1, 584(sp)
    sd   s2, 592(sp)
    sd   s3, 600(sp)
    sd   s4, 608(sp)
    sd   s5, 616(sp)
    sd   s6, 624(sp)
    sd   s7, 632(sp)
    sd   s8, 640(sp)
    sd   s9, 648(sp)
    sd   s10, 656(sp)
    sd   s11, 664(sp)

    vsetivli t1, 8, e32, m1, ta, ma

    // Vector initial state: constants, key, counter placeholder, nonce.
    li   t1, 0x61707865
    vmv.v.x v0, t1
    li   t1, 0x3320646e
    vmv.v.x v1, t1
    li   t1, 0x79622d32
    vmv.v.x v2, t1
    li   t1, 0x6b206574
    vmv.v.x v3, t1
    lw   t1, 0(a3)
    vmv.v.x v4, t1
    lw   t1, 4(a3)
    vmv.v.x v5, t1
    lw   t1, 8(a3)
    vmv.v.x v6, t1
    lw   t1, 12(a3)
    vmv.v.x v7, t1
    lw   t1, 16(a3)
    vmv.v.x v8, t1
    lw   t1, 20(a3)
    vmv.v.x v9, t1
    lw   t1, 24(a3)
    vmv.v.x v10, t1
    lw   t1, 28(a3)
    vmv.v.x v11, t1
    vmv.v.i v12, 0
    lw   t1, 0(a5)
    vmv.v.x v13, t1
    lw   t1, 4(a5)
    vmv.v.x v14, t1
    lw   t1, 8(a5)
    vmv.v.x v15, t1

    mv   t5, sp
    vse32.v v0, (t5)
    addi t5, t5, 32
    vse32.v v1, (t5)
    addi t5, t5, 32
    vse32.v v2, (t5)
    addi t5, t5, 32
    vse32.v v3, (t5)
    addi t5, t5, 32
    vse32.v v4, (t5)
    addi t5, t5, 32
    vse32.v v5, (t5)
    addi t5, t5, 32
    vse32.v v6, (t5)
    addi t5, t5, 32
    vse32.v v7, (t5)
    addi t5, t5, 32
    vse32.v v8, (t5)
    addi t5, t5, 32
    vse32.v v9, (t5)
    addi t5, t5, 32
    vse32.v v10, (t5)
    addi t5, t5, 32
    vse32.v v11, (t5)
    addi t5, t5, 32
    vse32.v v12, (t5)
    addi t5, t5, 32
    vse32.v v13, (t5)
    addi t5, t5, 32
    vse32.v v14, (t5)
    addi t5, t5, 32
    vse32.v v15, (t5)

    // Scalar initial state from lane 0 of the vector state; x12 = counter + 8.
    lw   a5, 0(sp)
    lw   a6, 32(sp)
    lw   a7, 64(sp)
    lw   s0, 96(sp)
    lw   s1, 128(sp)
    lw   s2, 160(sp)
    lw   s3, 192(sp)
    lw   s4, 224(sp)
    lw   s5, 256(sp)
    lw   s6, 288(sp)
    lw   s7, 320(sp)
    lw   s8, 352(sp)
    addi s9, a4, 8
    lw   s10, 416(sp)
    lw   s11, 448(sp)
    lw   t0, 480(sp)
    sw   a5, 512(sp)
    sw   a6, 516(sp)
    sw   a7, 520(sp)
    sw   s0, 524(sp)
    sw   s1, 528(sp)
    sw   s2, 532(sp)
    sw   s3, 536(sp)
    sw   s4, 540(sp)
    sw   s5, 544(sp)
    sw   s6, 548(sp)
    sw   s7, 552(sp)
    sw   s8, 556(sp)
    sw   s9, 560(sp)
    sw   s10, 564(sp)
    sw   s11, 568(sp)
    sw   t0, 572(sp)

    mv   a3, a4

.Lgroup9:
    vid.v v12
    vadd.vx v12, v12, a3

    // Rebuild the vector state from the saved scalar words.
    lw   t1, 0(sp)
    vmv.v.x v0, t1
    lw   t1, 32(sp)
    vmv.v.x v1, t1
    lw   t1, 64(sp)
    vmv.v.x v2, t1
    lw   t1, 96(sp)
    vmv.v.x v3, t1
    lw   t1, 128(sp)
    vmv.v.x v4, t1
    lw   t1, 160(sp)
    vmv.v.x v5, t1
    lw   t1, 192(sp)
    vmv.v.x v6, t1
    lw   t1, 224(sp)
    vmv.v.x v7, t1
    lw   t1, 256(sp)
    vmv.v.x v8, t1
    lw   t1, 288(sp)
    vmv.v.x v9, t1
    lw   t1, 320(sp)
    vmv.v.x v10, t1
    lw   t1, 352(sp)
    vmv.v.x v11, t1
    lw   t1, 416(sp)
    vmv.v.x v13, t1
    lw   t1, 448(sp)
    vmv.v.x v14, t1
    lw   t1, 480(sp)
    vmv.v.x v15, t1

    addi t1, a3, 8
    sw   t1, 560(sp)
    lw   a5, 512(sp)
    lw   a6, 516(sp)
    lw   a7, 520(sp)
    lw   s0, 524(sp)
    lw   s1, 528(sp)
    lw   s2, 532(sp)
    lw   s3, 536(sp)
    lw   s4, 540(sp)
    lw   s5, 544(sp)
    lw   s6, 548(sp)
    lw   s7, 552(sp)
    lw   s8, 556(sp)
    lw   s9, 560(sp)
    lw   s10, 564(sp)
    lw   s11, 568(sp)
    lw   t0, 572(sp)

    // Preload the top-half input while the round loop runs.
    li   t4, 64
    vlsseg8e32.v v16, (a1), t4

    li   a4, 10
    .p2align 4
.Lround9:
    QRGROUP_COLUMN
    QRGROUP_DIAG
    addi a4, a4, -1
    bnez a4, .Lround9

    // Vector feed-forward from the saved scalar words.
    lw   t1, 0(sp)
    vadd.vx v0, v0, t1
    lw   t1, 32(sp)
    vadd.vx v1, v1, t1
    lw   t1, 64(sp)
    vadd.vx v2, v2, t1
    lw   t1, 96(sp)
    vadd.vx v3, v3, t1
    lw   t1, 128(sp)
    vadd.vx v4, v4, t1
    lw   t1, 160(sp)
    vadd.vx v5, v5, t1
    lw   t1, 192(sp)
    vadd.vx v6, v6, t1
    lw   t1, 224(sp)
    vadd.vx v7, v7, t1
    lw   t1, 256(sp)
    vadd.vx v8, v8, t1
    lw   t1, 288(sp)
    vadd.vx v9, v9, t1
    lw   t1, 320(sp)
    vadd.vx v10, v10, t1
    lw   t1, 352(sp)
    vadd.vx v11, v11, t1
    vid.v v24
    vadd.vx v24, v24, a3
    vadd.vv v12, v12, v24
    lw   t1, 416(sp)
    vadd.vx v13, v13, t1
    lw   t1, 448(sp)
    vadd.vx v14, v14, t1
    lw   t1, 480(sp)
    vadd.vx v15, v15, t1

    addi t1, a1, 32
    vlsseg8e32.v v24, (t1), t4

    // Scalar feed-forward of block 8, overlapping the bottom-half load.
    lw   t3, 512(sp)
    add  a5, a5, t3
    lw   t3, 516(sp)
    add  a6, a6, t3
    lw   t3, 520(sp)
    add  a7, a7, t3
    lw   t3, 524(sp)
    add  s0, s0, t3
    lw   t3, 528(sp)
    add  s1, s1, t3
    lw   t3, 532(sp)
    add  s2, s2, t3
    lw   t3, 536(sp)
    add  s3, s3, t3
    lw   t3, 540(sp)
    add  s4, s4, t3
    lw   t3, 544(sp)
    add  s5, s5, t3
    lw   t3, 548(sp)
    add  s6, s6, t3
    lw   t3, 552(sp)
    add  s7, s7, t3
    lw   t3, 556(sp)
    add  s8, s8, t3
    lw   t3, 560(sp)
    add  s9, s9, t3
    lw   t3, 564(sp)
    add  s10, s10, t3
    lw   t3, 568(sp)
    add  s11, s11, t3
    lw   t3, 572(sp)
    add  t0, t0, t3

    // Interleave the scalar state stores with the vector XOR.
    vxor.vv v16, v16, v0
    sw   a5, 672(sp)
    vxor.vv v17, v17, v1
    sw   a6, 676(sp)
    vxor.vv v18, v18, v2
    sw   a7, 680(sp)
    vxor.vv v19, v19, v3
    sw   s0, 684(sp)
    vxor.vv v20, v20, v4
    sw   s1, 688(sp)
    vxor.vv v21, v21, v5
    sw   s2, 692(sp)
    vxor.vv v22, v22, v6
    sw   s3, 696(sp)
    vxor.vv v23, v23, v7
    sw   s4, 700(sp)
    vxor.vv v24, v24, v8
    sw   s5, 704(sp)
    vxor.vv v25, v25, v9
    sw   s6, 708(sp)
    vxor.vv v26, v26, v10
    sw   s7, 712(sp)
    vxor.vv v27, v27, v11
    sw   s8, 716(sp)
    vxor.vv v28, v28, v12
    sw   s9, 720(sp)
    vxor.vv v29, v29, v13
    sw   s10, 724(sp)
    vxor.vv v30, v30, v14
    sw   s11, 728(sp)
    vxor.vv v31, v31, v15
    sw   t0, 732(sp)

    vssseg8e32.v v16, (a0), t4
    addi t2, a0, 32
    vssseg8e32.v v24, (t2), t4

    addi t1, a1, 512
    addi t2, sp, 672
    addi t3, a0, 512
    li   t5, 16
    vsetvli t5, t5, e32, m2, ta, ma
    vle32.v v8, (t1)
    vle32.v v16, (t2)
    vxor.vv v8, v8, v16
    vse32.v v8, (t3)
    vsetivli t5, 8, e32, m1, ta, ma

    addi a0, a0, 576
    addi a1, a1, 576
    addi a3, a3, 9
    addi a2, a2, -1
    bnez a2, .Lgroup9

    ld   s0, 576(sp)
    ld   s1, 584(sp)
    ld   s2, 592(sp)
    ld   s3, 600(sp)
    ld   s4, 608(sp)
    ld   s5, 616(sp)
    ld   s6, 624(sp)
    ld   s7, 632(sp)
    ld   s8, 640(sp)
    ld   s9, 648(sp)
    ld   s10, 656(sp)
    ld   s11, 664(sp)
    addi sp, sp, 736

.Ldone9:
    ret
