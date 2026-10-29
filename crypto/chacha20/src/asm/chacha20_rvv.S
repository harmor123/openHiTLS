// ChaCha20 eight-block RVV 1.0 kernel.
//
// void chacha20_rvv_blocks(uint8_t *out, const uint8_t *in, size_t groups,
//                          const uint32_t key[8], uint32_t counter,
//                          const uint8_t nonce[12]);
//
// Each group is eight consecutive 64-byte blocks.  Lane b of every state
// vector belongs to block b.  VL is fixed at 8 (VLEN=256); no vsetvli is
// needed inside the round loop.  The initial state is kept on the stack so
// that the feed-forward add stays cheap.

    .text
    .p2align 2
    .globl chacha20_rvv_blocks
    .type chacha20_rvv_blocks, @function

// Four independent quarter rounds, interleaved to cover vector latency.
// v16-v19 are scratch registers for the three-instruction rotates.
.macro QR4 xa0, xb0, xc0, xd0, xa1, xb1, xc1, xd1, xa2, xb2, xc2, xd2, xa3, xb3, xc3, xd3
    vadd.vv \xa0, \xa0, \xb0
    vadd.vv \xa1, \xa1, \xb1
    vadd.vv \xa2, \xa2, \xb2
    vadd.vv \xa3, \xa3, \xb3
    vxor.vv \xd0, \xd0, \xa0
    vxor.vv \xd1, \xd1, \xa1
    vxor.vv \xd2, \xd2, \xa2
    vxor.vv \xd3, \xd3, \xa3
    vsll.vi v16, \xd0, 16
    vsll.vi v17, \xd1, 16
    vsll.vi v18, \xd2, 16
    vsll.vi v19, \xd3, 16
    vsrl.vi \xd0, \xd0, 16
    vsrl.vi \xd1, \xd1, 16
    vsrl.vi \xd2, \xd2, 16
    vsrl.vi \xd3, \xd3, 16
    vor.vv  \xd0, \xd0, v16
    vor.vv  \xd1, \xd1, v17
    vor.vv  \xd2, \xd2, v18
    vor.vv  \xd3, \xd3, v19

    vadd.vv \xc0, \xc0, \xd0
    vadd.vv \xc1, \xc1, \xd1
    vadd.vv \xc2, \xc2, \xd2
    vadd.vv \xc3, \xc3, \xd3
    vxor.vv \xb0, \xb0, \xc0
    vxor.vv \xb1, \xb1, \xc1
    vxor.vv \xb2, \xb2, \xc2
    vxor.vv \xb3, \xb3, \xc3
    vsll.vi v16, \xb0, 12
    vsll.vi v17, \xb1, 12
    vsll.vi v18, \xb2, 12
    vsll.vi v19, \xb3, 12
    vsrl.vi \xb0, \xb0, 20
    vsrl.vi \xb1, \xb1, 20
    vsrl.vi \xb2, \xb2, 20
    vsrl.vi \xb3, \xb3, 20
    vor.vv  \xb0, \xb0, v16
    vor.vv  \xb1, \xb1, v17
    vor.vv  \xb2, \xb2, v18
    vor.vv  \xb3, \xb3, v19

    vadd.vv \xa0, \xa0, \xb0
    vadd.vv \xa1, \xa1, \xb1
    vadd.vv \xa2, \xa2, \xb2
    vadd.vv \xa3, \xa3, \xb3
    vxor.vv \xd0, \xd0, \xa0
    vxor.vv \xd1, \xd1, \xa1
    vxor.vv \xd2, \xd2, \xa2
    vxor.vv \xd3, \xd3, \xa3
    vsll.vi v16, \xd0, 8
    vsll.vi v17, \xd1, 8
    vsll.vi v18, \xd2, 8
    vsll.vi v19, \xd3, 8
    vsrl.vi \xd0, \xd0, 24
    vsrl.vi \xd1, \xd1, 24
    vsrl.vi \xd2, \xd2, 24
    vsrl.vi \xd3, \xd3, 24
    vor.vv  \xd0, \xd0, v16
    vor.vv  \xd1, \xd1, v17
    vor.vv  \xd2, \xd2, v18
    vor.vv  \xd3, \xd3, v19

    vadd.vv \xc0, \xc0, \xd0
    vadd.vv \xc1, \xc1, \xd1
    vadd.vv \xc2, \xc2, \xd2
    vadd.vv \xc3, \xc3, \xd3
    vxor.vv \xb0, \xb0, \xc0
    vxor.vv \xb1, \xb1, \xc1
    vxor.vv \xb2, \xb2, \xc2
    vxor.vv \xb3, \xb3, \xc3
    vsll.vi v16, \xb0, 7
    vsll.vi v17, \xb1, 7
    vsll.vi v18, \xb2, 7
    vsll.vi v19, \xb3, 7
    vsrl.vi \xb0, \xb0, 25
    vsrl.vi \xb1, \xb1, 25
    vsrl.vi \xb2, \xb2, 25
    vsrl.vi \xb3, \xb3, 25
    vor.vv  \xb0, \xb0, v16
    vor.vv  \xb1, \xb1, v17
    vor.vv  \xb2, \xb2, v18
    vor.vv  \xb3, \xb3, v19
.endm

chacha20_rvv_blocks:
    beqz a2, .Lblocks_done

    addi sp, sp, -576
    sd   ra, 568(sp)
    sd   s0, 512(sp)
    sd   s1, 520(sp)
    sd   s2, 528(sp)
    sd   s3, 536(sp)
    sd   s4, 544(sp)
    sd   s5, 552(sp)
    sd   s6, 560(sp)
    sd   s7, 0(sp)
    sd   s8, 8(sp)

    mv   s0, a2                 // group count
    mv   s2, a4                 // counter
    mv   s4, a0                 // out
    mv   s5, a1                 // in

    vsetivli t0, 8, e32, m1, ta, ma

    // Keep the fifteen broadcast state words in integer registers; the
    // feed-forward then uses vadd.vx and needs no state memory traffic.
    lw   s3, 0(a5)
    lw   s7, 4(a5)
    lw   s8, 8(a5)

    lw   a4, 0(a3)
    lw   a5, 4(a3)
    lw   a6, 8(a3)
    lw   a7, 12(a3)
    lw   t4, 16(a3)
    lw   t5, 20(a3)
    lw   t6, 24(a3)
    lw   s1, 28(a3)

    li   a0, 0x61707865
    li   a1, 0x3320646e
    li   a2, 0x79622d32
    li   a3, 0x6b206574

    vmv.v.x v0, a0
    vmv.v.x v1, a1
    vmv.v.x v2, a2
    vmv.v.x v3, a3
    vmv.v.x v4, a4
    vmv.v.x v5, a5
    vmv.v.x v6, a6
    vmv.v.x v7, a7
    vmv.v.x v8, t4
    vmv.v.x v9, t5
    vmv.v.x v10, t6
    vmv.v.x v11, s1
    vmv.v.x v13, s3
    vmv.v.x v14, s7
    vmv.v.x v15, s8

    // First group: the initial state is already live in v0..v15.
    vid.v v12
    vadd.vx v12, v12, s2

.Lround_body:
    li   s6, 10
.Lround:
    QR4 v0, v4, v8, v12, v1, v5, v9, v13, v2, v6, v10, v14, v3, v7, v11, v15
    QR4 v0, v5, v10, v15, v1, v6, v11, v12, v2, v7, v8, v13, v3, v4, v9, v14
    addi s6, s6, -1
    bnez s6, .Lround

    // Feed-forward with the scalar register copies: vadd.vx adds the
    // broadcast input word to every lane without touching memory.
    vadd.vx v0, v0, a0
    vadd.vx v1, v1, a1
    vadd.vx v2, v2, a2
    vadd.vx v3, v3, a3
    vadd.vx v4, v4, a4
    vadd.vx v5, v5, a5
    vadd.vx v6, v6, a6
    vadd.vx v7, v7, a7
    vadd.vx v8, v8, t4
    vadd.vx v9, v9, t5
    vadd.vx v10, v10, t6
    vadd.vx v11, v11, s1
    vid.v v16
    vadd.vx v16, v16, s2
    vadd.vv v12, v12, v16
    vadd.vx v13, v13, s3
    vadd.vx v14, v14, s7
    vadd.vx v15, v15, s8

    li   t0, 64

    // Strided segment loads: eight 16-word blocks, stride 64 bytes.
    // The first segment covers words 0..7 of every block, the second
    // covers words 8..15.  v16..v23 and v24..v31 are input/XOR temps.
    vlsseg8e32.v v16, (s5), t0
    addi t1, s5, 32
    vlsseg8e32.v v24, (t1), t0

    vxor.vv v16, v16, v0
    vxor.vv v17, v17, v1
    vxor.vv v18, v18, v2
    vxor.vv v19, v19, v3
    vxor.vv v20, v20, v4
    vxor.vv v21, v21, v5
    vxor.vv v22, v22, v6
    vxor.vv v23, v23, v7

    vxor.vv v24, v24, v8
    vxor.vv v25, v25, v9
    vxor.vv v26, v26, v10
    vxor.vv v27, v27, v11
    vxor.vv v28, v28, v12
    vxor.vv v29, v29, v13
    vxor.vv v30, v30, v14
    vxor.vv v31, v31, v15

    vssseg8e32.v v16, (s4), t0
    addi t2, s4, 32
    vssseg8e32.v v24, (t2), t0

    addi s4, s4, 512
    addi s5, s5, 512
    addi s2, s2, 8
    addi s0, s0, -1
    beqz s0, .Ldone_loop

    // Prepare the next group from the scalar register copies.
    vmv.v.x v0, a0
    vmv.v.x v1, a1
    vmv.v.x v2, a2
    vmv.v.x v3, a3
    vmv.v.x v4, a4
    vmv.v.x v5, a5
    vmv.v.x v6, a6
    vmv.v.x v7, a7
    vmv.v.x v8, t4
    vmv.v.x v9, t5
    vmv.v.x v10, t6
    vmv.v.x v11, s1
    vid.v v12
    vadd.vx v12, v12, s2
    vmv.v.x v13, s3
    vmv.v.x v14, s7
    vmv.v.x v15, s8
    j .Lround_body

.Ldone_loop:
    ld   ra, 568(sp)
    ld   s0, 512(sp)
    ld   s1, 520(sp)
    ld   s2, 528(sp)
    ld   s3, 536(sp)
    ld   s4, 544(sp)
    ld   s5, 552(sp)
    ld   s6, 560(sp)
    ld   s7, 0(sp)
    ld   s8, 8(sp)
    addi sp, sp, 576

.Lblocks_done:
    ret
