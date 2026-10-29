/*
 * Copyright (c) 2026 harmor.
 *
 * Fixed 3-limb Montgomery multiply for P-192.  Generated from the CIOS C
 * reference with the project RV64 flags; the loop bounds are compile-time
 * constants so the body is fully unrolled.
 */
	.file	"gen_mont_mul_fixed.c"
	.option nopic
	.attribute arch, "rv64i2p1_m2p0_a2p1_f2p2_d2p2_c2p0_b1p0_zicsr2p0_zifencei2p0_zmmul1p0_zaamo1p0_zalrsc1p0_zca1p0_zcd1p0_zba1p0_zbb1p0_zbc1p0_zbs1p0"
	.attribute unaligned_access, 0
	.attribute stack_align, 16
	.text
	.align	1
	.globl	hitls_mont_mul_3_rv64
	.type	hitls_mont_mul_3_rv64, @function
hitls_mont_mul_3_rv64:
.LFB0:
	.cfi_startproc
	ld	t2,0(a1)
	ld	t3,0(a2)
	addi	sp,sp,-32
	.cfi_def_cfa_offset 32
	sd	s1,16(sp)
	mul	a5,t2,t3
	ld	t6,0(a3)
	ld	t1,8(a3)
	ld	a6,16(a3)
	ld	t0,8(a1)
	ld	t4,16(a1)
	sd	s0,24(sp)
	sd	s2,8(sp)
	sd	s3,0(sp)
	.cfi_offset 9, -16
	.cfi_offset 8, -8
	.cfi_offset 18, -24
	.cfi_offset 19, -32
	ld	a7,8(a2)
	mul	s1,a4,a5
	ld	t5,16(a2)
	mul	a3,s1,t6
	mulhu	s0,t2,t3
	add	a5,a3,a5
	sltu	a5,a5,a3
	mul	a1,t3,t0
	mul	s2,s1,t1
	add	a1,s0,a1
	sltu	s0,a1,s0
	mulhu	a3,s1,t6
	add	a1,s2,a1
	add	a5,a5,a3
	add	a3,a5,a1
	sltu	a1,a1,s2
	mulhu	s2,s1,t1
	sltu	a5,a3,a5
	add	a1,a1,s2
	mulhu	s2,t3,t0
	add	a5,a5,a1
	add	s0,s0,s2
	mul	s2,t3,t4
	mulhu	t3,t3,t4
	add	s0,s2,s0
	sltu	s2,s0,s2
	mul	s3,s1,a6
	add	s2,s2,t3
	mulhu	t3,s1,a6
	add	s0,s3,s0
	sltu	s3,s0,s3
	add	s0,a5,s0
	sltu	a5,s0,a5
	add	s3,s3,t3
	add	a5,a5,s3
	mul	s3,t4,a7
	add	a5,s2,a5
	sltu	s1,a5,s2
	mulhu	t3,t4,a7
	add	s3,a5,s3
	sltu	a5,s3,a5
	mul	a2,t2,a7
	add	t3,a5,t3
	mulhu	a5,t2,a7
	add	a2,a3,a2
	sltu	a3,a2,a3
	add	a5,a3,a5
	mul	a3,t0,a7
	add	a5,s0,a5
	sltu	s0,a5,s0
	mulhu	s2,t0,a7
	add	a5,a3,a5
	sltu	a3,a5,a3
	mul	a1,a4,a2
	add	a7,s2,s0
	add	a3,a3,a7
	add	s3,a3,s3
	sltu	a3,s3,a3
	add	a3,a3,t3
	add	t3,s1,a3
	sltu	a7,t3,s1
	mul	s0,a6,a1
	mulhu	a3,a6,a1
	add	s1,s0,s3
	sltu	s0,s1,s0
	ld	s3,0(sp)
	.cfi_restore 19
	mul	s2,t1,a1
	add	s0,s0,a3
	mulhu	a3,t1,a1
	add	a5,s2,a5
	sltu	s2,a5,s2
	add	s2,s2,a3
	mul	a3,t6,a1
	mulhu	a1,t6,a1
	add	a2,a3,a2
	sltu	a2,a2,a3
	mul	a3,t2,t5
	add	a2,a2,a1
	add	a5,a2,a5
	sltu	a2,a5,a2
	add	a2,a2,s2
	add	s1,a2,s1
	sltu	a2,s1,a2
	add	a2,a2,s0
	add	s0,t3,a2
	sltu	a2,s0,t3
	add	a3,a5,a3
	mul	t3,a4,a3
	sltu	a5,a3,a5
	add	a7,a2,a7
	ld	s2,8(sp)
	.cfi_restore 18
	mulhu	t2,t2,t5
	mul	a4,t6,t3
	add	a5,a5,t2
	add	a5,s1,a5
	sltu	s1,a5,s1
	mul	a1,t0,t5
	add	a2,a4,a3
	sltu	a2,a2,a4
	mulhu	t0,t0,t5
	add	a5,a1,a5
	sltu	a1,a5,a1
	mulhu	t2,t6,t3
	add	t0,t0,s1
	add	a1,a1,t0
	add	a1,s0,a1
	sltu	t0,a1,s0
	mul	a4,t1,t3
	add	a2,a2,t2
	mul	a3,t4,t5
	add	a5,a4,a5
	sltu	a4,a5,a4
	add	a5,a2,a5
	sltu	a2,a5,a2
	sub	t6,a5,t6
	sgtu	s0,t6,a5
	mulhu	s1,t1,t3
	add	a1,a3,a1
	sltu	a3,a1,a3
	mul	t2,a6,t3
	add	a4,a4,s1
	add	a2,a2,a4
	ld	s1,16(sp)
	.cfi_restore 9
	mulhu	t4,t4,t5
	add	a1,t2,a1
	add	t5,a2,a1
	sub	t1,t5,t1
	sltu	a1,a1,t2
	sltu	a4,t5,a2
	sub	a2,t1,s0
	sgtu	t2,t1,t5
	sgtu	t1,a2,t1
	sub	t1,t2,t1
	mulhu	t3,a6,t3
	add	t4,t4,t0
	add	a3,a3,t4
	add	a3,a7,a3
	andi	t1,t1,1
	sltu	a7,a3,a7
	ld	s0,24(sp)
	.cfi_restore 8
	add	a1,a1,t3
	add	a4,a4,a1
	add	a4,a3,a4
	sub	a6,a4,a6
	sub	a1,a6,t1
	sltu	a3,a4,a3
	sgtu	t1,a6,a4
	sgtu	a6,a1,a6
	add	a3,a3,a7
	sub	a6,t1,a6
	li	a7,1
	andn	a6,a7,a6
	snez	a3,a3
	or	a3,a3,a6
	addi	a6,a3,-1
	neg	a3,a3
	and	a5,a6,a5
	and	t6,t6,a3
	and	a2,a3,a2
	and	t5,a6,t5
	and	a3,a3,a1
	and	a6,a6,a4
	or	t6,t6,a5
	or	a2,a2,t5
	or	a3,a3,a6
	sd	t6,0(a0)
	sd	a2,8(a0)
	sd	a3,16(a0)
	addi	sp,sp,32
	.cfi_def_cfa_offset 0
	jr	ra
	.cfi_endproc
.LFE0:
	.size	hitls_mont_mul_3_rv64, .-hitls_mont_mul_3_rv64
	.ident	"GCC: (g5115c7e44) 15.2.0"
	.section	.note.GNU-stack,"",@progbits
