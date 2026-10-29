/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */

/*
 * P-256 (secp256r1) field, point and ECDSA/ECDH operations.
 *
 * Field and scalar arithmetic use four 64-bit limbs in Montgomery form.
 * Points use Jacobian coordinates; scalar multiplication uses a 4-bit fixed
 * window with a full-table constant-time selection.
 */

#include "p256_accel_api.h"

#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;

typedef struct {
    uint64_t v[4];
} u256;

typedef u256 fe;
typedef u256 sc;

typedef struct {
    fe x, y, z;
} jac_point;

#if defined(__riscv) && defined(P256_USE_ASM)
void p256_mont_mul_p_rv64(uint64_t r[4], const uint64_t a[4],
                          const uint64_t b[4]);
void p256_mont_mul_n_rv64(uint64_t r[4], const uint64_t a[4],
                          const uint64_t b[4]);
#endif

/* P-256 domain parameters, little-endian 64-bit limbs. */
static const uint64_t P[4] = {
    0xffffffffffffffffull, 0x00000000ffffffffull,
    0x0000000000000000ull, 0xffffffff00000001ull
};
static const uint64_t N[4] = {
    0xf3b9cac2fc632551ull, 0xbce6faada7179e84ull,
    0xffffffffffffffffull, 0xffffffff00000000ull
};
static const uint64_t B[4] = {
    0x3bce3c3e27d2604bull, 0x651d06b0cc53b0f6ull,
    0xb3ebbd55769886bcull, 0x5ac635d8aa3a93e7ull
};
static const uint64_t GX[4] = {
    0xf4a13945d898c296ull, 0x77037d812deb33a0ull,
    0xf8bce6e563a440f2ull, 0x6b17d1f2e12c4247ull
};
static const uint64_t GY[4] = {
    0xcbb6406837bf51f5ull, 0x2bce33576b315eceull,
    0x8ee7eb4a7c0f9e16ull, 0x4fe342e2fe1a7f9bull
};
static const uint64_t P_MINUS_2[4] = {
    0xfffffffffffffffdull, 0x00000000ffffffffull,
    0x0000000000000000ull, 0xffffffff00000001ull
};
static const uint64_t N_MINUS_2[4] = {
    0xf3b9cac2fc63254full, 0xbce6faada7179e84ull,
    0xffffffffffffffffull, 0xffffffff00000000ull
};
static const uint64_t P_R[4] = {
    0x0000000000000001ull, 0xffffffff00000000ull,
    0xffffffffffffffffull, 0x00000000fffffffeull
};
static const uint64_t P_R2[4] = {
    0x0000000000000003ull, 0xfffffffbffffffffull,
    0xfffffffffffffffeull, 0x00000004fffffffdull
};
static const uint64_t N_R[4] = {
    0x0c46353d039cdaafull, 0x4319055258e8617bull,
    0x0000000000000000ull, 0x00000000ffffffffull
};
static const uint64_t N_R2[4] = {
    0x83244c95be79eea2ull, 0x4699799c49bd6fa6ull,
    0x2845b2392b6bec59ull, 0x66e12d94f3d95620ull
};

#define P_N0 0x0000000000000001ull
#define N_N0 0xccd1c8aaee00bc4full

static void u256_from_bytes(u256 *r, const uint8_t b[32])
{
    int i;
    int j;

    for (i = 0; i < 4; i++) {
        uint64_t v = 0;
        const uint8_t *p = b + 8 * (3 - i);

        for (j = 0; j < 8; j++) {
            v = (v << 8) | p[j];
        }
        r->v[i] = v;
    }
}

static void u256_to_bytes(uint8_t b[32], const u256 *a)
{
    int i;
    int j;

    for (i = 0; i < 4; i++) {
        uint64_t v = a->v[3 - i];
        for (j = 7; j >= 0; j--) {
            b[i * 8 + j] = (uint8_t)(v & 0xff);
            v >>= 8;
        }
    }
}

static int u256_cmp_arr(const u256 *a, const uint64_t b[4])
{
    uint64_t borrow = 0;
    uint64_t diff = 0;
    uint64_t lt;
    uint64_t eq;
    int i;

    for (i = 0; i < 4; i++) {
        u128 uv = (u128)a->v[i] - b[i] - borrow;

        diff |= a->v[i] ^ b[i];
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    lt = borrow;
    eq = ((diff | (0 - diff)) >> 63) ^ 1;
    return (int)((uint64_t)(1 - lt - eq) - lt);
}

static int u256_is_zero(const u256 *a)
{
    return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0;
}

static int u256_equal(const u256 *a, const u256 *b)
{
    return ((a->v[0] ^ b->v[0]) | (a->v[1] ^ b->v[1]) |
            (a->v[2] ^ b->v[2]) | (a->v[3] ^ b->v[3])) == 0;
}

/*
 * Montgomery multiplication (CIOS style, 512-bit intermediate).  Inputs and
 * outputs are normalized to [0, mod).  N0 is -mod^-1 mod 2^64.
 */
#define DEFINE_MONT_MUL(fn, MOD, N0)                                         \
static void fn(u256 *r, const u256 *a, const u256 *b)                        \
{                                                                            \
    uint64_t t[12] = { 0 };                                                  \
    int i;                                                                   \
    int j;                                                                   \
    int k;                                                                   \
                                                                             \
    for (i = 0; i < 4; i++) {                                                \
        u128 carry = 0;                                                      \
        for (j = 0; j < 4; j++) {                                            \
            u128 uv = (u128)a->v[i] * b->v[j] + t[i + j] + carry;            \
            t[i + j] = (uint64_t)uv;                                         \
            carry = uv >> 64;                                                \
        }                                                                    \
        for (k = i + 4; k < 12; k++) {                                       \
            u128 uv = (u128)t[k] + (uint64_t)carry;                          \
            t[k] = (uint64_t)uv;                                             \
            carry = uv >> 64;                                                \
        }                                                                    \
    }                                                                        \
                                                                             \
    for (i = 0; i < 4; i++) {                                                \
        uint64_t m = t[i] * (uint64_t)(N0);                                  \
        u128 carry = 0;                                                      \
        for (j = 0; j < 4; j++) {                                            \
            u128 uv = (u128)m * (MOD)[j] + t[i + j] + carry;                 \
            t[i + j] = (uint64_t)uv;                                         \
            carry = uv >> 64;                                                \
        }                                                                    \
        for (k = i + 4; k < 12; k++) {                                       \
            u128 uv = (u128)t[k] + (uint64_t)carry;                          \
            t[k] = (uint64_t)uv;                                             \
            carry = uv >> 64;                                                \
        }                                                                    \
    }                                                                        \
                                                                             \
    {                                                                        \
        uint64_t d[4];                                                       \
        uint64_t borrow = 0;                                                 \
        uint64_t need;                                                       \
        uint64_t mask;                                                       \
                                                                             \
        for (i = 0; i < 4; i++) {                                            \
            u128 sub = (u128)t[4 + i] - (MOD)[i] - borrow;                   \
            d[i] = (uint64_t)sub;                                            \
            borrow = (uint64_t)(sub >> 64) & 1;                              \
        }                                                                    \
        need = (uint64_t)(t[8] != 0) | (borrow ^ 1);                         \
        mask = 0 - need;                                                     \
        for (i = 0; i < 4; i++) {                                            \
            r->v[i] = (d[i] & mask) | (t[4 + i] & ~mask);                    \
        }                                                                    \
    }                                                                        \
}

#if !(defined(__riscv) && defined(P256_USE_ASM))
DEFINE_MONT_MUL(mont_mul_p, P, P_N0)
#endif
#if defined(__riscv) && defined(P256_USE_ASM)
static void mont_mul_n(sc *r, const sc *a, const sc *b)
{
    p256_mont_mul_n_rv64(r->v, a->v, b->v);
}
#else
DEFINE_MONT_MUL(mont_mul_n, N, N_N0)
#endif

static void fe_add(fe *r, const fe *a, const fe *b)
{
    uint64_t s[4];
    uint64_t d[4];
    uint64_t carry;
    uint64_t borrow = 0;
    u128 uv;
    int i;
    uint64_t need;
    uint64_t mask;

    uv = (u128)a->v[0] + b->v[0];
    s[0] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);
    uv = (u128)a->v[1] + b->v[1] + carry;
    s[1] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);
    uv = (u128)a->v[2] + b->v[2] + carry;
    s[2] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);
    uv = (u128)a->v[3] + b->v[3] + carry;
    s[3] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);

    for (i = 0; i < 4; i++) {
        uv = (u128)s[i] - P[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    need = carry | (borrow ^ 1);
    mask = 0 - need;
    for (i = 0; i < 4; i++) {
        r->v[i] = (d[i] & mask) | (s[i] & ~mask);
    }
}

static void fe_sub(fe *r, const fe *a, const fe *b)
{
    uint64_t d[4];
    uint64_t borrow = 0;
    uint64_t carry;
    uint64_t mask;
    u128 uv;
    int i;

    for (i = 0; i < 4; i++) {
        uv = (u128)a->v[i] - b->v[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    mask = 0 - borrow;
    carry = 0;
    for (i = 0; i < 4; i++) {
        uv = (u128)d[i] + (P[i] & mask) + carry;
        r->v[i] = (uint64_t)uv;
        carry = (uint64_t)(uv >> 64);
    }
}

static void fe_double(fe *r, const fe *a)
{
    fe_add(r, a, a);
}

static void fe_triple(fe *r, const fe *a)
{
    fe t;

    fe_double(&t, a);
    fe_add(r, &t, a);
}

static void fe_mul4(fe *r, const fe *a)
{
    fe t;

    fe_double(&t, a);
    fe_double(r, &t);
}

static void fe_mul8(fe *r, const fe *a)
{
    fe t;

    fe_mul4(&t, a);
    fe_double(r, &t);
}

static void fe_mul(fe *r, const fe *a, const fe *b)
{
#if defined(__riscv) && defined(P256_USE_ASM)
    p256_mont_mul_p_rv64(r->v, a->v, b->v);
#else
    mont_mul_p(r, a, b);
#endif
}

static void fe_sqr(fe *r, const fe *a)
{
    fe_mul(r, a, a);
}

static void fe_to_mont(fe *r, const fe *a)
{
    fe r2;
    memcpy(r2.v, P_R2, sizeof(r2.v));
    fe_mul(r, a, &r2);
}

static void fe_from_mont(fe *r, const fe *a)
{
    fe one = { { 1, 0, 0, 0 } };
    fe_mul(r, a, &one);
}

static void fe_inv(fe *r, const fe *a)
{
    fe res;
    fe base = *a;
    int i;

    memcpy(res.v, P_R, sizeof(res.v));
    for (i = 255; i >= 0; i--) {
        fe_sqr(&res, &res);
        if ((P_MINUS_2[i >> 6] >> (i & 63)) & 1) {
            fe_mul(&res, &res, &base);
        }
    }
    *r = res;
}

static int fe_is_zero(const fe *a)
{
    return u256_is_zero(a);
}

static int fe_is_equal(const fe *a, const fe *b)
{
    return u256_equal(a, b);
}

/* All-ones mask if a == 0, zero otherwise. */
static uint64_t ct_zero_mask(const fe *a)
{
    uint64_t x = a->v[0] | a->v[1] | a->v[2] | a->v[3];

    x = (x | (0 - x)) >> 63;
    return x - 1;
}

/* r = mask ? b : a */
static void point_select_mask(jac_point *r, const jac_point *a,
                              const jac_point *b, uint64_t mask)
{
    int j;

    for (j = 0; j < 4; j++) {
        r->x.v[j] = (a->x.v[j] & ~mask) | (b->x.v[j] & mask);
        r->y.v[j] = (a->y.v[j] & ~mask) | (b->y.v[j] & mask);
        r->z.v[j] = (a->z.v[j] & ~mask) | (b->z.v[j] & mask);
    }
}

static int fe_parse(fe *r, const uint8_t b[32])
{
    u256_from_bytes(r, b);
    return u256_cmp_arr(r, P) < 0 ? 0 : -1;
}

static void sc_add(sc *r, const sc *a, const sc *b)
{
    uint64_t s[4];
    uint64_t d[4];
    uint64_t carry;
    uint64_t borrow = 0;
    u128 uv;
    int i;
    uint64_t need;
    uint64_t mask;

    uv = (u128)a->v[0] + b->v[0];
    s[0] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);
    uv = (u128)a->v[1] + b->v[1] + carry;
    s[1] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);
    uv = (u128)a->v[2] + b->v[2] + carry;
    s[2] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);
    uv = (u128)a->v[3] + b->v[3] + carry;
    s[3] = (uint64_t)uv;
    carry = (uint64_t)(uv >> 64);

    for (i = 0; i < 4; i++) {
        uv = (u128)s[i] - N[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    need = carry | (borrow ^ 1);
    mask = 0 - need;
    for (i = 0; i < 4; i++) {
        r->v[i] = (d[i] & mask) | (s[i] & ~mask);
    }
}

static void sc_mul(sc *r, const sc *a, const sc *b)
{
    mont_mul_n(r, a, b);
}

static void sc_sqr(sc *r, const sc *a)
{
    mont_mul_n(r, a, a);
}

static void sc_to_mont(sc *r, const sc *a)
{
    sc r2;
    memcpy(r2.v, N_R2, sizeof(r2.v));
    mont_mul_n(r, a, &r2);
}

static void sc_from_mont(sc *r, const sc *a)
{
    sc one = { { 1, 0, 0, 0 } };
    mont_mul_n(r, a, &one);
}

static void sc_inv(sc *r, const sc *a)
{
    sc res;
    sc base = *a;
    int i;

    memcpy(res.v, N_R, sizeof(res.v));
    for (i = 255; i >= 0; i--) {
        sc_sqr(&res, &res);
        if ((N_MINUS_2[i >> 6] >> (i & 63)) & 1) {
            sc_mul(&res, &res, &base);
        }
    }
    *r = res;
}

static void sc_reduce_once(sc *r)
{
    uint64_t d[4];
    uint64_t borrow = 0;
    uint64_t mask;
    int i;

    for (i = 0; i < 4; i++) {
        u128 uv = (u128)r->v[i] - N[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    mask = 0 - (borrow ^ 1);
    for (i = 0; i < 4; i++) {
        r->v[i] = (d[i] & mask) | (r->v[i] & ~mask);
    }
}

static int sc_parse(sc *r, const uint8_t b[32], int nonzero)
{
    u256_from_bytes(r, b);
    if (u256_cmp_arr(r, N) >= 0) {
        return -1;
    }
    if (nonzero && u256_is_zero(r)) {
        return -1;
    }
    return 0;
}

static void point_set_infinity(jac_point *p)
{
    memset(p, 0, sizeof(*p));
}

static int point_is_infinity(const jac_point *p)
{
    return fe_is_zero(&p->z);
}

static void point_double(jac_point *r, const jac_point *p)
{
    fe delta, gamma, beta, alpha, t1, t2;
    fe x3, y3, z3;

    fe_sqr(&delta, &p->z);
    fe_sqr(&gamma, &p->y);
    fe_mul(&beta, &p->x, &gamma);

    fe_sub(&t1, &p->x, &delta);
    fe_add(&t2, &p->x, &delta);
    fe_mul(&alpha, &t1, &t2);
    fe_triple(&alpha, &alpha);

    fe_sqr(&t1, &alpha);
    fe_mul8(&t2, &beta);
    fe_sub(&x3, &t1, &t2);

    fe_add(&t1, &p->y, &p->z);
    fe_sqr(&t1, &t1);
    fe_sub(&t1, &t1, &gamma);
    fe_sub(&z3, &t1, &delta);

    fe_mul4(&t1, &beta);
    fe_sub(&t1, &t1, &x3);
    fe_mul(&t1, &alpha, &t1);
    fe_sqr(&t2, &gamma);
    fe_mul8(&t2, &t2);
    fe_sub(&y3, &t1, &t2);

    r->x = x3;
    r->y = y3;
    r->z = z3;
}

static void point_add(jac_point *r, const jac_point *p, const jac_point *q)
{
    fe z1z1, z2z2, u1, u2, s1, s2;
    fe h, rr, i, j, v, t1, t2;
    fe x3, y3, z3;
    jac_point add_out, dbl_out, inf_out, out;
    uint64_t mp = ct_zero_mask(&p->z);
    uint64_t mq = ct_zero_mask(&q->z);
    uint64_t hz;
    uint64_t rz;
    uint64_t same;
    uint64_t opp;

    fe_sqr(&z1z1, &p->z);
    fe_sqr(&z2z2, &q->z);
    fe_mul(&u1, &p->x, &z2z2);
    fe_mul(&u2, &q->x, &z1z1);

    fe_mul(&t1, &q->z, &z2z2);
    fe_mul(&s1, &p->y, &t1);
    fe_mul(&t1, &p->z, &z1z1);
    fe_mul(&s2, &q->y, &t1);

    fe_sub(&h, &u2, &u1);
    fe_sub(&t1, &s2, &s1);
    hz = ct_zero_mask(&h);
    rz = ct_zero_mask(&t1);

    fe_double(&t2, &h);
    fe_sqr(&i, &t2);
    fe_mul(&j, &h, &i);
    fe_double(&rr, &t1);
    fe_mul(&v, &u1, &i);

    fe_sqr(&x3, &rr);
    fe_sub(&x3, &x3, &j);
    fe_double(&t1, &v);
    fe_sub(&x3, &x3, &t1);

    fe_sub(&t1, &v, &x3);
    fe_mul(&t1, &rr, &t1);
    fe_mul(&t2, &s1, &j);
    fe_double(&t2, &t2);
    fe_sub(&y3, &t1, &t2);

    fe_add(&t1, &p->z, &q->z);
    fe_sqr(&t1, &t1);
    fe_sub(&t1, &t1, &z1z1);
    fe_sub(&t1, &t1, &z2z2);
    fe_mul(&z3, &t1, &h);

    add_out.x = x3;
    add_out.y = y3;
    add_out.z = z3;
    point_double(&dbl_out, p);
    point_set_infinity(&inf_out);

    same = hz & rz & ~mp & ~mq;
    opp = hz & ~rz & ~mp & ~mq;
    out = add_out;
    point_select_mask(&out, &out, &dbl_out, same);
    point_select_mask(&out, &out, &inf_out, opp);
    point_select_mask(&out, &out, q, mp);
    point_select_mask(&out, &out, p, mq);
    *r = out;
}

static void point_add_fast(jac_point *r, const jac_point *p, const jac_point *q)
{
    fe z1z1, z2z2, u1, u2, s1, s2;
    fe h, rr, i, j, v, t1, t2;
    fe x3, y3, z3;

    if (point_is_infinity(p)) {
        *r = *q;
        return;
    }
    if (point_is_infinity(q)) {
        *r = *p;
        return;
    }
    fe_sqr(&z1z1, &p->z);
    fe_sqr(&z2z2, &q->z);
    fe_mul(&u1, &p->x, &z2z2);
    fe_mul(&u2, &q->x, &z1z1);
    fe_mul(&t1, &q->z, &z2z2);
    fe_mul(&s1, &p->y, &t1);
    fe_mul(&t1, &p->z, &z1z1);
    fe_mul(&s2, &q->y, &t1);
    fe_sub(&h, &u2, &u1);
    fe_sub(&rr, &s2, &s1);
    if (fe_is_zero(&h)) {
        if (fe_is_zero(&rr)) {
            point_double(r, p);
        } else {
            point_set_infinity(r);
        }
        return;
    }
    fe_double(&t2, &h);
    fe_sqr(&i, &t2);
    fe_mul(&j, &h, &i);
    fe_double(&rr, &rr);
    fe_mul(&v, &u1, &i);
    fe_sqr(&x3, &rr);
    fe_sub(&x3, &x3, &j);
    fe_double(&t1, &v);
    fe_sub(&x3, &x3, &t1);
    fe_sub(&t1, &v, &x3);
    fe_mul(&t1, &rr, &t1);
    fe_mul(&t2, &s1, &j);
    fe_double(&t2, &t2);
    fe_sub(&y3, &t1, &t2);
    fe_add(&t1, &p->z, &q->z);
    fe_sqr(&t1, &t1);
    fe_sub(&t1, &t1, &z1z1);
    fe_sub(&t1, &t1, &z2z2);
    fe_mul(&z3, &t1, &h);
    r->x = x3;
    r->y = y3;
    r->z = z3;
}

static uint64_t ct_eq_mask(uint64_t a, uint64_t b)
{
    uint64_t x = a ^ b;

    x = (x | (0 - x)) >> 63;
    return x - 1;
}

static void point_select(jac_point *r, const jac_point table[16], unsigned idx)
{
    unsigned i;
    int j;

    point_set_infinity(r);
    for (i = 0; i < 16; i++) {
        uint64_t mask = ct_eq_mask(i, idx);

        for (j = 0; j < 4; j++) {
            r->x.v[j] |= table[i].x.v[j] & mask;
            r->y.v[j] |= table[i].y.v[j] & mask;
            r->z.v[j] |= table[i].z.v[j] & mask;
        }
    }
}

static void build_table(jac_point table[16], const jac_point *p)
{
    int i;

    point_set_infinity(&table[0]);
    table[1] = *p;
    point_double(&table[2], p);
    for (i = 3; i < 16; i++) {
        point_add(&table[i], &table[i - 1], &table[1]);
    }
}

static unsigned scalar_nibble(const uint64_t v[4], int idx)
{
    return (unsigned)((v[idx >> 4] >> ((idx & 15) * 4)) & 0xf);
}

static void scalar_mul_table(jac_point *r, const jac_point table[16],
                             const uint64_t k[4])
{
    jac_point acc;
    jac_point t;
    int i;

    point_set_infinity(&acc);
    for (i = 63; i >= 0; i--) {
        point_double(&acc, &acc);
        point_double(&acc, &acc);
        point_double(&acc, &acc);
        point_double(&acc, &acc);
        point_select(&t, table, scalar_nibble(k, i));
        point_add(&acc, &acc, &t);
    }
    *r = acc;
}

static int sc_is_zero(const uint64_t k[5])
{
    return (k[0] | k[1] | k[2] | k[3] | k[4]) == 0;
}

static void sc_sub_small(uint64_t k[5], int d)
{
    uint64_t v;
    int i;

    if (d < 0) {
        v = (uint64_t)(-d);
        for (i = 0; i < 5; i++) {
            uint64_t old = k[i];
            k[i] = old + v;
            v = (k[i] < old);
        }
        return;
    }
    v = (uint64_t)d;
    for (i = 0; i < 5; i++) {
        uint64_t old = k[i];
        k[i] = old - v;
        v = (old < v);
    }
}

static void sc_halve(uint64_t k[5])
{
    k[0] = (k[0] >> 1) | (k[1] << 63);
    k[1] = (k[1] >> 1) | (k[2] << 63);
    k[2] = (k[2] >> 1) | (k[3] << 63);
    k[3] = (k[3] >> 1) | (k[4] << 63);
    k[4] >>= 1;
}

static int sc_wnaf5(int8_t out[258], const uint64_t in[4])
{
    uint64_t k[5] = { in[0], in[1], in[2], in[3], 0 };
    int n = 0;
    int i;

    while (!sc_is_zero(k)) {
        int d = 0;
        if (k[0] & 1u) {
            d = (int)(k[0] & 31u);
            if (d > 16) {
                d -= 32;
            }
            sc_sub_small(k, d);
        }
        out[n++] = (int8_t)d;
        sc_halve(k);
    }
    for (i = n; i < 258; i++) {
        out[i] = 0;
    }
    return n;
}

static void point_neg(jac_point *r, const jac_point *p)
{
    fe zero;

    memset(&zero, 0, sizeof(zero));
    r->x = p->x;
    r->z = p->z;
    fe_sub(&r->y, &zero, &p->y);
}

static void scalar_mul_joint(jac_point *r, const jac_point gt[16],
                             const uint64_t u1[4], const jac_point qt[16],
                             const uint64_t u2[4])
{
    jac_point acc;
    int8_t n1[258];
    int8_t n2[258];
    int len1 = sc_wnaf5(n1, u1);
    int len2 = sc_wnaf5(n2, u2);
    int len = len1 > len2 ? len1 : len2;
    int i;

    /* ECDSA verification inputs are public, so use direct table indexing,
       signed digits and variable-time addition here; signing remains on the
       constant-time scalar_mul_table path. */
    point_set_infinity(&acc);
    for (i = len - 1; i >= 0; i--) {
        int d1 = n1[i];
        int d2 = n2[i];

        point_double(&acc, &acc);
        if (d1 != 0) {
            jac_point t;
            int a = d1 < 0 ? -d1 : d1;
            if (d1 < 0) {
                point_neg(&t, &gt[a]);
            } else {
                t = gt[a];
            }
            point_add_fast(&acc, &acc, &t);
        }
        if (d2 != 0) {
            jac_point t;
            int a = d2 < 0 ? -d2 : d2;
            if (d2 < 0) {
                point_neg(&t, &qt[a]);
            } else {
                t = qt[a];
            }
            point_add_fast(&acc, &acc, &t);
        }
    }
    *r = acc;
}

static void point_to_affine(fe *x, fe *y, const jac_point *p)
{
    fe zi, zi2, zi3;

    fe_inv(&zi, &p->z);
    fe_sqr(&zi2, &zi);
    fe_mul(&zi3, &zi2, &zi);
    fe_mul(x, &p->x, &zi2);
    fe_mul(y, &p->y, &zi3);
}

static u256 g_one_mont;
static u256 g_b_mont;
static jac_point g_table[16];
static int g_ready;

static void init_constants(void)
{
    u256 b;

    if (g_ready) {
        return;
    }
    memcpy(g_one_mont.v, P_R, sizeof(g_one_mont.v));
    memcpy(b.v, B, sizeof(b.v));
    fe_to_mont(&g_b_mont, &b);

    {
        u256 gx;
        u256 gy;
        jac_point g;

        memcpy(gx.v, GX, sizeof(gx.v));
        memcpy(gy.v, GY, sizeof(gy.v));
        fe_to_mont(&g.x, &gx);
        fe_to_mont(&g.y, &gy);
        g.z = g_one_mont;
        build_table(g_table, &g);
    }
    g_ready = 1;
}

static int point_is_on_curve(const fe *x, const fe *y)
{
    fe lhs, rhs, t;

    fe_sqr(&lhs, y);
    fe_sqr(&t, x);
    fe_mul(&rhs, &t, x);
    fe_triple(&t, x);
    fe_sub(&rhs, &rhs, &t);
    fe_add(&rhs, &rhs, &g_b_mont);
    return fe_is_equal(&lhs, &rhs);
}

static int parse_point(jac_point *p, const uint8_t qx[32], const uint8_t qy[32])
{
    fe x;
    fe y;

    if (fe_parse(&x, qx) != 0 || fe_parse(&y, qy) != 0) {
        return -1;
    }
    fe_to_mont(&p->x, &x);
    fe_to_mont(&p->y, &y);
    if (!point_is_on_curve(&p->x, &p->y)) {
        return -1;
    }
    p->z = g_one_mont;
    return 0;
}

int p256_public_key(uint8_t qx[32], uint8_t qy[32], const uint8_t d[32])
{
    sc dv;
    jac_point r;
    fe x;
    fe y;

    if (sc_parse(&dv, d, 1) != 0) {
        return -1;
    }
    init_constants();
    scalar_mul_table(&r, g_table, dv.v);
    if (point_is_infinity(&r)) {
        return -1;
    }
    point_to_affine(&x, &y, &r);
    fe_from_mont(&x, &x);
    fe_from_mont(&y, &y);
    u256_to_bytes(qx, &x);
    u256_to_bytes(qy, &y);
    return 0;
}

int p256_validate_public_key(const uint8_t qx[32], const uint8_t qy[32])
{
    jac_point q;
    jac_point table[16];
    jac_point r;

    init_constants();
    if (parse_point(&q, qx, qy) != 0) {
        return 0;
    }
    build_table(table, &q);
    scalar_mul_table(&r, table, N);
    return point_is_infinity(&r) ? 1 : 0;
}

int p256_ecdsa_sign(uint8_t r_out[32], uint8_t s_out[32], const uint8_t d[32],
                    const uint8_t z[32], const uint8_t k[32])
{
    sc dv, kv, zv;
    sc km, dm, zm, rm;
    sc rd, t, kinv, s, r;
    jac_point kg;
    fe x, y;
    fe xc;

    if (sc_parse(&dv, d, 1) != 0 || sc_parse(&kv, k, 1) != 0) {
        return -1;
    }
    u256_from_bytes(&zv, z);
    sc_reduce_once(&zv);
    init_constants();

    scalar_mul_table(&kg, g_table, kv.v);
    if (point_is_infinity(&kg)) {
        return -1;
    }
    point_to_affine(&x, &y, &kg);
    fe_from_mont(&xc, &x);
    r = xc;
    sc_reduce_once(&r);
    if (u256_is_zero(&r)) {
        return -1;
    }

    sc_to_mont(&km, &kv);
    sc_to_mont(&dm, &dv);
    sc_to_mont(&zm, &zv);
    sc_to_mont(&rm, &r);
    sc_mul(&rd, &rm, &dm);
    sc_add(&t, &zm, &rd);
    sc_inv(&kinv, &km);
    sc_mul(&s, &kinv, &t);
    sc_from_mont(&s, &s);
    if (u256_is_zero(&s)) {
        return -1;
    }

    u256_to_bytes(r_out, &r);
    u256_to_bytes(s_out, &s);
    return 0;
}

int p256_ecdsa_verify(const uint8_t qx[32], const uint8_t qy[32],
                      const uint8_t z[32], const uint8_t r_in[32],
                      const uint8_t s_in[32])
{
    sc r, s, zv, zm, sm, rm, w, u1, u2;
    jac_point q;
    jac_point qt[16];
    jac_point xp;
    fe x, y;
    sc rx;

    if (sc_parse(&r, r_in, 1) != 0 || sc_parse(&s, s_in, 1) != 0) {
        return 0;
    }
    u256_from_bytes(&zv, z);
    sc_reduce_once(&zv);
    init_constants();
    if (parse_point(&q, qx, qy) != 0) {
        return 0;
    }

    sc_to_mont(&sm, &s);
    sc_inv(&w, &sm);
    sc_to_mont(&zm, &zv);
    sc_to_mont(&rm, &r);
    sc_mul(&u1, &zm, &w);
    sc_mul(&u2, &rm, &w);
    sc_from_mont(&u1, &u1);
    sc_from_mont(&u2, &u2);

    build_table(qt, &q);
    scalar_mul_joint(&xp, g_table, u1.v, qt, u2.v);
    if (point_is_infinity(&xp)) {
        return 0;
    }

    point_to_affine(&x, &y, &xp);
    fe_from_mont(&x, &x);
    rx = x;
    sc_reduce_once(&rx);
    return u256_equal(&rx, &r) ? 1 : 0;
}

int p256_ecdh(uint8_t out[32], const uint8_t d[32], const uint8_t qx[32],
              const uint8_t qy[32])
{
    sc dv;
    jac_point q;
    jac_point table[16];
    jac_point r;
    fe x, y;

    if (sc_parse(&dv, d, 1) != 0) {
        return -1;
    }
    init_constants();
    if (parse_point(&q, qx, qy) != 0) {
        return -1;
    }
    build_table(table, &q);
    scalar_mul_table(&r, table, dv.v);
    if (point_is_infinity(&r)) {
        return -1;
    }
    point_to_affine(&x, &y, &r);
    fe_from_mont(&x, &x);
    u256_to_bytes(out, &x);
    return 0;
}
