/*
 * Copyright (c) 2026 harmor.
 *
 * RISC-V ECC acceleration implementation for openHiTLS.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 */

/*
 * Generic prime-curve implementation template for P-384/P-521.
 *
 * The including translation unit defines:
 *   CURVE_PREFIX, CURVE_LIMBS, CURVE_BYTES, CURVE_BITS,
 *   CURVE_WINDOW_BITS, CURVE_WINDOWS
 *   P, N, B, GX, GY, P_MINUS_2, N_MINUS_2, P_R, P_R2, N_R, N_R2,
 *   P_N0, N_N0
 * and includes its public header first.
 */

#include <stdint.h>
#include <string.h>

typedef unsigned __int128 u128;
typedef struct { uint64_t v[CURVE_LIMBS]; } fe;
typedef fe sc;
typedef struct { fe x, y, z; } jac_point;
typedef struct { fe x, y; } affine_point;

#define CURVE_TABLE_SIZE (1 << CURVE_WINDOW_BITS)

#define CURVE_CAT2(a, b) a##b
#define CURVE_CAT(a, b) CURVE_CAT2(a, b)
#define FN(name) CURVE_CAT(CURVE_PREFIX, CURVE_CAT(_, name))

static void u_from_bytes(fe *r, const uint8_t b[CURVE_BYTES])
{
    int i;

    memset(r, 0, sizeof(*r));
    for (i = 0; i < CURVE_BYTES; i++) {
        int bit = (CURVE_BYTES - 1 - i) * 8;
        r->v[bit >> 6] |= (uint64_t)b[i] << (bit & 63);
    }
}

static void u_to_bytes(uint8_t b[CURVE_BYTES], const fe *a)
{
    int i;

    for (i = 0; i < CURVE_BYTES; i++) {
        int bit = (CURVE_BYTES - 1 - i) * 8;
        b[i] = (uint8_t)(a->v[bit >> 6] >> (bit & 63));
    }
}

static int u_cmp_arr(const fe *a, const uint64_t b[CURVE_LIMBS])
{
    uint64_t borrow = 0;
    uint64_t diff = 0;
    uint64_t lt;
    uint64_t eq;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)a->v[i] - b[i] - borrow;

        diff |= a->v[i] ^ b[i];
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    lt = borrow;
    eq = ((diff | (0 - diff)) >> 63) ^ 1;
    return (int)((uint64_t)(1 - lt - eq) - lt);
}

static int u_is_zero(const fe *a)
{
    uint64_t x = 0;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        x |= a->v[i];
    }
    return x == 0;
}

static int u_equal(const fe *a, const fe *b)
{
    uint64_t x = 0;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        x |= a->v[i] ^ b->v[i];
    }
    return x == 0;
}

#if defined(CURVE_USE_RISCV_ASM)
extern void hitls_mont_mul_6_rv64(uint64_t *r, const uint64_t *a,
    const uint64_t *b, const uint64_t *mod, uint64_t n0);
#define DEFINE_MONT_MUL(fn, MOD, N0)                                         \
static void fn(fe *r, const fe *a, const fe *b)                              \
{                                                                            \
    hitls_mont_mul_6_rv64(r->v, a->v, b->v, (const uint64_t *)(MOD),         \
                          (uint64_t)(N0));                                   \
}
#else
#define DEFINE_MONT_MUL(fn, MOD, N0)                                         \
static void fn(fe *r, const fe *a, const fe *b)                              \
{                                                                            \
    uint64_t t[CURVE_LIMBS + 2] = { 0 };                                     \
    int i;                                                                   \
    int j;                                                                   \
                                                                             \
    for (i = 0; i < CURVE_LIMBS; i++) {                                      \
        uint64_t carry = 0;                                                  \
        u128 uv;                                                             \
                                                                             \
        for (j = 0; j < CURVE_LIMBS; j++) {                                  \
            uv = (u128)a->v[j] * b->v[i] + t[j] + carry;                     \
            t[j] = (uint64_t)uv;                                             \
            carry = (uint64_t)(uv >> 64);                                    \
        }                                                                    \
        uv = (u128)t[CURVE_LIMBS] + carry;                                   \
        t[CURVE_LIMBS] = (uint64_t)uv;                                       \
        t[CURVE_LIMBS + 1] = (uint64_t)(uv >> 64);                           \
                                                                             \
        {                                                                    \
            uint64_t m = t[0] * (uint64_t)(N0);                              \
                                                                             \
            uv = (u128)m * (MOD)[0] + t[0];                                  \
            carry = (uint64_t)(uv >> 64);                                    \
            for (j = 1; j < CURVE_LIMBS; j++) {                              \
                uv = (u128)m * (MOD)[j] + t[j] + carry;                      \
                t[j - 1] = (uint64_t)uv;                                     \
                carry = (uint64_t)(uv >> 64);                                \
            }                                                                \
            uv = (u128)t[CURVE_LIMBS] + carry;                               \
            t[CURVE_LIMBS - 1] = (uint64_t)uv;                               \
            t[CURVE_LIMBS] = t[CURVE_LIMBS + 1] +                          \
                             (uint64_t)(uv >> 64);                           \
            t[CURVE_LIMBS + 1] = 0;                                          \
        }                                                                    \
    }                                                                        \
                                                                             \
    {                                                                        \
        uint64_t d[CURVE_LIMBS];                                             \
        uint64_t borrow = 0;                                                 \
        uint64_t need;                                                       \
        uint64_t mask;                                                       \
                                                                             \
        for (i = 0; i < CURVE_LIMBS; i++) {                                  \
            u128 sub = (u128)t[i] - (MOD)[i] - borrow;                       \
            d[i] = (uint64_t)sub;                                            \
            borrow = (uint64_t)(sub >> 64) & 1;                              \
        }                                                                    \
        need = (uint64_t)(t[CURVE_LIMBS] != 0) | (borrow ^ 1);               \
        mask = 0 - need;                                                     \
        for (i = 0; i < CURVE_LIMBS; i++) {                                  \
            r->v[i] = (d[i] & mask) | (t[i] & ~mask);                        \
        }                                                                    \
    }                                                                        \
}

#endif

DEFINE_MONT_MUL(mont_mul_p, P, P_N0)
DEFINE_MONT_MUL(mont_mul_n, N, N_N0)

static void fe_add(fe *r, const fe *a, const fe *b)
{
    uint64_t s[CURVE_LIMBS];
    uint64_t d[CURVE_LIMBS];
    uint64_t carry = 0;
    uint64_t borrow = 0;
    uint64_t need;
    uint64_t mask;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)a->v[i] + b->v[i] + carry;
        s[i] = (uint64_t)uv;
        carry = (uint64_t)(uv >> 64);
    }
    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)s[i] - P[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    need = carry | (borrow ^ 1);
    mask = 0 - need;
    for (i = 0; i < CURVE_LIMBS; i++) {
        r->v[i] = (d[i] & mask) | (s[i] & ~mask);
    }
}

static void fe_sub(fe *r, const fe *a, const fe *b)
{
    uint64_t d[CURVE_LIMBS];
    uint64_t borrow = 0;
    uint64_t carry = 0;
    uint64_t mask;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)a->v[i] - b->v[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    mask = 0 - borrow;
    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)d[i] + (P[i] & mask) + carry;
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
    mont_mul_p(r, a, b);
}

static void fe_sqr(fe *r, const fe *a)
{
    mont_mul_p(r, a, a);
}

static void fe_to_mont(fe *r, const fe *a)
{
    fe r2;
    memcpy(r2.v, P_R2, sizeof(r2.v));
    mont_mul_p(r, a, &r2);
}

static void fe_from_mont(fe *r, const fe *a)
{
    fe one = { { 1 } };
    mont_mul_p(r, a, &one);
}

static void fe_inv(fe *r, const fe *a)
{
    fe res;
    fe base = *a;
    int i;

    memcpy(res.v, P_R, sizeof(res.v));
    for (i = CURVE_BITS - 1; i >= 0; i--) {
        fe_sqr(&res, &res);
        if ((P_MINUS_2[i >> 6] >> (i & 63)) & 1) {
            fe_mul(&res, &res, &base);
        }
    }
    *r = res;
}

static int fe_cmp_var(const fe *a, const fe *b)
{
    int i;

    for (i = CURVE_LIMBS - 1; i >= 0; i--) {
        if (a->v[i] != b->v[i]) {
            return a->v[i] > b->v[i] ? 1 : -1;
        }
    }
    return 0;
}

static int fe_is_one(const fe *a)
{
    int i;

    if (a->v[0] != 1) {
        return 0;
    }
    for (i = 1; i < CURVE_LIMBS; i++) {
        if (a->v[i] != 0) {
            return 0;
        }
    }
    return 1;
}

static void fe_shr1(fe *r, const fe *a)
{
    int i;

    for (i = 0; i < CURVE_LIMBS - 1; i++) {
        r->v[i] = (a->v[i] >> 1) | (a->v[i + 1] << 63);
    }
    r->v[CURVE_LIMBS - 1] = a->v[CURVE_LIMBS - 1] >> 1;
}

/* Variable-time binary extended GCD inverse.  Only used on public verify
 * inputs; secret-scalar operations keep the constant-time Fermat inverse. */
static void fe_inv_var(fe *r, const fe *a)
{
    fe u = *a;
    fe v;
    fe x1 = { { 1 } };
    fe x2 = { { 0 } };

    memcpy(v.v, P, sizeof(v.v));
    while (!fe_is_one(&u) && !fe_is_one(&v)) {
        while ((u.v[0] & 1u) == 0) {
            fe_shr1(&u, &u);
            if ((x1.v[0] & 1u) != 0) {
                uint64_t carry = 0;
                int i;

                for (i = 0; i < CURVE_LIMBS; i++) {
                    u128 uv = (u128)x1.v[i] + P[i] + carry;

                    x1.v[i] = (uint64_t)uv;
                    carry = (uint64_t)(uv >> 64);
                }
                for (i = 0; i < CURVE_LIMBS - 1; i++) {
                    x1.v[i] = (x1.v[i] >> 1) | (x1.v[i + 1] << 63);
                }
                x1.v[CURVE_LIMBS - 1] =
                    (x1.v[CURVE_LIMBS - 1] >> 1) | (carry << 63);
            } else {
                fe_shr1(&x1, &x1);
            }
        }
        while ((v.v[0] & 1u) == 0) {
            fe_shr1(&v, &v);
            if ((x2.v[0] & 1u) != 0) {
                uint64_t carry = 0;
                int i;

                for (i = 0; i < CURVE_LIMBS; i++) {
                    u128 uv = (u128)x2.v[i] + P[i] + carry;

                    x2.v[i] = (uint64_t)uv;
                    carry = (uint64_t)(uv >> 64);
                }
                for (i = 0; i < CURVE_LIMBS - 1; i++) {
                    x2.v[i] = (x2.v[i] >> 1) | (x2.v[i + 1] << 63);
                }
                x2.v[CURVE_LIMBS - 1] =
                    (x2.v[CURVE_LIMBS - 1] >> 1) | (carry << 63);
            } else {
                fe_shr1(&x2, &x2);
            }
        }
        if (fe_cmp_var(&u, &v) >= 0) {
            uint64_t borrow = 0;
            uint64_t carry;
            uint64_t mask;
            int i;

            for (i = 0; i < CURVE_LIMBS; i++) {
                u128 uv = (u128)u.v[i] - v.v[i] - borrow;

                u.v[i] = (uint64_t)uv;
                borrow = (uint64_t)(uv >> 64) & 1;
            }
            borrow = 0;
            for (i = 0; i < CURVE_LIMBS; i++) {
                u128 uv = (u128)x1.v[i] - x2.v[i] - borrow;

                x1.v[i] = (uint64_t)uv;
                borrow = (uint64_t)(uv >> 64) & 1;
            }
            mask = 0 - borrow;
            carry = 0;
            for (i = 0; i < CURVE_LIMBS; i++) {
                u128 uv = (u128)x1.v[i] + (P[i] & mask) + carry;

                x1.v[i] = (uint64_t)uv;
                carry = (uint64_t)(uv >> 64);
            }
        } else {
            uint64_t borrow = 0;
            uint64_t carry;
            uint64_t mask;
            int i;

            for (i = 0; i < CURVE_LIMBS; i++) {
                u128 uv = (u128)v.v[i] - u.v[i] - borrow;

                v.v[i] = (uint64_t)uv;
                borrow = (uint64_t)(uv >> 64) & 1;
            }
            borrow = 0;
            for (i = 0; i < CURVE_LIMBS; i++) {
                u128 uv = (u128)x2.v[i] - x1.v[i] - borrow;

                x2.v[i] = (uint64_t)uv;
                borrow = (uint64_t)(uv >> 64) & 1;
            }
            mask = 0 - borrow;
            carry = 0;
            for (i = 0; i < CURVE_LIMBS; i++) {
                u128 uv = (u128)x2.v[i] + (P[i] & mask) + carry;

                x2.v[i] = (uint64_t)uv;
                carry = (uint64_t)(uv >> 64);
            }
        }
    }

    {
        fe t;

        if (fe_is_one(&u)) {
            t = x1;
        } else {
            t = x2;
        }
        /* t = a^-1 * R^-1; convert back to Montgomery form a^-1 * R. */
        fe r2;
        memcpy(r2.v, P_R2, sizeof(r2.v));
        mont_mul_p(&t, &t, &r2);   /* a^-1 (plain) */
        fe_to_mont(r, &t);         /* a^-1 * R */
    }
}

static int fe_is_zero(const fe *a)
{
    return u_is_zero(a);
}

static int fe_is_equal(const fe *a, const fe *b)
{
    return u_equal(a, b);
}

static uint64_t ct_zero_mask(const fe *a)
{
    uint64_t x = 0;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        x |= a->v[i];
    }
    x = (x | (0 - x)) >> 63;
    return x - 1;
}

static void point_select_mask(jac_point *r, const jac_point *a,
                              const jac_point *b, uint64_t mask)
{
    int j;

    for (j = 0; j < CURVE_LIMBS; j++) {
        r->x.v[j] = (a->x.v[j] & ~mask) | (b->x.v[j] & mask);
        r->y.v[j] = (a->y.v[j] & ~mask) | (b->y.v[j] & mask);
        r->z.v[j] = (a->z.v[j] & ~mask) | (b->z.v[j] & mask);
    }
}

static int fe_parse(fe *r, const uint8_t b[CURVE_BYTES])
{
    u_from_bytes(r, b);
    return u_cmp_arr(r, P) < 0 ? 0 : -1;
}

static void sc_add(sc *r, const sc *a, const sc *b)
{
    uint64_t s[CURVE_LIMBS];
    uint64_t d[CURVE_LIMBS];
    uint64_t carry = 0;
    uint64_t borrow = 0;
    uint64_t need;
    uint64_t mask;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)a->v[i] + b->v[i] + carry;
        s[i] = (uint64_t)uv;
        carry = (uint64_t)(uv >> 64);
    }
    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)s[i] - N[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    need = carry | (borrow ^ 1);
    mask = 0 - need;
    for (i = 0; i < CURVE_LIMBS; i++) {
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
    sc one = { { 1 } };
    mont_mul_n(r, a, &one);
}

static void sc_inv(sc *r, const sc *a)
{
    sc res;
    sc base = *a;
    int i;

    memcpy(res.v, N_R, sizeof(res.v));
    for (i = CURVE_BITS - 1; i >= 0; i--) {
        sc_sqr(&res, &res);
        if ((N_MINUS_2[i >> 6] >> (i & 63)) & 1) {
            sc_mul(&res, &res, &base);
        }
    }
    *r = res;
}

static void sc_reduce_once(sc *r)
{
    uint64_t d[CURVE_LIMBS];
    uint64_t borrow = 0;
    uint64_t mask;
    int i;

    for (i = 0; i < CURVE_LIMBS; i++) {
        u128 uv = (u128)r->v[i] - N[i] - borrow;
        d[i] = (uint64_t)uv;
        borrow = (uint64_t)(uv >> 64) & 1;
    }
    mask = 0 - (borrow ^ 1);
    for (i = 0; i < CURVE_LIMBS; i++) {
        r->v[i] = (d[i] & mask) | (r->v[i] & ~mask);
    }
}

static int sc_parse(sc *r, const uint8_t b[CURVE_BYTES], int nonzero)
{
    u_from_bytes(r, b);
    if (u_cmp_arr(r, N) >= 0) {
        return -1;
    }
    if (nonzero && u_is_zero(r)) {
        return -1;
    }
    return 0;
}

static int sc_from_z(sc *r, const uint8_t z[CURVE_BYTES])
{
    u_from_bytes(r, z);
#if (CURVE_BITS % 64) != 0
    if ((r->v[CURVE_LIMBS - 1] >> (CURVE_BITS % 64)) != 0) {
        return -1;
    }
#endif
    sc_reduce_once(r);
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

static void point_set_affine(jac_point *r, const affine_point *p)
{
    fe one;

    memcpy(one.v, P_R, sizeof(one.v));
    r->x = p->x;
    r->y = p->y;
    r->z = one;
}

/* Variable-time Jacobian + affine addition, used by public-input verify. */
static void point_add_mixed(jac_point *r, const jac_point *p,
                            const affine_point *q)
{
    fe z1z1, u2, s2, h, rr, i, j, v, t, t2;
    fe x3, y3, z3;

    if (fe_is_zero(&p->z)) {
        point_set_affine(r, q);
        return;
    }
    fe_sqr(&z1z1, &p->z);
    fe_mul(&u2, &q->x, &z1z1);
    fe_mul(&s2, &q->y, &p->z);
    fe_mul(&s2, &s2, &z1z1);
    fe_sub(&h, &u2, &p->x);
    fe_sub(&rr, &s2, &p->y);
    if (fe_is_zero(&h)) {
        if (fe_is_zero(&rr)) {
            point_double(r, p);
        } else {
            point_set_infinity(r);
        }
        return;
    }
    fe_double(&i, &h);
    fe_sqr(&i, &i);
    fe_mul(&j, &h, &i);
    fe_double(&rr, &rr);
    fe_mul(&v, &p->x, &i);
    fe_sqr(&x3, &rr);
    fe_sub(&x3, &x3, &j);
    fe_double(&t, &v);
    fe_sub(&x3, &x3, &t);
    fe_sub(&t, &v, &x3);
    fe_mul(&t, &rr, &t);
    fe_mul(&t2, &p->y, &j);
    fe_double(&t2, &t2);
    fe_sub(&y3, &t, &t2);
    fe_add(&t, &p->z, &h);
    fe_sqr(&t, &t);
    fe_sub(&t, &t, &z1z1);
    fe_sqr(&t2, &h);
    fe_sub(&z3, &t, &t2);
    r->x = x3;
    r->y = y3;
    r->z = z3;
}

static void point_add_signed(jac_point *acc, const affine_point *q, int neg)
{
    affine_point t = *q;

    if (neg) {
        fe zero = { { 0 } };
        fe_sub(&t.y, &zero, &t.y);
    }
    point_add_mixed(acc, acc, &t);
}

/* Width-w non-adjacent form.  Returns the number of digits; naf[] holds
 * signed digits in [-2^(w-1)+1, 2^(w-1)-1]. */
static int scalar_wnaf(int8_t naf[CURVE_BITS + 2],
                       const uint64_t k[CURVE_LIMBS], int w)
{
    uint64_t t[CURVE_LIMBS + 1];
    int i;
    int len = 0;

    memcpy(t, k, sizeof(uint64_t) * CURVE_LIMBS);
    t[CURVE_LIMBS] = 0;
    for (;;) {
        int d = 0;
        int nonzero = 0;

        for (i = 0; i <= CURVE_LIMBS; i++) {
            nonzero |= t[i] != 0;
        }
        if (!nonzero) {
            break;
        }
        if (t[0] & 1u) {
            d = (int)(t[0] & ((1u << w) - 1u));
            if (d >= (1 << (w - 1))) {
                d -= (1 << w);
            }
            if (d > 0) {
                uint64_t borrow = (uint64_t)d;

                for (i = 0; i <= CURVE_LIMBS; i++) {
                    u128 uv = (u128)t[i] - borrow;

                    t[i] = (uint64_t)uv;
                    borrow = (uint64_t)(uv >> 64) & 1;
                    if (borrow == 0) {
                        break;
                    }
                }
            } else if (d < 0) {
                uint64_t carry = (uint64_t)(-d);

                for (i = 0; i <= CURVE_LIMBS; i++) {
                    u128 uv = (u128)t[i] + carry;

                    t[i] = (uint64_t)uv;
                    carry = (uint64_t)(uv >> 64);
                    if (carry == 0) {
                        break;
                    }
                }
            }
        }
        naf[len++] = (int8_t)d;
        for (i = 0; i < CURVE_LIMBS; i++) {
            t[i] = (t[i] >> 1) | (t[i + 1] << 63);
        }
        t[CURVE_LIMBS] >>= 1;
    }
    return len;
}

/* Batch-convert Jacobian points to affine with one field inversion. */
static void points_to_affine(affine_point out[], const jac_point in[], int n)
{
    fe prod[8];
    fe inv;
    fe zi2, zi3;
    int i;

    prod[0] = in[0].z;
    for (i = 1; i < n; i++) {
        fe_mul(&prod[i], &prod[i - 1], &in[i].z);
    }
    fe_inv_var(&inv, &prod[n - 1]);
    for (i = n - 1; i >= 0; i--) {
        fe zinv;

        if (i == 0) {
            zinv = inv;
        } else {
            fe_mul(&zinv, &inv, &prod[i - 1]);
            fe_mul(&inv, &inv, &in[i].z);
        }
        fe_sqr(&zi2, &zinv);
        fe_mul(&zi3, &zi2, &zinv);
        fe_mul(&out[i].x, &in[i].x, &zi2);
        fe_mul(&out[i].y, &in[i].y, &zi3);
    }
}

/* Odd multiples 1P, 3P, ..., 15P in affine form (width-5 wNAF). */
static void build_naf_table(affine_point table[8], const jac_point *p)
{
    jac_point jt[8];
    jac_point two;
    int i;

    point_double(&two, p);
    jt[0] = *p;
    for (i = 1; i < 8; i++) {
        point_add(&jt[i], &jt[i - 1], &two);
    }
    points_to_affine(table, jt, 8);
}

static uint64_t ct_eq_mask(uint64_t a, uint64_t b)
{
    uint64_t x = a ^ b;

    x = (x | (0 - x)) >> 63;
    return x - 1;
}

static void point_select(jac_point *r, const jac_point *table, unsigned idx)
{
    unsigned i;
    int j;

    point_set_infinity(r);
    for (i = 0; i < CURVE_TABLE_SIZE; i++) {
        uint64_t mask = ct_eq_mask(i, idx);

        for (j = 0; j < CURVE_LIMBS; j++) {
            r->x.v[j] |= table[i].x.v[j] & mask;
            r->y.v[j] |= table[i].y.v[j] & mask;
            r->z.v[j] |= table[i].z.v[j] & mask;
        }
    }
}

static void build_table(jac_point table[CURVE_TABLE_SIZE], const jac_point *p)
{
    int i;

    point_set_infinity(&table[0]);
    table[1] = *p;
    point_double(&table[2], p);
    for (i = 3; i < CURVE_TABLE_SIZE; i++) {
        point_add(&table[i], &table[i - 1], &table[1]);
    }
}

static unsigned scalar_window(const uint64_t v[CURVE_LIMBS], int idx)
{
#if CURVE_WINDOW_BITS == 4
    return (unsigned)((v[idx >> 4] >> ((idx & 15) * 4)) & 0xf);
#else
    int bit = idx * CURVE_WINDOW_BITS;
    int limb = bit >> 6;
    int shift = bit & 63;
    uint64_t x = v[limb] >> shift;

    if (shift + CURVE_WINDOW_BITS > 64 && limb + 1 < CURVE_LIMBS) {
        x |= v[limb + 1] << (64 - shift);
    }
    return (unsigned)(x & (CURVE_TABLE_SIZE - 1));
#endif
}

static void scalar_mul_table(jac_point *r, const jac_point *table,
                             const uint64_t k[CURVE_LIMBS])
{
    jac_point acc;
    jac_point t;
    int i;
    int j;

    (void)j;
    point_set_infinity(&acc);
    for (i = CURVE_WINDOWS - 1; i >= 0; i--) {
#if CURVE_WINDOW_BITS == 4
        point_double(&acc, &acc);
        point_double(&acc, &acc);
        point_double(&acc, &acc);
        point_double(&acc, &acc);
#else
        for (j = 0; j < CURVE_WINDOW_BITS; j++) {
            point_double(&acc, &acc);
        }
#endif
        point_select(&t, table, scalar_window(k, i));
        point_add(&acc, &acc, &t);
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

static fe g_one_mont;
static fe g_b_mont;
static jac_point g_table[CURVE_TABLE_SIZE];
static affine_point g_naf[8];
static int g_ready;

static void init_constants(void)
{
    CURVE_INIT_ENTER(g_ready);
    fe b;
    fe gx;
    fe gy;
    jac_point g;

    memcpy(g_one_mont.v, P_R, sizeof(g_one_mont.v));
    memcpy(b.v, B, sizeof(b.v));
    fe_to_mont(&g_b_mont, &b);
    memcpy(gx.v, GX, sizeof(gx.v));
    memcpy(gy.v, GY, sizeof(gy.v));
    fe_to_mont(&g.x, &gx);
    fe_to_mont(&g.y, &gy);
    g.z = g_one_mont;
    build_table(g_table, &g);
    build_naf_table(g_naf, &g);
    CURVE_INIT_LEAVE(g_ready);
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

static int parse_point(jac_point *p, const uint8_t qx[CURVE_BYTES],
                       const uint8_t qy[CURVE_BYTES])
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

int FN(public_key)(uint8_t qx[CURVE_BYTES], uint8_t qy[CURVE_BYTES],
                   const uint8_t d[CURVE_BYTES])
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
    u_to_bytes(qx, &x);
    u_to_bytes(qy, &y);
    return 0;
}

int FN(validate_public_key)(const uint8_t qx[CURVE_BYTES],
                            const uint8_t qy[CURVE_BYTES])
{
    jac_point q;
    jac_point table[CURVE_TABLE_SIZE];
    jac_point r;

    init_constants();
    if (parse_point(&q, qx, qy) != 0) {
        return 0;
    }
    build_table(table, &q);
    scalar_mul_table(&r, table, N);
    return point_is_infinity(&r) ? 1 : 0;
}

int FN(ecdsa_sign)(uint8_t r_out[CURVE_BYTES], uint8_t s_out[CURVE_BYTES],
                   const uint8_t d[CURVE_BYTES], const uint8_t z[CURVE_BYTES],
                   const uint8_t k[CURVE_BYTES])
{
    sc dv, kv, zv;
    sc km, dm, zm, rm;
    sc rd, t, kinv, s, r;
    jac_point kg;
    fe x, y;
    fe xc;

    if (sc_parse(&dv, d, 1) != 0 || sc_parse(&kv, k, 1) != 0 ||
        sc_from_z(&zv, z) != 0) {
        return -1;
    }
    init_constants();
    scalar_mul_table(&kg, g_table, kv.v);
    if (point_is_infinity(&kg)) {
        return -1;
    }
    point_to_affine(&x, &y, &kg);
    fe_from_mont(&xc, &x);
    r = xc;
    sc_reduce_once(&r);
    if (u_is_zero(&r)) {
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
    if (u_is_zero(&s)) {
        return -1;
    }
    u_to_bytes(r_out, &r);
    u_to_bytes(s_out, &s);
    return 0;
}

int FN(ecdsa_verify)(const uint8_t qx[CURVE_BYTES],
                     const uint8_t qy[CURVE_BYTES],
                     const uint8_t z[CURVE_BYTES],
                     const uint8_t r_in[CURVE_BYTES],
                     const uint8_t s_in[CURVE_BYTES])
{
    sc r, s, zv, zm, sm, rm, w, u1, u2;
    jac_point q;
    jac_point xp;

    if (sc_parse(&r, r_in, 1) != 0 || sc_parse(&s, s_in, 1) != 0 ||
        sc_from_z(&zv, z) != 0) {
        return 0;
    }
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

    {
        int8_t naf1[CURVE_BITS + 2];
        int8_t naf2[CURVE_BITS + 2];
        int len1;
        int len2;
        int len;
        int i;
        affine_point qtab[8];

        len1 = scalar_wnaf(naf1, u1.v, 5);
        len2 = scalar_wnaf(naf2, u2.v, 5);
        build_naf_table(qtab, &q);
        len = len1 > len2 ? len1 : len2;
        point_set_infinity(&xp);
        for (i = len - 1; i >= 0; i--) {
            if (!point_is_infinity(&xp)) {
                point_double(&xp, &xp);
            }
            if (i < len1 && naf1[i] != 0) {
                int d = naf1[i];
                int idx = ((d < 0 ? -d : d) - 1) / 2;

                point_add_signed(&xp, &g_naf[idx], d < 0);
            }
            if (i < len2 && naf2[i] != 0) {
                int d = naf2[i];
                int idx = ((d < 0 ? -d : d) - 1) / 2;

                point_add_signed(&xp, &qtab[idx], d < 0);
            }
        }
    }
    if (point_is_infinity(&xp)) {
        return 0;
    }
    {
        fe z2;
        fe rhs;
        fe pf;
        sc rn;
        uint64_t carry = 0;
        int i;

        /* x_affine = X/Z^2.  Avoid the inversion: compare X with r*Z^2 in
         * Montgomery form.  Since p < 2n, x mod n is either x or x-n, so at
         * most one extra candidate r+n needs checking. */
        fe rfe;

        fe_sqr(&z2, &xp.z);
        memcpy(rfe.v, r.v, sizeof(rfe.v));
        fe_to_mont(&rfe, &rfe);
        fe_mul(&rhs, &rfe, &z2);
        if (fe_is_equal(&xp.x, &rhs)) {
            return 1;
        }
        for (i = 0; i < CURVE_LIMBS; i++) {
            u128 uv = (u128)r.v[i] + N[i] + carry;

            rn.v[i] = (uint64_t)uv;
            carry = (uint64_t)(uv >> 64);
        }
        memcpy(pf.v, P, sizeof(pf.v));
        if (carry != 0 || fe_cmp_var(&rn, &pf) >= 0) {
            return 0;
        }
        fe_to_mont(&rn, &rn);
        fe_mul(&rhs, &rn, &z2);
        return fe_is_equal(&xp.x, &rhs) ? 1 : 0;
    }
}

int FN(ecdh)(uint8_t out[CURVE_BYTES], const uint8_t d[CURVE_BYTES],
             const uint8_t qx[CURVE_BYTES], const uint8_t qy[CURVE_BYTES])
{
    sc dv;
    jac_point q;
    jac_point table[CURVE_TABLE_SIZE];
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
    u_to_bytes(out, &x);
    return 0;
}
