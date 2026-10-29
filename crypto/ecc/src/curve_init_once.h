/*
 * Copyright (c) 2026 harmor.
 *
 * One-time, thread-safe initialization helpers for the project ECC tables.
 * Licensed under the Mulan PSL v2; see LICENSE in the project root.
 *
 * The state is 0 (uninitialized), 1 (initializing) or 2 (ready).  A losing
 * thread waits for the winner to publish the tables before returning.
 */
#ifndef CURVE_INIT_ONCE_H
#define CURVE_INIT_ONCE_H

#if defined(__GNUC__) || defined(__clang__)
#define CURVE_INIT_ENTER(state)                                                     \
    do {                                                                            \
        int curve_expected = 0;                                                     \
        if ((state) == 2) {                                                         \
            return;                                                                 \
        }                                                                           \
        if (!__atomic_compare_exchange_n(&(state), &curve_expected, 1, 0,           \
                __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {                              \
            while (__atomic_load_n(&(state), __ATOMIC_ACQUIRE) != 2) {              \
            }                                                                       \
            return;                                                                 \
        }                                                                           \
    } while (0)

#define CURVE_INIT_LEAVE(state) __atomic_store_n(&(state), 2, __ATOMIC_RELEASE)
#else
#define CURVE_INIT_ENTER(state) do { if ((state) != 0) return; (state) = 1; } while (0)
#define CURVE_INIT_LEAVE(state) do { (state) = 2; } while (0)
#endif

#endif /* CURVE_INIT_ONCE_H */
