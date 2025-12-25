/*
 * smcPrime - Ultra-Fast Primality Testing Library (32-bit and 64-bit)
 * 
 * A unified primality testing library combining:
 * - smcPrime32: Fast 32-bit primality using native 64-bit math
 * - smcPrime64: Ultra-fast 64-bit primality using Montgomery arithmetic
 * 
 * Adopts optimizations from machine-prime (fastest Rust prime library):
 * - Montgomery arithmetic for Miller-Rabin (avoids modular division)
 * - Prime-inverse trial division (multiplication instead of %)
 * - Hensel lifting for modular inverses (faster than extended GCD)
 * - Optimal witness selection (minimal Miller-Rabin rounds)
 * 
 * Copyright 2025 ScaleCode Solutions
 * Released under MIT License
 */

#ifndef SMCPRIME_H
#define SMCPRIME_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef SMC_INLINE
  #if defined(_MSC_VER)
    #define SMC_INLINE static __forceinline
  #elif defined(__GNUC__) || defined(__clang__)
    #define SMC_INLINE static inline __attribute__((always_inline))
  #else
    #define SMC_INLINE static inline
  #endif
#endif

/* ===========================================================================
 * 32-BIT PRIMALITY TESTING
 * 
 * Uses native 64-bit math for fast modular arithmetic.
 * Only 3 witnesses needed: {2, 7, 61} (Jaeschke 1993)
 * =========================================================================== */

SMC_INLINE uint32_t smc_mulmod32(uint32_t a, uint32_t b, uint32_t m) {
    return (uint32_t)(((uint64_t)a * b) % m);
}

SMC_INLINE uint32_t smc_powmod32(uint32_t a, uint32_t b, uint32_t m) {
    uint32_t r = 1;
    a %= m;
    while (b) {
        if (b & 1) r = smc_mulmod32(r, a, m);
        b >>= 1;
        if (b) a = smc_mulmod32(a, a, m);
    }
    return r;
}

SMC_INLINE bool smc_sprp32(uint32_t n, uint32_t a) {
    if (a % n == 0) return true;
    uint32_t d = n - 1;
    uint32_t s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    uint32_t x = smc_powmod32(a, d, n);
    if (x == 1 || x == n - 1) return true;
    for (uint32_t r = 1; r < s; r++) {
        x = smc_mulmod32(x, x, n);
        if (x == n - 1) return true;
        if (x == 1) return false;
    }
    return false;
}

/*
 * Deterministic primality test for 32-bit integers
 * 
 * Uses witnesses {2, 7, 61} which are sufficient for all n < 2^32
 * Special case: 3215031751 is the only pseudoprime to all three bases
 */
SMC_INLINE bool smc_is_prime32(uint32_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;
    if (n < 9) return true;
    if (n % 3 == 0) return false;
    if (n % 5 == 0) return false;
    if (n % 7 == 0) return false;
    
    if (n == 3215031751U) return false;  /* Special pseudoprime */
    
    return smc_sprp32(n, 2) && smc_sprp32(n, 7) && smc_sprp32(n, 61);
}

SMC_INLINE uint32_t smc_next_prime32(uint32_t n) {
    if (n <= 2) return 2;
    if (n == 3) return 3;
    if ((n & 1) == 0) n++;
    while (!smc_is_prime32(n)) { n += 2; if (n < 2) return 0; }
    return n;
}

SMC_INLINE uint32_t smc_prev_prime32(uint32_t n) {
    if (n < 2) return 0;
    if (n == 2) return 2;
    if ((n & 1) == 0) n--;
    while (!smc_is_prime32(n)) { if (n < 3) return 2; n -= 2; }
    return n;
}

/* ===========================================================================
 * 64-BIT PRIMALITY TESTING (Montgomery Arithmetic)
 * 
 * Uses Montgomery multiplication to avoid expensive modular division.
 * Prime-inverse trial division for fast composite rejection.
 * 
 * ALERT from machine-prime: When using witness tables, ensure the strong
 * Fermat test checks that N is not a non-zero multiple of the witness.
 * This is handled by: witness = witness % N (or Montgomery transform)
 * =========================================================================== */

/* Prime inverses mod 2^64 for trial division (primes 3 to 331) */
static const uint64_t SMC_PRIME_INV64[] = {
    0xAAAAAAAAAAAAAAABULL, 0xCCCCCCCCCCCCCCCDULL, 0x6DB6DB6DB6DB6DB7ULL, 0x2E8BA2E8BA2E8BA3ULL,
    0x4EC4EC4EC4EC4EC5ULL, 0xF0F0F0F0F0F0F0F1ULL, 0x86BCA1AF286BCA1BULL, 0xD37A6F4DE9BD37A7ULL,
    0x34F72C234F72C235ULL, 0xEF7BDEF7BDEF7BDFULL, 0x14C1BACF914C1BADULL, 0x8F9C18F9C18F9C19ULL,
    0x82FA0BE82FA0BE83ULL, 0x51B3BEA3677D46CFULL, 0x21CFB2B78C13521DULL, 0xCBEEA4E1A08AD8F3ULL,
    0x4FBCDA3AC10C9715ULL, 0xF0B7672A07A44C6BULL, 0x193D4BB7E327A977ULL, 0x7E3F1F8FC7E3F1F9ULL,
    0x9B8B577E613716AFULL, 0xA3784A062B2E43DBULL, 0xF47E8FD1FA3F47E9ULL, 0xA3A0FD5C5F02A3A1ULL,
    0x3A4C0A237C32B16DULL, 0xDAB7EC1DD3431B57ULL, 0x77A04C8F8D28AC43ULL, 0xA6C0964FDA6C0965ULL,
    0x90FDBC090FDBC091ULL, 0x7EFDFBF7EFDFBF7FULL, 0x03E88CB3C9484E2BULL, 0xE21A291C077975B9ULL,
    0x3AEF6CA970586723ULL, 0xDF5B0F768CE2CABDULL, 0x6FE4DFC9BF937F27ULL, 0x5B4FE5E92C0685B5ULL,
    0x1F693A1C451AB30BULL, 0x8D07AA27DB35A717ULL, 0x882383B30D516325ULL, 0xED6866F8D962AE7BULL,
    0x3454DCA410F8ED9DULL, 0x1D7CA632EE936F3FULL, 0x70BF015390948F41ULL, 0xC96BDB9D3D137E0DULL,
    0x2697CC8AEF46C0F7ULL, 0xC0E8F2A76E68575BULL, 0x687763DFDB43BB1FULL, 0x1B10EA929BA144CBULL,
    0x1D10C4C0478BBCEDULL, 0x63FB9AEB1FDCD759ULL, 0x64AFAA4F437B2E0FULL, 0xF010FEF010FEF011ULL,
    0x28CBFBEB9A020A33ULL, 0xFF00FF00FF00FF01ULL, 0xD624FD1470E99CB7ULL, 0x8FB3DDBD6205B5C5ULL,
    0xD57DA36CA27ACDEFULL, 0xEE70C03B25E4463DULL, 0xC5B1A6B80749CB29ULL, 0x47768073C9B97113ULL,
    0x2591E94884CE32ADULL, 0xF02806ABC74BE1FBULL, 0x7EC3E8F3A7198487ULL, 0x58550F8A39409D09ULL,
    0xEC9E48AE6F71DE15ULL, 0x2FF3A018BFCE8063ULL,
};
#define SMC_NUM_PRIME_INV64 66

/* Montgomery inverse via Hensel lifting (Newton-Raphson iteration) */
SMC_INLINE uint64_t smc_mont_inv64(uint64_t n) {
    uint64_t est = (3 * n) ^ 2;
    est = (2 - est * n) * est;
    est = (2 - est * n) * est;
    est = (2 - est * n) * est;
    est = (2 - est * n) * est;
    return est;
}

/* Montgomery reduction */
SMC_INLINE uint64_t smc_mont_reduce64(uint64_t x_lo, uint64_t x_hi, uint64_t n, uint64_t n_inv) {
    uint64_t m = x_lo * n_inv;
#if defined(__SIZEOF_INT128__)
    uint64_t t = (uint64_t)(((__uint128_t)m * n) >> 64);
#else
    uint64_t a_lo = (uint32_t)m, a_hi = m >> 32;
    uint64_t b_lo = (uint32_t)n, b_hi = n >> 32;
    uint64_t p0 = a_lo * b_lo, p1 = a_lo * b_hi, p2 = a_hi * b_lo, p3 = a_hi * b_hi;
    uint64_t cy = ((p0 >> 32) + (uint32_t)p1 + (uint32_t)p2) >> 32;
    uint64_t t = p3 + (p1 >> 32) + (p2 >> 32) + cy;
#endif
    return (x_hi < t) ? x_hi - t + n : x_hi - t;
}

/* Montgomery multiplication */
SMC_INLINE uint64_t smc_mont_mul64(uint64_t a, uint64_t b, uint64_t n, uint64_t n_inv) {
#if defined(__SIZEOF_INT128__)
    __uint128_t prod = (__uint128_t)a * b;
    return smc_mont_reduce64((uint64_t)prod, (uint64_t)(prod >> 64), n, n_inv);
#else
    uint64_t a_lo = (uint32_t)a, a_hi = a >> 32;
    uint64_t b_lo = (uint32_t)b, b_hi = b >> 32;
    uint64_t p0 = a_lo * b_lo, p1 = a_lo * b_hi, p2 = a_hi * b_lo, p3 = a_hi * b_hi;
    uint64_t lo = p0 + (p1 << 32) + (p2 << 32);
    uint64_t cy = (lo < p0) + ((p0 >> 32) + (uint32_t)p1 + (uint32_t)p2 >= 0x100000000ULL);
    uint64_t hi = p3 + (p1 >> 32) + (p2 >> 32) + cy;
    return smc_mont_reduce64(lo, hi, n, n_inv);
#endif
}

/* Convert to Montgomery form */
SMC_INLINE uint64_t smc_to_mont64(uint64_t x, uint64_t n) {
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((__uint128_t)x << 64) % n);
#else
    uint64_t r = x % n;
    for (int i = 0; i < 64; i++) { r <<= 1; if (r >= n) r -= n; }
    return r;
#endif
}

/* One in Montgomery form */
SMC_INLINE uint64_t smc_mont_one64(uint64_t n) {
    return (UINT64_MAX % n) + 1;
}

/* Montgomery exponentiation */
SMC_INLINE uint64_t smc_mont_pow64(uint64_t base, uint64_t exp, uint64_t n, uint64_t n_inv, uint64_t one) {
    uint64_t result = one;
    while (exp > 0) {
        if (exp & 1) result = smc_mont_mul64(result, base, n, n_inv);
        base = smc_mont_mul64(base, base, n, n_inv);
        exp >>= 1;
    }
    return result;
}

/*
 * Strong Fermat test in Montgomery form
 * 
 * IMPORTANT: The witness 'a' is taken mod n via the Montgomery transform
 * (smc_to_mont64 does a % n internally). This ensures we don't have issues
 * when n is a multiple of the witness, as warned in machine-prime.
 */
SMC_INLINE bool smc_mont_sprp64(uint64_t n, uint64_t a, uint64_t n_inv, uint64_t one) {
    uint64_t d = n - 1;
    uint32_t s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    
    /* Convert witness to Montgomery form (handles a % n automatically) */
    uint64_t a_mont = smc_to_mont64(a % n, n);
    if (a_mont == 0) return true;  /* a is multiple of n, trivially passes */
    
    uint64_t x = smc_mont_pow64(a_mont, d, n, n_inv, one);
    uint64_t neg_one = n - one;
    if (neg_one >= n) neg_one -= n;
    
    if (x == one || x == neg_one) return true;
    
    for (uint32_t r = 1; r < s; r++) {
        x = smc_mont_mul64(x, x, n, n_inv);
        if (x == neg_one) return true;
        if (x == one) return false;
    }
    return false;
}

/*
 * Deterministic primality test for 64-bit integers
 * 
 * Uses prime-inverse trial division for fast composite rejection,
 * then Montgomery-based Miller-Rabin with optimal witnesses.
 */
SMC_INLINE bool smc_is_prime64(uint64_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;
    if (n < 9) return true;
    
    /* 
     * Fast trial division using prime inverses (from machine-prime)
     * Only valid for n < 55730344633563600 to avoid overflow issues
     */
    if (n < 55730344633563600ULL) {
        for (size_t i = 0; i < SMC_NUM_PRIME_INV64; i++) {
            uint64_t prod = n * SMC_PRIME_INV64[i];
            if (prod == 1) return true;   /* n IS this prime */
            if (prod < n) return false;   /* n is divisible by this prime */
        }
        if (n < 109561) return true;  /* 331^2, passed all trial divisions */
    } else {
        /* For large n, use traditional trial division */
        if (n % 3 == 0) return false;
        if (n % 5 == 0) return false;
        if (n % 7 == 0) return false;
        if (n % 11 == 0) return false;
        if (n % 13 == 0) return false;
    }
    
    /* Montgomery setup */
    uint64_t n_inv = smc_mont_inv64(n);
    uint64_t one = smc_mont_one64(n);
    
    /* Miller-Rabin with deterministic witnesses for 64-bit */
    if (!smc_mont_sprp64(n, 2, n_inv, one)) return false;
    if (n < 2047) return true;
    
    if (!smc_mont_sprp64(n, 3, n_inv, one)) return false;
    if (n < 1373653) return true;
    
    if (!smc_mont_sprp64(n, 5, n_inv, one)) return false;
    if (n < 25326001) return true;
    
    if (!smc_mont_sprp64(n, 7, n_inv, one)) return false;
    if (n < 3215031751ULL) return true;
    if (n == 3215031751ULL) return false;  /* Special pseudoprime */
    
    if (!smc_mont_sprp64(n, 11, n_inv, one)) return false;
    if (n < 2152302898747ULL) return true;
    
    if (!smc_mont_sprp64(n, 13, n_inv, one)) return false;
    if (n < 3474749660383ULL) return true;
    
    if (!smc_mont_sprp64(n, 17, n_inv, one)) return false;
    if (n < 341550071728321ULL) return true;
    
    if (!smc_mont_sprp64(n, 19, n_inv, one)) return false;
    if (!smc_mont_sprp64(n, 23, n_inv, one)) return false;
    if (n < 3825123056546413051ULL) return true;
    
    if (!smc_mont_sprp64(n, 29, n_inv, one)) return false;
    if (!smc_mont_sprp64(n, 31, n_inv, one)) return false;
    if (!smc_mont_sprp64(n, 37, n_inv, one)) return false;
    
    return true;
}

/*
 * Worst-case optimized version (for numbers likely to be prime)
 * Skips trial division, goes straight to Miller-Rabin
 */
SMC_INLINE bool smc_is_prime64_wc(uint64_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;
    if (n < 9) return true;
    if (n == 3215031751ULL) return false;
    
    uint64_t n_inv = smc_mont_inv64(n);
    uint64_t one = smc_mont_one64(n);
    
    static const uint8_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) {
        if (!smc_mont_sprp64(n, witnesses[i], n_inv, one)) return false;
    }
    return true;
}

SMC_INLINE uint64_t smc_next_prime64(uint64_t n) {
    if (n <= 2) return 2;
    if (n == 3) return 3;
    if ((n & 1) == 0) n++;
    while (!smc_is_prime64(n)) { n += 2; if (n < 2) return 0; }
    return n;
}

SMC_INLINE uint64_t smc_prev_prime64(uint64_t n) {
    if (n < 2) return 0;
    if (n == 2) return 2;
    if ((n & 1) == 0) n--;
    while (!smc_is_prime64(n)) { if (n < 3) return 2; n -= 2; }
    return n;
}

/* ===========================================================================
 * CONVENIENCE ALIASES
 * =========================================================================== */

/* Default to 64-bit for smc_is_prime */
#define smc_is_prime     smc_is_prime64
#define smc_is_prime_wc  smc_is_prime64_wc
#define smc_next_prime   smc_next_prime64
#define smc_prev_prime   smc_prev_prime64

#ifdef __cplusplus
}
#endif

#endif /* SMCPRIME_H */
