//! smcPrime - Ultra-fast primality testing
//!
//! Deterministic primality testing for 32-bit and 64-bit integers using
//! Montgomery arithmetic and optimal Miller-Rabin witnesses.
//!
//! # Example
//! ```
//! use smcprime::{is_prime, is_prime32, is_prime64, next_prime, prev_prime};
//!
//! assert!(is_prime(17));
//! assert!(!is_prime(18));
//! assert_eq!(next_prime(100), 101);
//! assert_eq!(prev_prime(100), 97);
//! ```

#![no_std]

#[cfg(feature = "std")]
extern crate std;

// ============================================================================
// 32-BIT PRIMALITY TESTING
// ============================================================================

/// Modular multiplication for 32-bit integers
#[inline(always)]
fn mulmod32(a: u32, b: u32, m: u32) -> u32 {
    ((a as u64 * b as u64) % m as u64) as u32
}

/// Modular exponentiation for 32-bit integers
#[inline(always)]
fn powmod32(mut a: u32, mut b: u32, m: u32) -> u32 {
    let mut r = 1u32;
    a %= m;
    while b > 0 {
        if b & 1 == 1 {
            r = mulmod32(r, a, m);
        }
        b >>= 1;
        if b > 0 {
            a = mulmod32(a, a, m);
        }
    }
    r
}

/// Strong probable prime test for 32-bit integers
#[inline(always)]
fn sprp32(n: u32, a: u32) -> bool {
    if a % n == 0 {
        return true;
    }
    let mut d = n - 1;
    let mut s = 0u32;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }
    let mut x = powmod32(a, d, n);
    if x == 1 || x == n - 1 {
        return true;
    }
    for _ in 1..s {
        x = mulmod32(x, x, n);
        if x == n - 1 {
            return true;
        }
        if x == 1 {
            return false;
        }
    }
    false
}

/// Deterministic primality test for 32-bit integers
///
/// Uses witnesses {2, 7, 61} which are sufficient for all n < 2^32.
/// Special case: 3215031751 is the only pseudoprime to all three bases.
#[inline]
pub fn is_prime32(n: u32) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n & 1 == 0 {
        return false;
    }
    if n < 9 {
        return true;
    }
    if n % 3 == 0 || n % 5 == 0 || n % 7 == 0 {
        return false;
    }
    if n == 3215031751 {
        return false;
    }
    sprp32(n, 2) && sprp32(n, 7) && sprp32(n, 61)
}

/// Find the next prime >= n (32-bit)
#[inline]
pub fn next_prime32(n: u32) -> u32 {
    if n <= 2 {
        return 2;
    }
    let mut n = if n & 1 == 0 { n + 1 } else { n };
    while !is_prime32(n) {
        n = n.wrapping_add(2);
        if n < 2 {
            return 0;
        }
    }
    n
}

/// Find the previous prime <= n (32-bit)
#[inline]
pub fn prev_prime32(n: u32) -> u32 {
    if n < 2 {
        return 0;
    }
    if n == 2 {
        return 2;
    }
    let mut n = if n & 1 == 0 { n - 1 } else { n };
    while !is_prime32(n) {
        if n < 3 {
            return 2;
        }
        n -= 2;
    }
    n
}

// ============================================================================
// 64-BIT PRIMALITY TESTING (Montgomery Arithmetic)
// ============================================================================

/// Montgomery inverse via Hensel lifting
#[inline(always)]
fn mont_inv64(n: u64) -> u64 {
    let mut est = (3u64.wrapping_mul(n)) ^ 2;
    est = (2u64.wrapping_sub(est.wrapping_mul(n))).wrapping_mul(est);
    est = (2u64.wrapping_sub(est.wrapping_mul(n))).wrapping_mul(est);
    est = (2u64.wrapping_sub(est.wrapping_mul(n))).wrapping_mul(est);
    est = (2u64.wrapping_sub(est.wrapping_mul(n))).wrapping_mul(est);
    est
}

/// Montgomery reduction
#[inline(always)]
fn mont_reduce64(x_lo: u64, x_hi: u64, n: u64, n_inv: u64) -> u64 {
    let m = x_lo.wrapping_mul(n_inv);
    let t = ((m as u128 * n as u128) >> 64) as u64;
    if x_hi < t {
        x_hi.wrapping_sub(t).wrapping_add(n)
    } else {
        x_hi - t
    }
}

/// Montgomery multiplication
#[inline(always)]
fn mont_mul64(a: u64, b: u64, n: u64, n_inv: u64) -> u64 {
    let prod = a as u128 * b as u128;
    mont_reduce64(prod as u64, (prod >> 64) as u64, n, n_inv)
}

/// Convert to Montgomery form
#[inline(always)]
fn to_mont64(x: u64, n: u64) -> u64 {
    ((x as u128) << 64 % n as u128) as u64
}

/// One in Montgomery form
#[inline(always)]
fn mont_one64(n: u64) -> u64 {
    (u64::MAX % n).wrapping_add(1)
}

/// Montgomery exponentiation
#[inline(always)]
fn mont_pow64(mut base: u64, mut exp: u64, n: u64, n_inv: u64, one: u64) -> u64 {
    let mut result = one;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mont_mul64(result, base, n, n_inv);
        }
        base = mont_mul64(base, base, n, n_inv);
        exp >>= 1;
    }
    result
}

/// Strong Fermat test in Montgomery form
#[inline(always)]
fn mont_sprp64(n: u64, a: u64, n_inv: u64, one: u64) -> bool {
    let mut d = n - 1;
    let mut s = 0u32;
    while d & 1 == 0 {
        d >>= 1;
        s += 1;
    }

    let a_mont = to_mont64(a % n, n);
    if a_mont == 0 {
        return true;
    }

    let x = mont_pow64(a_mont, d, n, n_inv, one);
    let neg_one = if n > one { n - one } else { n.wrapping_sub(one) };

    if x == one || x == neg_one {
        return true;
    }

    let mut x = x;
    for _ in 1..s {
        x = mont_mul64(x, x, n, n_inv);
        if x == neg_one {
            return true;
        }
        if x == one {
            return false;
        }
    }
    false
}

/// Deterministic primality test for 64-bit integers
///
/// Uses Montgomery-based Miller-Rabin with optimal witnesses.
#[inline]
pub fn is_prime64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 {
        return true;
    }
    if n & 1 == 0 {
        return false;
    }
    if n < 9 {
        return true;
    }
    if n % 3 == 0 || n % 5 == 0 || n % 7 == 0 {
        return false;
    }
    if n == 3215031751 {
        return false;
    }

    // For small n, use 32-bit version
    if n <= u32::MAX as u64 {
        return is_prime32(n as u32);
    }

    let n_inv = mont_inv64(n);
    let one = mont_one64(n);

    // Witnesses sufficient for all 64-bit integers
    const WITNESSES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &w in &WITNESSES {
        if !mont_sprp64(n, w, n_inv, one) {
            return false;
        }
    }
    true
}

/// Find the next prime >= n (64-bit)
#[inline]
pub fn next_prime64(n: u64) -> u64 {
    if n <= 2 {
        return 2;
    }
    let mut n = if n & 1 == 0 { n + 1 } else { n };
    while !is_prime64(n) {
        n = n.wrapping_add(2);
        if n < 2 {
            return 0;
        }
    }
    n
}

/// Find the previous prime <= n (64-bit)
#[inline]
pub fn prev_prime64(n: u64) -> u64 {
    if n < 2 {
        return 0;
    }
    if n == 2 {
        return 2;
    }
    let mut n = if n & 1 == 0 { n - 1 } else { n };
    while !is_prime64(n) {
        if n < 3 {
            return 2;
        }
        n -= 2;
    }
    n
}

// ============================================================================
// DEFAULT ALIASES (64-bit)
// ============================================================================

/// Primality test (defaults to 64-bit)
#[inline]
pub fn is_prime(n: u64) -> bool {
    is_prime64(n)
}

/// Find next prime (defaults to 64-bit)
#[inline]
pub fn next_prime(n: u64) -> u64 {
    next_prime64(n)
}

/// Find previous prime (defaults to 64-bit)
#[inline]
pub fn prev_prime(n: u64) -> u64 {
    prev_prime64(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_primes() {
        let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        for &p in &primes {
            assert!(is_prime(p), "{} should be prime", p);
        }
    }

    #[test]
    fn test_small_composites() {
        let composites = [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25];
        for &c in &composites {
            assert!(!is_prime(c), "{} should not be prime", c);
        }
    }

    #[test]
    fn test_special_pseudoprime() {
        assert!(!is_prime32(3215031751));
        assert!(!is_prime64(3215031751));
    }

    #[test]
    fn test_next_prime() {
        assert_eq!(next_prime(0), 2);
        assert_eq!(next_prime(2), 2);
        assert_eq!(next_prime(3), 3);
        assert_eq!(next_prime(4), 5);
        assert_eq!(next_prime(100), 101);
    }

    #[test]
    fn test_prev_prime() {
        assert_eq!(prev_prime(2), 2);
        assert_eq!(prev_prime(3), 3);
        assert_eq!(prev_prime(4), 3);
        assert_eq!(prev_prime(100), 97);
    }

    #[test]
    fn test_large_primes() {
        // Some known large primes
        assert!(is_prime64(1000000007));
        assert!(is_prime64(1000000009));
        assert!(is_prime64(18446744073709551557)); // Largest 64-bit prime
    }
}
