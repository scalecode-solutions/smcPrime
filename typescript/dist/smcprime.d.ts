/**
 * smcPrime - Ultra-fast primality testing
 *
 * Deterministic primality testing for 32-bit and 64-bit integers using
 * Montgomery arithmetic and optimal Miller-Rabin witnesses.
 */
/**
 * Deterministic primality test for 32-bit integers
 * Uses witnesses {2, 7, 61} which are sufficient for all n < 2^32
 */
export declare function isPrime32(n: number | bigint): boolean;
/** Find the next prime >= n (32-bit) */
export declare function nextPrime32(n: number | bigint): bigint;
/** Find the previous prime <= n (32-bit) */
export declare function prevPrime32(n: number | bigint): bigint;
/**
 * Deterministic primality test for 64-bit integers
 * Uses Montgomery-based Miller-Rabin with optimal witnesses
 */
export declare function isPrime64(n: number | bigint): boolean;
/** Find the next prime >= n (64-bit) */
export declare function nextPrime64(n: number | bigint): bigint;
/** Find the previous prime <= n (64-bit) */
export declare function prevPrime64(n: number | bigint): bigint;
/** Primality test (defaults to 64-bit) */
export declare function isPrime(n: number | bigint): boolean;
/** Find next prime (defaults to 64-bit) */
export declare function nextPrime(n: number | bigint): bigint;
/** Find previous prime (defaults to 64-bit) */
export declare function prevPrime(n: number | bigint): bigint;
