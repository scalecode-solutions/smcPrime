/**
 * smcPrime - Ultra-fast primality testing
 * 
 * Deterministic primality testing for 32-bit and 64-bit integers using
 * Montgomery arithmetic and optimal Miller-Rabin witnesses.
 */

const MASK32 = 0xFFFFFFFFn;
const MASK64 = 0xFFFFFFFFFFFFFFFFn;

// ============================================================================
// 32-BIT PRIMALITY TESTING
// ============================================================================

function mulmod32(a: bigint, b: bigint, m: bigint): bigint {
    return (a * b) % m;
}

function powmod32(a: bigint, b: bigint, m: bigint): bigint {
    let r = 1n;
    a = a % m;
    while (b > 0n) {
        if (b & 1n) r = mulmod32(r, a, m);
        b >>= 1n;
        if (b > 0n) a = mulmod32(a, a, m);
    }
    return r;
}

function sprp32(n: bigint, a: bigint): boolean {
    if (a % n === 0n) return true;
    let d = n - 1n;
    let s = 0n;
    while ((d & 1n) === 0n) { d >>= 1n; s++; }
    let x = powmod32(a, d, n);
    if (x === 1n || x === n - 1n) return true;
    for (let r = 1n; r < s; r++) {
        x = mulmod32(x, x, n);
        if (x === n - 1n) return true;
        if (x === 1n) return false;
    }
    return false;
}

/**
 * Deterministic primality test for 32-bit integers
 * Uses witnesses {2, 7, 61} which are sufficient for all n < 2^32
 */
export function isPrime32(n: number | bigint): boolean {
    const bn = BigInt(n);
    if (bn < 2n) return false;
    if (bn === 2n) return true;
    if ((bn & 1n) === 0n) return false;
    if (bn < 9n) return true;
    if (bn % 3n === 0n || bn % 5n === 0n || bn % 7n === 0n) return false;
    if (bn === 3215031751n) return false;
    return sprp32(bn, 2n) && sprp32(bn, 7n) && sprp32(bn, 61n);
}

/** Find the next prime >= n (32-bit) */
export function nextPrime32(n: number | bigint): bigint {
    let bn = BigInt(n);
    if (bn <= 2n) return 2n;
    if ((bn & 1n) === 0n) bn++;
    while (!isPrime32(bn)) {
        bn += 2n;
        if (bn > MASK32) return 0n;
    }
    return bn;
}

/** Find the previous prime <= n (32-bit) */
export function prevPrime32(n: number | bigint): bigint {
    let bn = BigInt(n);
    if (bn < 2n) return 0n;
    if (bn === 2n) return 2n;
    if ((bn & 1n) === 0n) bn--;
    while (!isPrime32(bn)) {
        if (bn < 3n) return 2n;
        bn -= 2n;
    }
    return bn;
}

// ============================================================================
// 64-BIT PRIMALITY TESTING (Montgomery Arithmetic)
// ============================================================================

function montInv64(n: bigint): bigint {
    let est = ((3n * n) ^ 2n) & MASK64;
    est = ((2n - est * n) * est) & MASK64;
    est = ((2n - est * n) * est) & MASK64;
    est = ((2n - est * n) * est) & MASK64;
    est = ((2n - est * n) * est) & MASK64;
    return est;
}

function montReduce64(xLo: bigint, xHi: bigint, n: bigint, nInv: bigint): bigint {
    const m = (xLo * nInv) & MASK64;
    const t = (m * n) >> 64n;
    if (xHi < t) {
        return ((xHi - t) & MASK64) + n;
    }
    return xHi - t;
}

function montMul64(a: bigint, b: bigint, n: bigint, nInv: bigint): bigint {
    const prod = a * b;
    return montReduce64(prod & MASK64, prod >> 64n, n, nInv);
}

function toMont64(x: bigint, n: bigint): bigint {
    return ((x % n) << 64n) % n;
}

function montOne64(n: bigint): bigint {
    return (MASK64 % n) + 1n;
}

function montPow64(base: bigint, exp: bigint, n: bigint, nInv: bigint, one: bigint): bigint {
    let result = one;
    while (exp > 0n) {
        if (exp & 1n) result = montMul64(result, base, n, nInv);
        base = montMul64(base, base, n, nInv);
        exp >>= 1n;
    }
    return result;
}

function montSprp64(n: bigint, a: bigint, nInv: bigint, one: bigint): boolean {
    let d = n - 1n;
    let s = 0n;
    while ((d & 1n) === 0n) { d >>= 1n; s++; }
    
    const aMont = toMont64(a % n, n);
    if (aMont === 0n) return true;
    
    let x = montPow64(aMont, d, n, nInv, one);
    const negOne = n - one;
    
    if (x === one || x === negOne) return true;
    
    for (let r = 1n; r < s; r++) {
        x = montMul64(x, x, n, nInv);
        if (x === negOne) return true;
        if (x === one) return false;
    }
    return false;
}

/**
 * Deterministic primality test for 64-bit integers
 * Uses Montgomery-based Miller-Rabin with optimal witnesses
 */
export function isPrime64(n: number | bigint): boolean {
    const bn = BigInt(n);
    if (bn < 2n) return false;
    if (bn === 2n) return true;
    if ((bn & 1n) === 0n) return false;
    if (bn < 9n) return true;
    if (bn % 3n === 0n || bn % 5n === 0n || bn % 7n === 0n) return false;
    if (bn === 3215031751n) return false;
    
    // For small n, use 32-bit version
    if (bn <= MASK32) {
        return isPrime32(bn);
    }
    
    const nInv = montInv64(bn);
    const one = montOne64(bn);
    
    const witnesses = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n];
    for (const w of witnesses) {
        if (!montSprp64(bn, w, nInv, one)) return false;
    }
    return true;
}

/** Find the next prime >= n (64-bit) */
export function nextPrime64(n: number | bigint): bigint {
    let bn = BigInt(n);
    if (bn <= 2n) return 2n;
    if ((bn & 1n) === 0n) bn++;
    while (!isPrime64(bn)) {
        bn += 2n;
        if (bn > MASK64) return 0n;
    }
    return bn;
}

/** Find the previous prime <= n (64-bit) */
export function prevPrime64(n: number | bigint): bigint {
    let bn = BigInt(n);
    if (bn < 2n) return 0n;
    if (bn === 2n) return 2n;
    if ((bn & 1n) === 0n) bn--;
    while (!isPrime64(bn)) {
        if (bn < 3n) return 2n;
        bn -= 2n;
    }
    return bn;
}

// ============================================================================
// DEFAULT ALIASES (64-bit)
// ============================================================================

/** Primality test (defaults to 64-bit) */
export function isPrime(n: number | bigint): boolean {
    return isPrime64(n);
}

/** Find next prime (defaults to 64-bit) */
export function nextPrime(n: number | bigint): bigint {
    return nextPrime64(n);
}

/** Find previous prime (defaults to 64-bit) */
export function prevPrime(n: number | bigint): bigint {
    return prevPrime64(n);
}

// Test when run directly
if (typeof require !== 'undefined' && require.main === module) {
    console.log("Testing smcPrime...");
    
    // Test small primes
    const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
    for (const p of primes) {
        if (!isPrime(p)) {
            console.error(`FAIL: ${p} should be prime`);
            process.exit(1);
        }
    }
    
    // Test composites
    const composites = [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25];
    for (const c of composites) {
        if (isPrime(c)) {
            console.error(`FAIL: ${c} should not be prime`);
            process.exit(1);
        }
    }
    
    // Test special pseudoprime
    if (isPrime(3215031751)) {
        console.error("FAIL: 3215031751 should not be prime");
        process.exit(1);
    }
    
    // Test next/prev prime
    if (nextPrime(100) !== 101n) {
        console.error("FAIL: nextPrime(100) should be 101");
        process.exit(1);
    }
    if (prevPrime(100) !== 97n) {
        console.error("FAIL: prevPrime(100) should be 97");
        process.exit(1);
    }
    
    // Test large prime
    if (!isPrime64(1000000007n)) {
        console.error("FAIL: 1000000007 should be prime");
        process.exit(1);
    }
    
    console.log("All tests passed!");
}
