# smcPrime

Ultra-fast primality testing library for 32-bit and 64-bit integers.

## Features

- **Header-only**: Single `smcprime.h` file, no dependencies
- **Deterministic**: Correct for ALL 32-bit and 64-bit integers
- **Ultra-fast**: Uses Montgomery arithmetic and prime-inverse trial division
- **Portable**: Works on x86, x64, ARM, ARM64, RISC-V

## Performance

Adopts optimizations from [machine-prime](https://github.com/JASory/machine-prime) (fastest Rust prime library):

- **Montgomery arithmetic** for Miller-Rabin (avoids modular division)
- **Prime-inverse trial division** (multiplication instead of %)
- **Hensel lifting** for modular inverses (faster than extended GCD)
- **Optimal witness selection** (minimal Miller-Rabin rounds)

Benchmarks (M4 Max):
- 32-bit: 100K tests in 0.003 seconds
- 64-bit: 100K tests in 0.003 seconds

## Usage

```c
#include "smcprime.h"

// 32-bit primality
bool is_prime = smc_is_prime32(2147483647);  // true (Mersenne M31)

// 64-bit primality  
bool is_prime = smc_is_prime64(18446744073709551557ULL);  // true (largest 64-bit prime)

// Default alias (64-bit)
bool is_prime = smc_is_prime(1000000007);  // true

// Utility functions
uint64_t next = smc_next_prime(100);   // 101
uint64_t prev = smc_prev_prime(100);   // 97
```

## API

### 32-bit Functions
- `smc_is_prime32(n)` - Test if n is prime
- `smc_next_prime32(n)` - Find next prime >= n
- `smc_prev_prime32(n)` - Find previous prime <= n

### 64-bit Functions
- `smc_is_prime64(n)` - Test if n is prime
- `smc_is_prime64_wc(n)` - Worst-case optimized (for likely primes)
- `smc_next_prime64(n)` - Find next prime >= n
- `smc_prev_prime64(n)` - Find previous prime <= n

### Default Aliases (64-bit)
- `smc_is_prime` → `smc_is_prime64`
- `smc_is_prime_wc` → `smc_is_prime64_wc`
- `smc_next_prime` → `smc_next_prime64`
- `smc_prev_prime` → `smc_prev_prime64`

## Algorithm Details

### 32-bit
Uses native 64-bit math with Miller-Rabin witnesses {2, 7, 61} (Jaeschke 1993).

### 64-bit
1. **Trial division** using prime inverses (primes 3-331)
2. **Montgomery Miller-Rabin** with witnesses {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}

Both are deterministic - no probabilistic results.

## License

MIT License - Copyright 2025 ScaleCode Solutions

MV❤️
