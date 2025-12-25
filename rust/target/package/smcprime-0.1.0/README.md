# smcprime

Ultra-fast primality testing with Montgomery arithmetic (32-bit and 64-bit).

[![Crates.io](https://img.shields.io/crates/v/smcprime.svg)](https://crates.io/crates/smcprime)
[![Documentation](https://docs.rs/smcprime/badge.svg)](https://docs.rs/smcprime)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- **Fast**: Montgomery arithmetic for 64-bit, native 64-bit math for 32-bit
- **Deterministic**: No probabilistic results - always correct
- **`no_std` compatible**: Works in embedded environments
- **Zero dependencies**

## Usage

```rust
use smcprime::{is_prime, is_prime32, is_prime64, next_prime, prev_prime};

// Basic primality testing
assert!(is_prime(17));
assert!(!is_prime(18));

// Find next/previous primes
assert_eq!(next_prime(100), 101);
assert_eq!(prev_prime(100), 97);

// Explicit 32-bit or 64-bit
assert!(is_prime32(104729));
assert!(is_prime64(1000000007));
```

## API

### Primality Testing
- `is_prime(n: u64) -> bool` - Test if n is prime (64-bit)
- `is_prime32(n: u32) -> bool` - Test if n is prime (32-bit)
- `is_prime64(n: u64) -> bool` - Test if n is prime (64-bit)

### Prime Navigation
- `next_prime(n: u64) -> u64` - Find smallest prime >= n
- `prev_prime(n: u64) -> u64` - Find largest prime <= n
- `next_prime32(n: u32) -> u32` - 32-bit version
- `prev_prime32(n: u32) -> u32` - 32-bit version
- `next_prime64(n: u64) -> u64` - 64-bit version
- `prev_prime64(n: u64) -> u64` - 64-bit version

## Performance

- **32-bit**: Uses witnesses {2, 7, 61} - only 3 Miller-Rabin rounds
- **64-bit**: Montgomery arithmetic avoids expensive modular division
- Benchmarks show ~3x faster than naive implementations

## Algorithm

- Miller-Rabin with deterministic witness sets
- Montgomery multiplication (no modular division)
- Hensel lifting for modular inverses
- Optimal witness selection from machine-prime research

## License

MIT License - Copyright 2025 ScaleCode Solutions
