# smcprime

Ultra-fast primality testing with Montgomery arithmetic (32-bit and 64-bit).

[![npm version](https://badge.fury.io/js/smcprime.svg)](https://www.npmjs.com/package/smcprime)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- **Fast**: Montgomery arithmetic for 64-bit, native math for 32-bit
- **Deterministic**: No probabilistic results - always correct
- **TypeScript**: Full type definitions included
- **Zero dependencies**

## Installation

```bash
npm install smcprime
```

## Usage

```typescript
import { isPrime, isPrime32, isPrime64, nextPrime, prevPrime } from 'smcprime';

// Basic primality testing
console.log(isPrime(17));      // true
console.log(isPrime(18));      // false

// Find next/previous primes
console.log(nextPrime(100));   // 101n
console.log(prevPrime(100));   // 97n

// Explicit 32-bit or 64-bit
console.log(isPrime32(104729));        // true
console.log(isPrime64(1000000007n));   // true
```

## API

### Primality Testing
- `isPrime(n: number | bigint): boolean` - Test if n is prime (64-bit)
- `isPrime32(n: number | bigint): boolean` - Test if n is prime (32-bit)
- `isPrime64(n: number | bigint): boolean` - Test if n is prime (64-bit)

### Prime Navigation
- `nextPrime(n: number | bigint): bigint` - Find smallest prime >= n
- `prevPrime(n: number | bigint): bigint` - Find largest prime <= n
- `nextPrime32(n: number | bigint): bigint` - 32-bit version
- `prevPrime32(n: number | bigint): bigint` - 32-bit version
- `nextPrime64(n: number | bigint): bigint` - 64-bit version
- `prevPrime64(n: number | bigint): bigint` - 64-bit version

## Requirements

- ES2020+ (for BigInt support)
- Node.js 12+ or modern browsers

## License

MIT License - Copyright 2025 ScaleCode Solutions
