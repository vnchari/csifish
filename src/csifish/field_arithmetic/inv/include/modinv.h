// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT-0

// ----------------------------------------------------------------------------
// C prototypes for s2n-bignum functions, so you can use them in C programs via
//
//  #include "s2n-bignum.h"
//
// The functions are listed in alphabetical order with a brief description
// in comments for each one. For more detailed documentation see the comment
// banner at the top of the corresponding assembly (.S) file, and
// for the last word in what properties it satisfies see the spec in the
// formal proof (the .ml file in the architecture-specific directory).
//
// For some functions there are additional variants with names ending in
// "_alt". These have the same core mathematical functionality as their
// non-"alt" versions, but can be better suited to some microarchitectures:
//
//      - On x86, the "_alt" forms avoid BMI and ADX instruction set
//        extensions, so will run on any x86_64 machine, even older ones
//
//      - On ARM, the "_alt" forms target machines with higher multiplier
//        throughput, generally offering higher performance there.
//        The "_neon" forms target machines with NEON instructions.
// ----------------------------------------------------------------------------

#ifdef __APPLE__
#define S2N_BN_SYMBOL(NAME) _##NAME
#else
#define S2N_BN_SYMBOL(name) name
#endif

#if defined(_MSC_VER) || !defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L || defined(__STDC_NO_VLA__)
#define S2N_BIGNUM_STATIC
#else
#define S2N_BIGNUM_STATIC static
#endif

// Invert modulo m, z = (1/a) mod b, assuming b is an odd number > 1, a coprime to b
// Inputs a[k], b[k]; output z[k]; temporary buffer t[>=3*k]
//extern void modinv (uint64_t k, uint64_t *z, uint64_t *a, uint64_t *b, uint64_t *t);
