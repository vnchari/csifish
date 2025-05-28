


#[inline(always)]
pub(crate) fn ct_pick64(which: u64, in1: u64, in2: u64) -> u64 {
    // Constant time select.
    // if pick == 0xFF..FF (out = in1)
    // if pick == 0 (out = in2)
    // else out is undefined.
    (in1 & which) | (in2 & !which)
}

pub(crate) trait ConstantTimeOps {
    fn ca(self, rhs: u64, carry: u64) -> (u64, u64);
    fn cs(self, rhs: u64, carry: u64) -> (u64, u64);

    fn wide_mul_sum(self, rhs: u64, add1: u64) -> (u64, u64);
    fn wide_mul_sum2(self, rhs: u64, add1: u64, add2: u64) -> (u64, u64);
}

impl ConstantTimeOps for u64 {
    #[inline(always)]
    fn ca(self, rhs: u64, carry: u64) -> (u64, u64) {
        let r = unsafe {
            (self as u128)
                .unchecked_add(rhs as u128)
                .unchecked_add(carry as u128)
        };
        (r as u64, (r >> 64) as u64)
    }
    #[inline(always)]
    fn cs(self, rhs: u64, carry: u64) -> (u64, u64) {
        let r = (self as u128)
            .wrapping_sub(rhs as u128)
            .wrapping_sub(carry as u128);
        (r as u64, ((r >> 64) as u64) >> 63)
    }

    // hi, lo
    #[inline(always)]
    fn wide_mul_sum(self, rhs: u64, add1: u64) -> (u64, u64) {
        //no overflow guaranteed
        let wide = unsafe {
            (self as u128)
                .unchecked_mul(rhs as u128)
                .unchecked_add(add1 as u128)
        };
        ((wide >> 64) as u64, wide as u64)
    }

    // hi, lo
    #[inline(always)]
    fn wide_mul_sum2(self, rhs: u64, add1: u64, add2: u64) -> (u64, u64) {
        //no overflow guaranteed
        let wide = unsafe {
            (self as u128)
                .unchecked_mul(rhs as u128)
                .unchecked_add(add1 as u128)
                .unchecked_add(add2 as u128)
        };
        ((wide >> 64) as u64, wide as u64)
    }
}


/// Compares two limb arrays in constant time.
/// Returns -1 if lhs < rhs, 0 if equal, 1 if lhs > rhs.
pub(crate) fn cmp_limbs_ct(lhs: &[u64; 8], rhs: &[u64; 8]) -> i8 {
    let mut w;
    let mut borrow = 0;
    let mut diff = 0;

    for i in 0..8 {
        (w, borrow) = rhs[i].cs(lhs[i], borrow);
        diff = diff | w;
    }
    let sgn = ((borrow as i32) | ((1 - borrow) as i32).wrapping_neg()) as i8;
    (ct_is_non_zero64(diff) as i8) * sgn
}

/// Checks if a 64-bit integer is non-zero in constant time.
/// Returns 0 if zero, 1 otherwise.
#[inline(always)]
pub(crate) fn ct_is_non_zero64(i: u64) -> i32 {
    // ct_is_non_zero64 returns 0 in case i == 0, otherwise it returns 1.
    // Constant-time.
    // In case i==0 then i-1 will set MSB. Only in such case (i OR ~(i-1))
    // will result in MSB being not set (logical implication: (i-1)=>i is
    // false iff (i-1)==0 and i==non-zero). In every other case MSB is
    // set and hence function returns 1.
    ((i | !i.wrapping_sub(1)) >> 63) as i32
}

#[link(name = "modinv")]
extern "C" {
    /// External function to compute the modular inverse in constant time.
    /// Calls into "s2n-bignum" library
    ///
    /// # Safety
    ///
    /// - `a` must be coprime with `modulus`.
    /// - `modulus` must be an odd number greater than 1.
    /// - The `buf` must have sufficient space (at least 3*k limbs).
    ///
    /// # Parameters
    ///
    /// - `k`: Number of limbs.
    /// - `out`: Pointer to the output limbs.
    /// - `a`: Pointer to the input limbs.
    /// - `modulus`: Pointer to the modulus limbs.
    /// - `buf`: Pointer to a temporary buffer.
    pub(crate) fn modinv(k: u64, out: *mut u64, a: *const u64, modulus: *const u64, buf: *mut u64);
}


/// Constant-time equality check for two limb arrays.
/// Returns true if equal, false otherwise.
pub(crate) fn ct_equal(v: &[u64; 8], in_fp: &[u64; 8]) -> bool {
    // equal checks if v is equal to in. Constant time.
    let mut r = 0u64;
    for i in 0..8 {
        r |= v[i] ^ in_fp[i];
    }
    ct_is_non_zero64(r) == 0
}