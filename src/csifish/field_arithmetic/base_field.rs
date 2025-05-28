use crate::csifish::field_arithmetic::helpers::{cmp_limbs_ct, ct_equal, modinv};
use std::arch::asm;
use std::cmp::Ordering;
use std::cmp::Ordering::{Equal, Greater, Less};
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};
use ark_ff_macros::unroll_for_loops;
use crypto_bigint::rand_core::{CryptoRng, RngCore};
use subtle::{Choice, CtOption};
use crate::csifish::field_arithmetic::arithmetic::{ModularArithmetic, MontgomeryArithmetic};
use crate::csifish::field_arithmetic::display::decode_hex_byte;
use crate::csifish::field_arithmetic::helpers::{ct_pick64, ConstantTimeOps, ct_is_non_zero64};
use crate::csifish::constants::{DeserializationError, ONE_MONTGOMERY, P_MINUS_ONE_OVER_TWO, PRIME, R2_MONTGOMERY, R3_MONTGOMERY};

#[derive(Copy, Clone, Debug)]
pub struct FieldElement {
    pub(crate) limbs: [u64; 8],
}

impl FieldElement {
    pub const ONE: Self = FieldElement { limbs: ONE_MONTGOMERY };
    pub const ZERO: Self = FieldElement { limbs: [0, 0, 0, 0, 0, 0, 0, 0] };


    /// Computes the field element raised to a bounded exponent in constant time.
    /// The exponent must be less than 2^10.
    pub fn constant_time_bounded_exp(&self, pow: &u64) -> FieldElement {
        let mut pow = pow.clone();
        let mut this = self.clone();
        let mut tmp = FieldElement::ONE;
        let mut res = FieldElement::ONE;
        for i in 0..11 {
            let done = (pow | !pow.wrapping_sub(1)) >> 63;
            res.conditional_move(1 - done, &tmp);
            if i == 10 {
                break;
            }
            tmp.conditional_move(pow % 2, &(tmp * this));
            this = this.square();
            pow >>= 1;
        }
        res
    }
}

impl ModularArithmetic for FieldElement {
    const LIMBS: usize = 8;
    const MODULUS: &'static [u64] = &PRIME;
    type Element = FieldElement;

    fn from_raw_limbs(l: [u64; Self::LIMBS]) -> Self::Element {
        FieldElement {
                limbs: l,
        }
    }

    fn from_u8(x: u8) -> Self::Element {
        let mut limbs = [0u64; Self::LIMBS];
        limbs[0] = x as u64;
        Self::from_limbs_into_montgomery(limbs)
    }

    fn from_u16(x: u16) -> Self::Element {
        let mut limbs = [0u64; Self::LIMBS];
        limbs[0] = x as u64;
        Self::from_limbs_into_montgomery(limbs)
    }

    fn from_be_hex(hex: &str) -> Self::Element {
        let bytes = hex.as_bytes();
        assert_eq!(bytes.len(), 128, "hex string is not the expected size");

        let mut res = [0u64; Self::LIMBS];
        let mut buf = [0u8; 8];
        let mut err = 0;

        for i in 0..Self::LIMBS {
            let mut j = 0;
            while j < Self::LIMBS {
                let offset = (i * Self::LIMBS + j) * 2;
                let (result, byte_err) = decode_hex_byte([bytes[offset], bytes[offset + 1]]);
                err |= byte_err;
                buf[j] = result;
                j += 1;
            }
            res[Self::LIMBS - i - 1] = <u64>::from_be_bytes(buf);
        }
        assert_eq!(err, 0, "invalid hex byte");
        Self::from_limbs_into_montgomery(res)
    }

    fn from_be_bytes(b: &[u8]) -> Result<Self, DeserializationError> {
        let v = b
            .chunks_exact(8)
            .map(|x| Ok(<u64>::from_be_bytes(x.try_into()?)))
            .collect::<Result<Vec<u64>, DeserializationError>>()?;
        let limbs = v[..8].try_into()?;
        Ok(FieldElement { limbs })
    }

    fn get_be_bytes(&self) -> [u8; 64] {
        // unimplemented!();
        self.get_standard()
            .iter()
            .flat_map(|x| x.to_be_bytes())
            .collect::<Vec<u8>>()
            .try_into()
            .unwrap()
    }

    fn random(rng: &mut (impl CryptoRng + RngCore)) -> Self::Element {
        let mut r = [0u64; Self::LIMBS];
        loop {
            for i in 0..8 {
                r[i] = rng.next_u64();
            }
            if Self::vartime_is_less(&r, &PRIME) {
                return FieldElement::from_limbs_into_montgomery(r);
            };
        }
    }

    fn random_under_half(rng: &mut (impl CryptoRng + RngCore)) -> Self::Element {
        let mut r = [0u64; Self::LIMBS];
        loop {
            for i in 0..Self::LIMBS {
                r[i] = rng.next_u64();
            }
            // Ensure the top bit of the top limb is not set to keep the element below half.
            r[7] >>= 1;
            if Self::vartime_is_less(&r, &P_MINUS_ONE_OVER_TWO) {
                return FieldElement::from_limbs_into_montgomery(r);
            };
        }
    }
    #[inline(always)]
    fn neg(self) -> Self::Element {
        Self::from_raw_limbs(Self::sub_reduce(&PRIME, &self.limbs))
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        //constant time
        let mut r = 0u64;
        for &item in self.limbs.iter() {
            r |= item;
        }
        ct_is_non_zero64(r) == 0
    }

    #[inline(always)]
    fn add_reduce(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS] {
        let mut r = [0u64; Self::LIMBS];
        let mut t = [0u64; Self::LIMBS];
        let mut c: u64 = 0;
        (r[0], c) = x[0].ca(y[0], c);
        (r[1], c) = x[1].ca(y[1], c);
        (r[2], c) = x[2].ca(y[2], c);
        (r[3], c) = x[3].ca(y[3], c);
        (r[4], c) = x[4].ca(y[4], c);
        (r[5], c) = x[5].ca(y[5], c);
        (r[6], c) = x[6].ca(y[6], c);
        (r[7], _) = x[7].ca(y[7], c);

        (t[0], c) = r[0].cs(PRIME[0], 0u64);
        (t[1], c) = r[1].cs(PRIME[1], c);
        (t[2], c) = r[2].cs(PRIME[2], c);
        (t[3], c) = r[3].cs(PRIME[3], c);
        (t[4], c) = r[4].cs(PRIME[4], c);
        (t[5], c) = r[5].cs(PRIME[5], c);
        (t[6], c) = r[6].cs(PRIME[6], c);
        (t[7], c) = r[7].cs(PRIME[7], c);
        Self::cmovz_array(c, &mut r, &t);
        r
    }

    #[inline(always)]
    fn sub_reduce(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS] {
        let mut r = [0u64; Self::LIMBS];
        let mut c = 0;
        // // Perform the subtraction x - y
        (r[0], c) = x[0].cs(y[0], c);
        (r[1], c) = x[1].cs(y[1], c);
        (r[2], c) = x[2].cs(y[2], c);
        (r[3], c) = x[3].cs(y[3], c);
        (r[4], c) = x[4].cs(y[4], c);
        (r[5], c) = x[5].cs(y[5], c);
        (r[6], c) = x[6].cs(y[6], c);
        (r[7], c) = x[7].cs(y[7], c);

        // we want x < y, so conditionally add p to the result
        let w = (1 - c).wrapping_sub(1);
        (r[0], c) = r[0].ca(ct_pick64(w, PRIME[0], 0), 0);
        (r[1], c) = r[1].ca(ct_pick64(w, PRIME[1], 0), c);
        (r[2], c) = r[2].ca(ct_pick64(w, PRIME[2], 0), c);
        (r[3], c) = r[3].ca(ct_pick64(w, PRIME[3], 0), c);
        (r[4], c) = r[4].ca(ct_pick64(w, PRIME[4], 0), c);
        (r[5], c) = r[5].ca(ct_pick64(w, PRIME[5], 0), c);
        (r[6], c) = r[6].ca(ct_pick64(w, PRIME[6], 0), c);
        (r[7], _) = r[7].ca(ct_pick64(w, PRIME[7], 0), c);
        r
    }

    #[cfg(target_arch = "aarch64")]
    #[inline(always)]
    /// Conditional move for ARM64 architecture using inline assembly.
    /// Moves `r` to `t` if `c` is zero.
    fn cmovz_array(c: u64, t: &mut [u64; Self::LIMBS], r: &[u64; Self::LIMBS]) {
        let c = c as u8;
        unsafe {
            asm! {
            "cmp {0:w}, 0",
            "csel {1:x}, {5:x}, {1:x}, EQ",
            "csel {2:x}, {6:x}, {2:x}, EQ",
            "csel {3:x}, {7:x}, {3:x}, EQ",
            "csel {4:x}, {8:x}, {4:x}, EQ",
            in(reg) c,
            inlateout(reg) t[0],
            inlateout(reg) t[1],
            inlateout(reg) t[2],
            inlateout(reg) t[3],
            in(reg) r[0],
            in(reg) r[1],
            in(reg) r[2],
            in(reg) r[3],
            options(pure, nomem, nostack),
            };
            asm! {
            "cmp {8:w}, 0",
            "csel {0:x}, {4:x}, {0:x}, EQ",
            "csel {1:x}, {5:x}, {1:x}, EQ",
            "csel {2:x}, {6:x}, {2:x}, EQ",
            "csel {3:x}, {7:x}, {3:x}, EQ",
            inlateout(reg) t[4],
            inlateout(reg) t[5],
            inlateout(reg) t[6],
            inlateout(reg) t[7],
            in(reg) r[4],
            in(reg) r[5],
            in(reg) r[6],
            in(reg) r[7],
            in(reg) c,
            options(pure, nomem, nostack),
            };
        }
    }

    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86", target_arch = "x86_64")))]
    #[inline(always)]
    /// Generic conditional move for architectures other than ARM64 and x86/x86_64.
    /// Moves `r` to `t` if `c` is zero.
    fn cmovz_array(c: u64, t: &mut [u64; Self::LIMBS], r: &[u64; Self::LIMBS]) {
        let m = ((c | (!c).wrapping_add(1)) >> 63) & 1;
        let mask = (1 ^ m).wrapping_sub(1);
        for i in 0..Self::LIMBS {
            t[i] = (t[i] & mask) | (r[i] & !mask);
        }
    }

    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86", target_arch = "x86_64")))]
    #[inline(always)]
    /// Generic conditional move for architectures other than ARM64 and x86/x86_64.
    /// Moves `r` to `t` if `c` is zero.
    fn cmovz_array(c: u64, t: &mut [u64; Self::LIMBS], r: &[u64; Self::LIMBS]) {
        let m = ((c | (!c).wrapping_add(1)) >> 63) & 1;
        let mask = (1 ^ m).wrapping_sub(1);
        for i in 0..Self::LIMBS {
            t[i] = (t[i] & mask) | (r[i] & !mask);
        }
    }


    /// Conditional swap of two field elements based on `do_swap`.
    /// If `do_swap` is 1, swaps `t` and `r`; otherwise, does nothing.
    #[inline(always)]
    fn cmovz_swap(c: u64, t: &mut [u64; Self::LIMBS], r: &mut [u64; Self::LIMBS]) {
        let m = ((c | (!c).wrapping_add(1)) >> 63) & 1;
        let mask = (1 ^ m).wrapping_sub(1);
        for i in 0..Self::LIMBS {
            let tmp1 = t[i];
            t[i] = (t[i] & mask) | (r[i] & !mask);
            r[i] = (r[i] & mask) | (tmp1 & !mask);
        }
    }


    /// Conditionally moves `b` into `self` based on `do_move`.
    /// If `do_move` is 1, `self` becomes `b`; otherwise, it remains unchanged.
    fn conditional_move(&mut self, do_move: u64, b: &Self) {
        //move on zero, so flip the condition
        Self::cmovz_array(1 - do_move, &mut self.limbs, &b.limbs);
    }

    /// Conditionally swaps `self` with `b` based on `do_swap`.
    /// If `do_swap` is 1, `self` and `b` are swapped; otherwise, they remain unchanged.
    fn conditional_swap(&mut self, do_swap: u64, b: &mut Self) {
        //move on zero, so flip the condition
        Self::cmovz_swap(1 - do_swap, &mut self.limbs, &mut b.limbs);
    }

    #[unroll_for_loops(10)]
    fn vartime_is_less(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> bool {
        // Returns result of x<y operation.
        for i in (0..Self::LIMBS).rev() {
            let (v, c) = y[i].overflowing_sub(x[i]);
            if c {
                return false;
            }
            if v != 0 {
                return true;
            }
        }
        // x == y
        false
    }
}


/// Internal representation is assumed to always be in Montgomery Form
impl MontgomeryArithmetic for FieldElement {
    const INV: u64 = 0x66c1301f632e294d;

    fn from_limbs_into_montgomery(limbs: [u64; Self::LIMBS]) -> Self::Element {
        Self::Element {
            limbs: Self::montgomery_mul(&limbs, &R2_MONTGOMERY),
        }
    }

    fn get_montgomery(&self) -> [u64; Self::LIMBS] {
        self.limbs
    }

    fn get_standard(&self) -> [u64; Self::LIMBS] {
        let mut one = [0u64; 8];
        one[0] = 1;
        Self::montgomery_mul(&self.limbs, &one)
    }

    #[inline(always)]
    #[unroll_for_loops(10)]
    fn montgomery_mul(lhs: &[u64; Self::LIMBS], rhs: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS] {
        // TODO: Add comments + citation
        let mut t = [0u64; 8];
        let (mut high, mut carry);
        for i in 0..8 {
            (high, t[0]) = lhs[0].wide_mul_sum(rhs[i], t[0]);
            let m = t[0].wrapping_mul(0x66c1301f632e294d); // -prime^-1 mod 2^64
            (carry, _) = m.wide_mul_sum(PRIME[0], t[0]);
            for j in 1..8 {
                (high, t[j]) = lhs[j].wide_mul_sum2(rhs[i], high, t[j]);
                (carry, t[j - 1]) = m.wide_mul_sum2(PRIME[j], carry, t[j]);
            }
            t[7] = unsafe { carry.unchecked_add(high) };
        }
        let mut c = 0;
        let mut r: [u64; 8] = [0; 8];
        (r[0], c) = t[0].cs(PRIME[0], c);
        (r[1], c) = t[1].cs(PRIME[1], c);
        (r[2], c) = t[2].cs(PRIME[2], c);
        (r[3], c) = t[3].cs(PRIME[3], c);
        (r[4], c) = t[4].cs(PRIME[4], c);
        (r[5], c) = t[5].cs(PRIME[5], c);
        (r[6], c) = t[6].cs(PRIME[6], c);
        (r[7], c) = t[7].cs(PRIME[7], c);
        Self::cmovz_array(c, &mut t, &r);
        t
    }

    #[inline(always)]
    #[unroll_for_loops(10)]
    fn montgomery_square(a: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS] {
        // TODO: Add comments + citation
        let mut r = [0u64; 16];
        let mut carry = 0;
        for i in 0..7 {
            for j in (i + 1)..8 {
                (carry, r[i + j]) = a[i].wide_mul_sum2(a[j], r[i + j], carry);
            }
            r[i + 8] = carry;
            carry = 0;
        }
        r[15] = r[14] >> 63;
        for i in 2..15 {
            r[16 - i] = (r[16 - i] << 1) | (r[15 - i] >> 63);
        }
        r[1] <<= 1;

        for i in 0..8 {
            (carry, r[2 * i]) = a[i].wide_mul_sum2(a[i], r[2 * i], carry);
            let tmp = unsafe { (r[2 * i + 1] as u128).unchecked_add(carry as u128) };
            carry = (tmp >> 64) as u64;
            r[2 * i + 1] = tmp as u64;
        }

        let mut carry2 = 0;
        for i in 0..8 {
            let k = r[i].wrapping_mul(0x66c1301f632e294d);
            (carry, _) = PRIME[0].wide_mul_sum(k, r[i]);
            for j in 1..8 {
                (carry, r[j + i]) = k.wide_mul_sum2(PRIME[j], r[j + i], carry);
            }
            let tmp = (r[i + 8] as u128) + (carry as u128) + (carry2 as u128);
            r[i + 8] = tmp as u64;
            carry2 = (tmp >> 64) as u64;
        }
        let mut c = 0;
        (r[0], c) = r[8].cs(PRIME[0], c);
        (r[1], c) = r[9].cs(PRIME[1], c);
        (r[2], c) = r[10].cs(PRIME[2], c);
        (r[3], c) = r[11].cs(PRIME[3], c);
        (r[4], c) = r[12].cs(PRIME[4], c);
        (r[5], c) = r[13].cs(PRIME[5], c);
        (r[6], c) = r[14].cs(PRIME[6], c);
        (r[7], c) = r[15].cs(PRIME[7], c);
        let mut t = unsafe { r[8..].try_into().unwrap_unchecked() };
        let r = unsafe { r[..8].try_into().unwrap_unchecked() };
        Self::cmovz_array(c, &mut t, &r);
        t
    }

    #[inline]
    fn square(self) -> Self::Element {
        Self::from_raw_limbs(Self::montgomery_square(&self.limbs))
    }

    /// Computes the multiplicative inverse using the external `modinv` function.
    /// The result is converted back to Montgomery representation.
    #[inline]
    fn inv(self) -> CtOption<FieldElement> {
        let mut result = [0; Self::LIMBS];
        let mut buf = [0; Self::LIMBS * 3]; // Buffer should be at least 3*k limbs.
        unsafe {
            modinv(
                8u64,
                result.as_mut_ptr(),
                self.limbs.as_ptr(),
                PRIME.as_ptr(),
                buf.as_mut_ptr(),
            );
        }

        CtOption::new(
            FieldElement::from_raw_limbs(Self::montgomery_mul(&result, &R3_MONTGOMERY)),
            Choice::from(!self.is_zero() as u8),
        )
    }

    /// Perform exponentiation by squaring, obviously variable time
    fn vartime_exp(&self, pow: &u64) -> FieldElement {
        // replace x with x^k
        let mut pow = pow.clone();
        let mut this = self.clone();
        let mut result = FieldElement::ONE;
        while pow != 0 {
            if pow % 2 == 1 {
                result *= this;
            }
            this = this.square();
            pow >>= 1;
        }
        result
    }
}

impl Add for FieldElement {
    type Output = FieldElement;
    #[inline]
    fn add(self, rhs: Self) -> FieldElement {
        FieldElement {
            limbs: Self::add_reduce(&self.limbs, &rhs.limbs),
        }
    }
}
impl Add for &FieldElement {
    type Output = FieldElement;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        FieldElement {
            limbs: FieldElement::add_reduce(&self.limbs, &rhs.limbs)
        }
    }
}
impl AddAssign for FieldElement {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.limbs = Self::add_reduce(&self.limbs, &rhs.limbs);
    }
}
impl AddAssign<&Self> for FieldElement {
    #[inline]
    fn add_assign(&mut self, rhs: &Self) {
        self.limbs = Self::add_reduce(&self.limbs, &rhs.limbs);
    }
}
impl Sub for FieldElement {
    type Output = FieldElement;
    #[inline]
    fn sub(self, rhs: Self) -> FieldElement {
        let limbs = Self::sub_reduce(&self.limbs, &rhs.limbs);
        FieldElement { limbs }
    }
}
impl Sub for &FieldElement {
    type Output = FieldElement;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        FieldElement {
            limbs: FieldElement::sub_reduce(&self.limbs, &rhs.limbs),
        }
    }
}
impl SubAssign for FieldElement {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.limbs = Self::sub_reduce(&self.limbs, &rhs.limbs);
    }
}
impl SubAssign<&FieldElement> for FieldElement {
    #[inline]
    fn sub_assign(&mut self, rhs: &Self) {
        self.limbs = FieldElement::sub_reduce(&self.limbs, &rhs.limbs);
    }
}
impl Mul for FieldElement {
    type Output = FieldElement;

    //assuming both self, rhs are in montgomery form
    fn mul(self, rhs: Self) -> FieldElement {
        FieldElement {
            limbs: Self::montgomery_mul(&self.limbs, &rhs.limbs),
        }
    }
}
impl Mul for &FieldElement {
    type Output = FieldElement;

    //assuming both self, rhs are in montgomery form
    fn mul(self, rhs: Self) -> FieldElement {
        FieldElement {
            limbs: FieldElement::montgomery_mul(&self.limbs, &rhs.limbs),
        }
    }
}
impl Mul<&Self> for FieldElement {
    type Output = FieldElement;

    //assuming both self, rhs are in montgomery form
    fn mul(self, rhs: &Self) -> FieldElement {
        FieldElement {
            limbs: Self::montgomery_mul(&self.limbs, &rhs.limbs),
        }
    }
}
impl MulAssign for FieldElement {
    //assuming both self, rhs are in montgomery form
    fn mul_assign(&mut self, rhs: Self) {
        self.limbs = Self::montgomery_mul(&self.limbs, &rhs.limbs)
    }
}
impl MulAssign<&FieldElement> for FieldElement {
    //assuming both self, rhs are in montgomery form
    fn mul_assign(&mut self, rhs: &Self) {
        self.limbs = FieldElement::montgomery_mul(&self.limbs, &rhs.limbs)
    }
}

impl PartialEq<Self> for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        //constant time
        ct_equal(&self.limbs, &other.limbs)
    }
}

impl PartialOrd for FieldElement {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        //jump table, should be constant time
        match cmp_limbs_ct(&self.limbs, &other.limbs) {
            -1 => Some(Less),
            0 => Some(Equal),
            1 => Some(Greater),
            _ => None,
        }
    }

    fn lt(&self, other: &Self) -> bool {
        cmp_limbs_ct(&self.limbs, &other.limbs) == -1i8
    }

    fn le(&self, other: &Self) -> bool {
        cmp_limbs_ct(&self.limbs, &other.limbs) <= 0i8
    }

    fn gt(&self, other: &Self) -> bool {
        cmp_limbs_ct(&self.limbs, &other.limbs) == 1i8
    }

    fn ge(&self, other: &Self) -> bool {
        cmp_limbs_ct(&self.limbs, &other.limbs) >= 0i8
    }
}


#[cfg(test)]
mod tests {
    use crypto_bigint::modular::{ConstMontyForm, ConstMontyParams};
    use crypto_bigint::{impl_modulus, Random, U512};
    use rand::thread_rng;

    use super::*;

    impl_modulus!(ModulusP, U512, "65b48e8f740f89bffc8ab0d15e3e4c4ab42d083aedc88c425afbfcc69322c9cda7aac6c567f35507516730cc1f0b4f25c2721bf457aca8351b81b90533c6c87b");
    pub type ModP = ConstMontyForm<ModulusP, { ModulusP::LIMBS }>;

    #[test]
    fn add() {
        for _ in 0..500 {
            let p = ModP::random(&mut thread_rng());
            let q = ModP::random(&mut thread_rng());
            let p1 = FieldElement::from_raw_limbs(p.as_montgomery().clone().to_words());
            let q1 = FieldElement::from_raw_limbs(q.as_montgomery().clone().to_words());
            assert_eq!(
                (p + q).as_montgomery().to_words(),
                (p1 + q1).get_montgomery()
            );
        }
    }

    #[test]
    fn add_assign() {
        for _ in 0..500 {
            let mut p = FieldElement::random(&mut thread_rng());
            let q = FieldElement::random(&mut thread_rng());
            let j =
                ModP::new(&U512::from(p.limbs.clone())) + ModP::new(&U512::from(q.limbs.clone()));
            p += q;
            assert_eq!(j.retrieve().to_words(), (p).limbs);
        }
    }

    #[test]
    fn sub() {
        for _ in 0..500 {
            let p = ModP::random(&mut thread_rng());
            let q = ModP::random(&mut thread_rng());
            let p1 = FieldElement::from_raw_limbs(p.as_montgomery().clone().to_words());
            let q1 = FieldElement::from_raw_limbs(q.as_montgomery().clone().to_words());
            assert_eq!(
                (p - q).as_montgomery().to_words(),
                (p1 - q1).get_montgomery()
            );
        }
    }

    #[test]
    fn sub_assign() {
        for _ in 0..500 {
            let mut p = FieldElement::random(&mut thread_rng());
            let q = FieldElement::random(&mut thread_rng());
            let j =
                ModP::new(&U512::from(p.limbs.clone())) - ModP::new(&U512::from(q.limbs.clone()));
            p -= q;
            assert_eq!(j.retrieve().to_words(), p.limbs);
        }
    }

    #[test]
    fn monty_mul() {
        assert_eq!(ModulusP::R2.to_words(), R2_MONTGOMERY);
        assert_eq!(ModulusP::MOD_NEG_INV.0, 0x66c1301f632e294d);
        for _ in 0..500 {
            let p = ModP::random(&mut thread_rng());
            let q = ModP::random(&mut thread_rng());
            let p1 = FieldElement {
                limbs: p.as_montgomery().clone().to_words(),
            };
            let q1 = FieldElement {
                limbs: q.as_montgomery().clone().to_words(),
            };
            assert_eq!((p1 * q1).limbs, (p * q).as_montgomery().to_words());
        }
    }

    #[test]
    fn equality() {
        for _ in 0..500 {
            let p = FieldElement::random(&mut thread_rng());
            let q = FieldElement::random(&mut thread_rng());
            assert!(p == p);
            assert!(p != q);
        }
    }

    #[test]
    fn cmp() {
        let p = U512::random(&mut thread_rng());
        assert_eq!(
            0,
            cmp_limbs_ct(&p.clone().to_words(), &p.clone().to_words())
        );
        for _ in 0..500 {
            let p = U512::random(&mut thread_rng());
            let q = U512::random(&mut thread_rng());
            let j = cmp_limbs_ct(&p.clone().to_words(), &q.clone().to_words());
            assert_eq!(p.cmp(&q) as i8, j);
        }
    }

    #[test]
    fn square() {
        for _ in 0..50 {
            let mut result = ModP::random(&mut thread_rng());
            let mut result_monty = FieldElement::from_limbs_into_montgomery(result.retrieve().to_words());

            for _ in 0..50 {
                // let mut result = ModP::new(&U512::from_u8(1));
                result = result.square();
                result_monty = result_monty.square();
            }
            assert_eq!(result_monty.get_standard(), result.retrieve().to_words());
        }
    }

    #[test]
    fn neg() {
        for _ in 0..500 {
            let result = ModP::random(&mut thread_rng());
            let result_monty = FieldElement::from_limbs_into_montgomery(result.retrieve().to_words());
            assert_eq!(
                result.neg().retrieve().to_words(),
                result_monty.neg().get_standard()
            );
        }
    }

    #[test]
    fn inv() {
        for _ in 0..500 {
            let result = ModP::random(&mut thread_rng());
            let result_monty = FieldElement::from_limbs_into_montgomery(result.retrieve().to_words());
            let inv_field = result_monty.inv().unwrap();
            let inv = result.inv().unwrap();
            assert_eq!(inv.retrieve().to_words(), inv_field.get_standard());
        }
    }

    #[test]
    fn var_exp() {
        for _ in 0..500 {
            let random = thread_rng().next_u32() >> 22;
            let result = ModP::random(&mut thread_rng());
            let result_monty = FieldElement::from_limbs_into_montgomery(result.retrieve().to_words());
            let pow = result_monty.vartime_exp(&(random as u64));
            let pow2 = result.pow_bounded_exp(&U512::from_u32(random), 10);
            assert_eq!(pow2.retrieve().to_words(), pow.get_standard());
        }
    }

    #[test]
    fn swap() {
        let mut a = [1u64; 8];
        let mut b = [2u64; 8];
        FieldElement::cmovz_swap(0, &mut a, &mut b);
        assert_eq!(a, [2u64; 8]);
        assert_eq!(b, [1u64; 8]);
        let mut a = [1u64; 8];
        let mut b = [2u64; 8];
        FieldElement::cmovz_swap(1, &mut a, &mut b);
        assert_eq!(a, [1u64; 8]);
        assert_eq!(b, [2u64; 8]);
    }

    #[test]
    fn cmov() {
        let p1 = FieldElement::random(&mut thread_rng());
        let p2 = FieldElement::random(&mut thread_rng());
        let mut p_test = p1.clone();
        p_test.conditional_move(0, &p2);
        assert_eq!(p1, p_test);
        let p1 = FieldElement::random(&mut thread_rng());
        let p2 = FieldElement::random(&mut thread_rng());
        let mut p_test = p1.clone();
        p_test.conditional_move(1, &p2);
        assert_eq!(p2, p_test);
    }

    #[test]
    fn ct_exp() {
        for _ in 0..500 {
            let random = thread_rng().next_u32() >> 22;
            let result = ModP::random(&mut thread_rng());
            let result_monty = FieldElement::from_limbs_into_montgomery(result.retrieve().to_words());
            let pow = result_monty.constant_time_bounded_exp(&(random as u64));
            let pow2 = result.pow_bounded_exp(&U512::from_u32(random), 10);
            assert_eq!(pow2.retrieve().to_words(), pow.get_standard());
        }
    }

}
