use std::ops::{Add, AddAssign, Sub, SubAssign};

use ark_ff_macros::unroll_for_loops;
use crypto_bigint::rand_core::{CryptoRng, RngCore};

use crate::csifish::field_arithmetic::display::decode_hex_byte;
use crate::csifish::field_arithmetic::arithmetic::ModularArithmetic;
use crate::csifish::field_arithmetic::helpers::{ConstantTimeOps, ct_is_non_zero64, ct_pick64};
use crate::csifish::constants::{CLASSGROUP_ORDER, DeserializationError};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ClassGroupElement {
    pub limbs: [u64; 5],
}

impl ModularArithmetic for ClassGroupElement {
    const LIMBS: usize = 5;
    const MODULUS: &'static [u64] = &CLASSGROUP_ORDER;
    type Element = ClassGroupElement;

    fn from_raw_limbs(l: [u64; Self::LIMBS]) -> Self::Element {
        ClassGroupElement {
            limbs: l,
        }
    }

    fn from_u8(x: u8) -> Self::Element {
        let mut limbs = [0u64; Self::LIMBS];
        limbs[0] = x as u64;
        Self::from_raw_limbs(limbs)
    }

    fn from_u16(x: u16) -> Self::Element {
        let mut limbs = [0u64; Self::LIMBS];
        limbs[0] = x as u64;
        Self::from_raw_limbs(limbs)
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
        Self::from_raw_limbs(res)
    }

    fn from_be_bytes(b: &[u8]) -> Result<Self, DeserializationError> {
        let v = b
            .chunks_exact(8)
            .map(|x| Ok(<u64>::from_be_bytes(x.try_into()?)))
            .collect::<Result<Vec<u64>, DeserializationError>>()?;
        let value = v[..Self::LIMBS].try_into()?;
        Ok(ClassGroupElement { limbs: value })
    }

    fn get_be_bytes(&self) -> [u8; Self::LIMBS * 8] {
        // self.limbs
        //     .iter()
        //     .flat_map(|x| x.to_be_bytes())
        //     .collect::<Vec<u8>>()
        //     .try_into()
        //     .unwrap()
        unimplemented!()
    }

    fn random(rng: &mut (impl CryptoRng + RngCore)) -> Self::Element {
        let r = [0u64; Self::LIMBS];
        loop {
            let mut r = unsafe { std::mem::transmute::<[u64; Self::LIMBS], [u8; 40]>(r) };
            rng.fill_bytes(&mut r[..33]);
            let mut r = unsafe { std::mem::transmute::<[u8; 40], [u64; Self::LIMBS]>(r) };
            r[4] >>= 6;
            if Self::vartime_is_less(&r, &CLASSGROUP_ORDER) {
                return ClassGroupElement { limbs: r };
            };
        }
    }

    fn random_under_half(rng: &mut (impl CryptoRng + RngCore)) -> Self::Element {
        unimplemented!()
    }

    #[inline(always)]
    fn neg(self) -> Self::Element {
        Self::Element {
            limbs: Self::sub_reduce(&CLASSGROUP_ORDER, &self.limbs),
        }
    }

    fn is_zero(&self) -> bool {
        //constant time
        let mut r = 0u64;
        for &item in self.limbs.iter() {
            r |= item;
        }
        ct_is_non_zero64(r) == 0
    }

    fn add_reduce(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS] {
        let mut r = [0u64; Self::LIMBS];
        let mut t = [0u64; Self::LIMBS];
        let mut c = 0;
        (r[0], c) = x[0].ca(y[0], c);
        (r[1], c) = x[1].ca(y[1], c);
        (r[2], c) = x[2].ca(y[2], c);
        (r[3], c) = x[3].ca(y[3], c);
        (r[4], _) = x[4].ca(y[4], c);

        (t[0], c) = r[0].cs(CLASSGROUP_ORDER[0], 0);
        (t[1], c) = r[1].cs(CLASSGROUP_ORDER[1], c);
        (t[2], c) = r[2].cs(CLASSGROUP_ORDER[2], c);
        (t[3], c) = r[3].cs(CLASSGROUP_ORDER[3], c);
        (t[4], c) = r[4].cs(CLASSGROUP_ORDER[4], c);
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

        // we want x < y, so conditionally add p to the result
        let w = (1 - c).wrapping_sub(1);
        (r[0], c) = r[0].ca(ct_pick64(w, CLASSGROUP_ORDER[0], 0), 0);
        (r[1], c) = r[1].ca(ct_pick64(w, CLASSGROUP_ORDER[1], 0), c);
        (r[2], c) = r[2].ca(ct_pick64(w, CLASSGROUP_ORDER[2], 0), c);
        (r[3], c) = r[3].ca(ct_pick64(w, CLASSGROUP_ORDER[3], 0), c);
        (r[4], _) = r[4].ca(ct_pick64(w, CLASSGROUP_ORDER[4], 0), c);
        r
    }

    #[inline(always)]
    fn cmovz_array(c: u64, t: &mut [u64; Self::LIMBS], r: &[u64; Self::LIMBS]) {
        let m = ((c | (!c).wrapping_add(1)) >> 63) & 1;
        let mask = (1 ^ m).wrapping_sub(1);
        for i in 0..Self::LIMBS {
            t[i] = (t[i] & mask) | (r[i] & !mask);
        }
    }

    fn cmovz_swap(c: u64, t: &mut [u64; Self::LIMBS], r: &mut [u64; Self::LIMBS]) {
        let m = ((c | (!c).wrapping_add(1)) >> 63) & 1;
        let mask = (1 ^ m).wrapping_sub(1);
        for i in 0..Self::LIMBS {
            let tmp1 = t[i];
            t[i] = (t[i] & mask) | (r[i] & !mask);
            r[i] = (r[i] & mask) | (tmp1 & !mask);
        }
    }

    fn conditional_move(&mut self, do_move: u64, b: &Self::Element) {
        Self::cmovz_array(1 - do_move, &mut self.limbs, &b.limbs);
    }

    fn conditional_swap(&mut self, do_swap: u64, b: &mut Self::Element) {
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

impl Add for ClassGroupElement {
    type Output = ClassGroupElement;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        ClassGroupElement {
            limbs: Self::add_reduce(&self.limbs, &rhs.limbs),
        }
    }
}

impl Add for &ClassGroupElement {
    type Output = ClassGroupElement;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        ClassGroupElement {
            limbs: ClassGroupElement::add_reduce(&self.limbs, &rhs.limbs),
        }
    }
}

impl AddAssign for ClassGroupElement {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.limbs = ClassGroupElement::add_reduce(&self.limbs, &rhs.limbs);
    }
}

impl AddAssign<&Self> for ClassGroupElement {
    #[inline]
    fn add_assign(&mut self, rhs: &Self) {
        self.limbs = ClassGroupElement::add_reduce(&self.limbs, &rhs.limbs);
    }
}

impl Sub for ClassGroupElement {
    type Output = ClassGroupElement;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        ClassGroupElement {
            limbs: ClassGroupElement::sub_reduce(&self.limbs, &rhs.limbs),
        }
    }
}

impl Sub for &ClassGroupElement {
    type Output = ClassGroupElement;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        ClassGroupElement {
            limbs: ClassGroupElement::sub_reduce(&self.limbs, &rhs.limbs),
        }
    }
}

impl SubAssign for ClassGroupElement {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.limbs = ClassGroupElement::sub_reduce(&self.limbs, &rhs.limbs);
    }
}

impl SubAssign<&ClassGroupElement> for ClassGroupElement {
    #[inline]
    fn sub_assign(&mut self, rhs: &Self) {
        self.limbs = ClassGroupElement::sub_reduce(&self.limbs, &rhs.limbs);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_bigint::modular::{ConstMontyForm, ConstMontyParams};
    use crypto_bigint::{impl_modulus, Random, U320};
    use rand::thread_rng;
    impl_modulus!(
        ModulusClassGroup,
        U320,
        "000000000000000233002CB20D405A4F0C6DBD5A6A941DF1DF68A8029B289F124291AA03CD95356F"
    );
    pub type ModClassGroup = ConstMontyForm<ModulusClassGroup, { ModulusClassGroup::LIMBS }>;

    #[test]
    fn add_mod_class_group() {
        for _ in 0..500 {
            let p = ModClassGroup::random(&mut thread_rng());
            let q = ModClassGroup::random(&mut thread_rng());
            let p1 = ClassGroupElement::from_raw_limbs(p.retrieve().clone().to_words());
            let q1 = ClassGroupElement::from_raw_limbs(q.retrieve().clone().to_words());
            assert_eq!((p + q).retrieve().to_words(), (p1 + q1).limbs);
        }
    }

    #[test]
    fn sub_mod_class_group() {
        for _ in 0..500 {
            let p = ModClassGroup::random(&mut thread_rng());
            let q = ModClassGroup::random(&mut thread_rng());
            let p1 = ClassGroupElement::from_raw_limbs(p.retrieve().clone().to_words());
            let q1 = ClassGroupElement::from_raw_limbs(q.retrieve().clone().to_words());
            assert_eq!((p - q).retrieve().to_words(), (p1 - q1).limbs);
        }
    }
}
