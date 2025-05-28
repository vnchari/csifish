use crypto_bigint::rand_core::{CryptoRng, RngCore};
use subtle::CtOption;

use crate::csifish::constants::DeserializationError;

pub trait ModularArithmetic: Sized + Clone + Copy {
    const LIMBS: usize;
    const MODULUS: &'static [u64];
    type Element;

    fn from_raw_limbs(l: [u64; Self::LIMBS]) -> Self::Element;
    fn from_u8(x: u8) -> Self::Element;
    fn from_u16(x: u16) -> Self::Element;
    fn from_be_hex(hex: &str) -> Self::Element;
    fn from_be_bytes(b: &[u8]) -> Result<Self, DeserializationError>;
    fn get_be_bytes(&self) -> [u8; Self::LIMBS * 8];

    fn random(rng: &mut (impl CryptoRng + RngCore)) -> Self::Element;
    fn random_under_half(rng: &mut (impl CryptoRng + RngCore)) -> Self::Element;

    fn neg(self) -> Self::Element;
    fn is_zero(&self) -> bool;

    fn add_reduce(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS];
    fn sub_reduce(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS];

    fn cmovz_array(c: u64, t: &mut [u64; Self::LIMBS], r: &[u64; Self::LIMBS]);
    fn cmovz_swap(c: u64, t: &mut [u64; Self::LIMBS], r: &mut [u64; Self::LIMBS]);

    fn conditional_move(&mut self, do_move: u64, b: &Self::Element);
    fn conditional_swap(&mut self, do_swap: u64, b: &mut Self::Element);

    fn vartime_is_less(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> bool;
}

pub trait MontgomeryArithmetic: ModularArithmetic {
    const INV: u64;

    fn from_limbs_into_montgomery(limbs: [u64; Self::LIMBS]) -> Self::Element;
    fn get_montgomery(&self) -> [u64; Self::LIMBS];
    fn get_standard(&self) -> [u64; Self::LIMBS];

    fn montgomery_mul(x: &[u64; Self::LIMBS], y: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS];
    fn montgomery_square(x: &[u64; Self::LIMBS]) -> [u64; Self::LIMBS];

    fn square(self) -> Self::Element;
    fn inv(self) -> CtOption<Self::Element>;

    fn vartime_exp(&self, pow: &u64) -> Self::Element;
}