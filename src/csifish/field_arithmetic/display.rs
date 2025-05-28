use std::fmt::{Display, LowerHex, UpperHex};
use crate::csifish::field_arithmetic::arithmetic::MontgomeryArithmetic;

use crate::csifish::field_arithmetic::base_field::FieldElement;

impl Display for FieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::UpperHex::fmt(self, f)
    }
}

impl LowerHex for FieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for limb in self.get_standard().iter().rev() {
            std::fmt::LowerHex::fmt(limb, f)?;
        }
        Ok(())
    }
}

impl UpperHex for FieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for limb in self.get_standard().iter().rev() {
            write!(f, "{:016X}", limb)?;
        }
        Ok(())
    }
}


/// Decodes two hex characters into a single byte.
/// Returns the byte and an error flag.
#[inline(always)]
pub(crate) const fn decode_hex_byte(bytes: [u8; 2]) -> (u8, u16) {
    let hi = decode_nibble(bytes[0]);
    let lo = decode_nibble(bytes[1]);
    let byte = (hi << 4) | lo;
    let err = byte >> 8;
    let result = byte as u8;
    (result, err)
}

/// Decodes a single hex character into its numerical value.
/// Returns 0xFFFF if the character is invalid.
#[inline(always)]
const fn decode_nibble(src: u8) -> u16 {
    let byte = src as i16;
    let mut ret: i16 = -1;

    // 0-9  0x30-0x39
    // if (byte > 0x2f && byte < 0x3a) ret += byte - 0x30 + 1; // -47
    ret += (((0x2fi16 - byte) & (byte - 0x3a)) >> 8) & (byte - 47);
    // A-F  0x41-0x46
    // if (byte > 0x40 && byte < 0x47) ret += byte - 0x41 + 10 + 1; // -54
    ret += (((0x40i16 - byte) & (byte - 0x47)) >> 8) & (byte - 54);
    // a-f  0x61-0x66
    // if (byte > 0x60 && byte < 0x67) ret += byte - 0x61 + 10 + 1; // -86
    ret += (((0x60i16 - byte) & (byte - 0x67)) >> 8) & (byte - 86);

    ret as u16
}

