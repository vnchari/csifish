use sha3::digest::core_api::XofReaderCoreWrapper;
use sha3::digest::{ExtendableOutput, ExtendableOutputReset, Update};
use sha3::{TurboShake128, TurboShake128Core, TurboShake128ReaderCore};

pub type HashType = [u8; 16];
pub(crate) const HASH_SIZE: usize = std::mem::size_of::<HashType>();

pub(crate) struct Hasher {
    num_hashes: u32,
}

impl Hasher {
    pub fn new(num_hashes: u32) -> Hasher {
        Hasher { num_hashes }
    }
    pub fn hash(&self, input: &Vec<u8>) -> HashType {
        let mut result: HashType = HashType::default();
        let mut hasher = TurboShake128::from_core(<TurboShake128Core>::new(0x01));
        hasher.update(input); //
        hasher.finalize_xof_reset_into(&mut result);
        for _ in 0..self.num_hashes {
            hasher.update(&result);
            hasher.finalize_xof_reset_into(&mut result);
        }
        result
    }

    pub fn hash_extendable(
        &self,
        input: &Vec<u8>,
    ) -> XofReaderCoreWrapper<TurboShake128ReaderCore> {
        let mut result: HashType = HashType::default();
        let mut hasher = TurboShake128::from_core(<TurboShake128Core>::new(0x01));
        hasher.update(input); //
        hasher.finalize_xof_reset_into(&mut result);
        for _ in 0..(self.num_hashes - 1) {
            hasher.update(&result);
            hasher.finalize_xof_reset_into(&mut result);
        }
        hasher.update(&result);
        hasher.finalize_xof()
    }
}
