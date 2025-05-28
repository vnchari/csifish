use std::io::Read;
use rand::thread_rng;
use rayon::prelude::*;
use signature::{Error, Signer, Verifier};

use crate::csifish::field_arithmetic::classgroup::ClassGroupElement;

use crate::csifish::constant_time::{OneTimeAction, OneTimeCurve};
use crate::csifish::field_arithmetic::arithmetic::ModularArithmetic;
use crate::csifish::hash::{Hasher, HashType};
use crate::csifish::merkle::{ClassGroupMerkleProof, ClassGroupMerkleTree};
use crate::csifish::montgomery::MontgomeryCurve;
use crate::csifish::variable_time::VariableTimeAction;
use crate::csifish::constants::BASE_CURVE;

pub struct SigningKey<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> {
    proof_tree: ClassGroupMerkleTree<CURVES, ROUNDS, HASHES>,
    public_curves: Vec<MontgomeryCurve>,
    secret_actions: Vec<ClassGroupElement>,
}

impl<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> SigningKey<CURVES, ROUNDS, HASHES> {
    pub fn generate() -> SigningKey<CURVES, ROUNDS, HASHES> {
        let (cge, curves): (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) = Self::one_time_curves(CURVES as usize);
        let _tree = ClassGroupMerkleTree::from_leaves(&curves);
        SigningKey {
            proof_tree: _tree,
            public_curves: curves,
            secret_actions: cge,
        }
    }

    pub fn verifying_key(&self) -> VerifyingKey {
        VerifyingKey {
            root: self.proof_tree.root(),
            merkle_key: self.proof_tree.merkle_key(),
        }
    }

    fn one_time_curves(num_curves: usize) -> (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) {
        (0..num_curves).into_par_iter().map(|_| {
            let r = ClassGroupElement::random(&mut thread_rng());
            let curve = r.reduce().one_time_blinded_action(&BASE_CURVE);
            (r, curve.normalize())
        }).unzip()
    }

    fn variable_time_curves(num_curves: usize) -> (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) {
        (0..num_curves).into_par_iter().map(|_| {
            let r = ClassGroupElement::random(&mut thread_rng());
            let curve = r.reduce().variable_time_action(&BASE_CURVE);
            (r, curve)
        }).unzip()
    }
}

pub struct VerifyingKey {
    root: HashType,
    merkle_key: HashType,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Signature<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> {
    num_curves: u32,
    challenges: Vec<u8>,
    ephemeral_cge: Vec<ClassGroupElement>,
    opened_curves: Vec<MontgomeryCurve>,
    proof: ClassGroupMerkleProof,
}

impl<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> Signer<Signature<CURVES, ROUNDS, HASHES>> for SigningKey<CURVES, ROUNDS, HASHES> {
    fn try_sign(&self, message: &[u8]) -> Result<Signature<CURVES, ROUNDS, HASHES>, Error> {
        let (b, ephemeral_curves): (Vec<ClassGroupElement>, Vec<MontgomeryCurve>) = Self::variable_time_curves(ROUNDS as usize);
        let mut v: Vec<u8> = ephemeral_curves.into_iter().map(|x| x.to_be_bytes()).flatten().collect();
        v.extend_from_slice(message);
        let mut hasher = Hasher::new(HASHES).hash_extendable(&v);
        let mut challenges = vec![0u8; (ROUNDS * 4) as usize];
        hasher.read_exact(&mut challenges).map_err(|_| Error::new())?;

        let mut ephemeral_cge: Vec<ClassGroupElement> = Vec::new();
        let mut opened_curve_indices: Vec<u32> = Vec::new();
        let mut opened_curves: Vec<MontgomeryCurve> = Vec::new();

        for (challenge_bytes, ephem_cge) in challenges.chunks_exact(4).zip(b) {
            let n = <i32>::from_be_bytes(challenge_bytes.try_into().unwrap());
            // this is uniform because CURVES is a power of two
            let curve_num = (n.unsigned_abs() % CURVES) as usize;
            let s = if n > 0 {
                let j = ephem_cge - self.secret_actions[curve_num].clone();
                // println!("{}", j.reduce().variable_time_action(&self.public_curves[curve_num]).a.x- ec.a.x);
                j
            } else {
                let j = ephem_cge + self.secret_actions[curve_num].clone();
                // println!("T: {}", j.reduce().variable_time_action(&self.public_curves[curve_num].twist()).a.x-  ec.a.x);
                j
            };
            opened_curves.push(self.public_curves[curve_num].clone());
            ephemeral_cge.push(s);
            opened_curve_indices.push(curve_num as u32);
        }
        Ok(Signature {
            num_curves: CURVES,
            challenges,
            ephemeral_cge,
            opened_curves,
            proof: self.proof_tree.proof_from_leaf_indices(&opened_curve_indices),
        })
    }
}

impl<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> Verifier<Signature<CURVES, ROUNDS, HASHES>> for VerifyingKey {
    fn verify(&self, msg: &[u8], signature: &Signature<CURVES, ROUNDS, HASHES>) -> Result<(), Error> {
        if signature.challenges.len() % 4 != 0 {
            return Err(Error::new());
        }
        let challenges = signature.challenges.chunks_exact(4).map(|chunk| i32::from_be_bytes(chunk.try_into().unwrap())).collect::<Vec<i32>>();

        if challenges.len() != signature.opened_curves.len() {
            return Err(Error::new());
        }

        let leaf_hashes = challenges.iter().zip(&signature.opened_curves).map(|(challenge, curve)| {
            let label = (challenge.unsigned_abs() % CURVES) + signature.num_curves;
            let v = [
                curve.a.x.get_be_bytes().to_vec(),
                label.to_be_bytes().to_vec(),
                self.merkle_key.to_vec(),
            ].concat();
            Ok((label, Hasher::new(HASHES).hash(&v)))
        }).collect::<Result<Vec<(u32, HashType)>, Error>>()?;
        signature.proof.verify(&self.root, leaf_hashes, &self.merkle_key).map_err(|_e| Error::new()).unwrap();

        let ephemeral_curves = signature.ephemeral_cge.iter().zip(&signature.opened_curves).zip(&challenges).map(|((ri, curve), challenge)| {
            Ok(if *challenge > 0 {
                ri.reduce().variable_time_action(curve)
            } else {
                ri.reduce().variable_time_action(&curve.twist())
            })
        }).collect::<Result<Vec<MontgomeryCurve>, Error>>()?;
        let mut v: Vec<u8> = ephemeral_curves.into_iter().map(|x| x.to_be_bytes()).flatten().collect();
        v.extend_from_slice(msg);
        let mut hasher = Hasher::new(HASHES).hash_extendable(&v);
        let mut derived_challenges_bytes = vec![0u8; (ROUNDS * 4) as usize];
        hasher.read_exact(&mut derived_challenges_bytes).map_err(|_| Error::new())?;
        let net_diff = derived_challenges_bytes.chunks_exact(4).zip(&challenges).map(|(x, c)| <i32>::from_be_bytes(x.try_into().unwrap()) - c).reduce(|x, acc| acc + x);
        match net_diff {
            Some(0) => { Ok(()) },
            _ => Err(Error::new())
        }
    }
}


mod tests {
    use rand::{RngCore, thread_rng};

    use super::*;

    #[test]
    fn verify() {
        let mut msg = [0u8; 1024];
        thread_rng().fill_bytes(&mut msg);
        let j = SigningKey::<256, 7, 11>::generate();
        // let signature = j.try_sign(&msg).unwrap();
        // j.verifying_key().verify(&msg, &signature).unwrap();
        for _ in 0..10 {
            let signature = j.try_sign(&msg).unwrap();
            j.verifying_key().verify(&msg, &signature).unwrap();
        }
    }
}
