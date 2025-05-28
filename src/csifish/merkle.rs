use std::collections::{HashMap, HashSet, VecDeque};

use rand::{RngCore, thread_rng};
use sha3::Digest;
use crate::csifish::field_arithmetic::arithmetic::ModularArithmetic;

use crate::csifish::hash::{Hasher, HashType};
use crate::csifish::montgomery::MontgomeryCurve;
use crate::csifish::constants::VerificationFailed;

#[derive(Debug, Clone, PartialEq)]
pub struct ClassGroupMerkleTree<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> {
    root: HashType,
    merkle_key: HashType,
    layers: Vec<Vec<HashType>>,
}

impl<const CURVES: u32, const ROUNDS: u32, const HASHES: u32> ClassGroupMerkleTree<CURVES, ROUNDS, HASHES> {
    pub fn from_leaves(leaves: &[MontgomeryCurve]) -> Self {
        let mut merkle_key = HashType::default();
        thread_rng().fill_bytes(&mut merkle_key);

        let curves = leaves.len() as u32;
        assert_eq!(CURVES, curves);
        let depth = leaves.len().ilog2();

        let mut layers: Vec<Vec<HashType>> = Vec::new();
        layers.push(
            leaves.iter().enumerate().map(|(i, curve)| {
                let label = curves + i as u32;
                let v = [
                    curve.a.x.get_be_bytes().to_vec(),
                    label.to_be_bytes().to_vec(),
                    merkle_key.to_vec(),
                ].concat();
                let result: HashType = Hasher::new(HASHES).hash(&v);
                result
            }).collect::<Vec<HashType>>(),
        );
        for _ in 0..depth {
            layers.push(
                layers[layers.len() - 1].chunks(2).enumerate().map(|(i, chunk)| {
                    let label = (curves >> (layers.len())) + i as u32;
                    //hash left and right sibling, label,
                    let v = [
                        &chunk[0][..],
                        &chunk[1][..],
                        &label.to_be_bytes()[..],
                        &merkle_key[..],
                    ].concat();
                    let result: HashType = Hasher::new(HASHES).hash(&v);
                    result
                }).collect::<Vec<HashType>>(),
            );
        }
        ClassGroupMerkleTree {
            root: layers.last().unwrap()[0],
            merkle_key,
            layers,
        }
    }

    pub fn proof_from_leaf_indices(&self, leaf_indices: &[u32]) -> ClassGroupMerkleProof {
        let mut level: Vec<u32> = leaf_indices.iter().map(|x| x + CURVES).collect();
        let mut unknown = HashSet::<u32>::new();
        let mut known = HashSet::<u32>::new();
        for _ in 0..self.depth() {
            known.extend(level.iter().clone());
            let mut next_level = Vec::new();
            for idx in level {
                let is_odd = idx % 2;
                let (sibling, parent) = (idx + 1 - 2 * (is_odd), idx / 2);
                if !known.contains(&sibling) {
                    unknown.insert(sibling);
                }
                next_level.push(parent);
            }
            level = next_level;
        }
        let proof_indices = unknown.into_iter().collect::<Vec<u32>>();
        let mut proof: Vec<(u32, HashType)> = Vec::new();
        for idx in proof_indices {
            let level = idx.ilog2();
            let pos_in_level = idx - (1 << level);
            proof.push((idx, self.layers[self.depth() - level as usize][pos_in_level as usize]))
        }
        ClassGroupMerkleProof {
            num_hashes: HASHES,
            proof,
        }
    }
    pub fn depth(&self) -> usize {
        self.layers.len() - 1
    }
    pub fn root(&self) -> HashType {
        self.root
    }
    pub fn merkle_key(&self) -> HashType {
        self.merkle_key
    }
    pub fn leaves(&self) -> Vec<HashType> {
        self.layers[0].clone()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ClassGroupMerkleProof {
    num_hashes: u32,
    proof: Vec<(u32, HashType)>,
}

impl ClassGroupMerkleProof {
    pub fn verify(
        &self,
        root: &HashType,
        mut leaf_hashes: Vec<(u32, HashType)>,
        merkle_key: &HashType,
    ) -> Result<(), VerificationFailed> {
        leaf_hashes.sort_by_key(|(a, _)| *a);
        leaf_hashes.dedup();
        let mut level: VecDeque<(u32, HashType)> = VecDeque::from(leaf_hashes);
        let mut tree: HashMap<u32, HashType> = HashMap::from_iter((&self.proof).clone());
        loop {
            assert!(level.len() > 0);
            let (label, hash) = level.pop_front().unwrap();
            if tree.contains_key(&(label / 2)) {
                continue
            }
            let hasher = Hasher::new(self.num_hashes);
            let result = if (label % 2) == 0 && level.front().is_some_and(|(i, _)| (label + 1) == *i) {
                Ok(hasher.hash(&[&hash[..], &level.front().unwrap().1[..], &(label / 2).to_be_bytes()[..], &merkle_key[..]].concat()))
            } else {
                match tree.get(&(label + 1 - (2 * (label % 2)))) {
                    Some(in_proof) => Ok(if (label % 2) == 0 {
                        hasher.hash(&[&hash[..], &in_proof[..], &(label / 2).to_be_bytes()[..], &merkle_key[..]].concat())
                    } else {
                        hasher.hash(&[&in_proof[..], &hash[..], &(label / 2).to_be_bytes()[..], &merkle_key[..]].concat())
                    }),
                    None => { println!("Cur: {} Des: {}", label, label + 1 - (2 * (label % 2))); Err(VerificationFailed) }
                }
            }?;
            let next = (label / 2, result);
            level.push_back(next);
            tree.insert(next.0, next.1);
            if next.0 == 1 {
                break;
            }
        }
        let putative_root = tree.get(&1).unwrap();
        if putative_root == root {
            Ok(())
        } else {
            Err(VerificationFailed)
        }
    }
}

mod tests {
    use crate::csifish::field_arithmetic::base_field::FieldElement;
    use crate::csifish::field_arithmetic::arithmetic::{ModularArithmetic};

    use super::*;

    #[test]
    fn test() {
        let j: Vec<MontgomeryCurve> = (0..16).map(|_| MontgomeryCurve::new(FieldElement::random(&mut thread_rng()))).collect();
        let mt = ClassGroupMerkleTree::<16, 7, 12>::from_leaves(&j);
        let proof = mt.proof_from_leaf_indices(&[0, 3, 14]);
        // assert_eq!(proof.proof[0].1,
        //            Hasher::new(12).hash(&[j[1].a.x.retrieve_to_be_bytes().to_vec(), 5u32.to_be_bytes().to_vec(), mt.merkle_key.to_vec()].concat())
        // );
        // assert_eq!(proof.proof[1].1,
        //            Hasher::new(12).hash(&[j[2].a.x.retrieve_to_be_bytes().to_vec(), 6u32.to_be_bytes().to_vec(), mt.merkle_key.to_vec()].concat())
        // );
        // let mut serialized = proof.serialize();
        // let proof = ClassGroupMerkleProof::deserialize(serialized.as_slice()).unwrap();
        let result = proof.verify(
            &mt.root(),
            Vec::from([(16, mt.leaves()[0]), (19, mt.leaves()[3]), (30, mt.leaves()[14])]),
            &mt.merkle_key,
        );
        result.unwrap();
    }

    // #[test]
    // fn serialize_merkle_tree() {
    //     let j: Vec<MontgomeryCurve> = (0..256).map(|_| MontgomeryCurve::new(FieldElement::random(&mut thread_rng()))).collect();
    //     let mt = ClassGroupMerkleTree::<256, 7, 12>::from_leaves(&j);
    //     let serialized = mt.serialize();
    //     let deserialized = ClassGroupMerkleTree::deserialize(&serialized).unwrap();
    //     assert_eq!(mt, deserialized);
    // }
    //
    // #[test]
    // #[should_panic]
    // fn serialize_failure() {
    //     let j: Vec<MontgomeryCurve> = (0..256).map(|_| MontgomeryCurve::new(FieldElement::random(&mut thread_rng()))).collect();
    //     let mt = ClassGroupMerkleTree::<256, 7, 12>::from_leaves(&j);
    //     let serialized = mt.serialize();
    //     ClassGroupMerkleTree::<256, 7, 12>::deserialize(&serialized[..serialized.len() - 1]).unwrap();
    // }
}
