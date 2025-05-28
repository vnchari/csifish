use rug::{Float, Integer};
use rug::integer::Order;
use rug::float::Round;
use rand::{Rng, thread_rng};

use crate::csifish::constants::{BASIS, FLOAT_PRECISION, NUM_PRIMES, ORTHO_BASIS, ORTHO_NORMS, POOL};
use crate::csifish::field_arithmetic::classgroup::ClassGroupElement;

pub fn dot(b: &[Integer; 74], basis_idx: usize) -> Float {
    let slice = &ORTHO_BASIS[NUM_PRIMES * basis_idx..NUM_PRIMES * (basis_idx + 1)];
    slice
        .iter()
        .zip(b)
        .fold(Float::with_val(FLOAT_PRECISION, 0.0), |accum, (s, b)| {
            accum + Float::with_val(FLOAT_PRECISION, s * b)
        })
}

pub fn sub_slice(a: &[i8], b: &[i8]) -> [i8; 74] {
    let mut result = [0i8; 74];
    for idx in 0..result.len() {
        assert!(idx < a.len());
        assert!(idx < b.len());
        result[idx] = a[idx] - b[idx];
    }
    result
}

pub fn add_slice(a: &[i8], b: &[i8]) -> [i8; 74] {
    let mut result = [0i8; 74];
    for idx in 0..result.len() {
        assert!(idx < a.len());
        assert!(idx < b.len());
        result[idx] = a[idx] + b[idx];
    }
    result
}

pub fn l1(a: &[i8; 74]) -> u16 {
    a.iter().fold(0u16, |s, x| s + (x.unsigned_abs() as u16))
}

pub fn dlw_reduce(e: [i8; 74], pool_size: usize) -> [i8; 74] {
    let mut e_prime = e;
    let mut stalled = false;
    let mut best_norm = l1(&e_prime);
    while !stalled {
        stalled = true;
        for idx in 0..pool_size {
            let s: &[i8] = &POOL[idx * NUM_PRIMES..idx * NUM_PRIMES + 74];
            let diff = &sub_slice(&e_prime, s);
            let sum = &add_slice(&e_prime, s);
            let l1diff = l1(diff);
            let l1sum = l1(sum);
            if l1sum < best_norm {
                best_norm = l1sum;
                e_prime = *sum;
                stalled = false;
            }
            if l1diff < best_norm {
                best_norm = l1diff;
                e_prime = *diff;
                stalled = false;
            }
        }
    }
    e_prime
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReducedClassGroupElement {
    // exponent for each ideal, separated into positive/negative exponents
    // exponents[0][i] is positive exponent for ideal i
    pub(crate) exponents: [i8; NUM_PRIMES],
    // // (p+1)/k, where k is product of primes which have non-zero exponent
    // // in the corresponding positive/negative_exponents
}

impl ClassGroupElement {
    pub fn reduce(&self) -> ReducedClassGroupElement {
        let pool_size = 7500usize;
        let mut b: [Integer; NUM_PRIMES] = [Integer::ZERO; 74];
        let p_mp: Integer = Integer::from_digits(&self.limbs.map(|x| x.to_le()), Order::LsfLe);
        b[0].clone_from(&p_mp);
        for basis_idx in (0..NUM_PRIMES).rev() {
            let (c, _) = (dot(&b, basis_idx) / &ORTHO_NORMS[basis_idx])
                .to_integer_round(Round::Nearest)
                .unwrap();
            let slice = &BASIS[basis_idx * NUM_PRIMES..(basis_idx + 1) * NUM_PRIMES];
            for dim_idx in 0..b.len() {
                b[dim_idx] -= &c * slice[dim_idx];
            }
        }

        let mut e_prime = dlw_reduce(b.map(|x| x.to_i8().unwrap()), pool_size);
        let mut best_len = l1(&e_prime);

        for _ in 0..2 {
            let shifted = {
                let ridx = thread_rng().gen_range(0..10000);
                let ridx2 = thread_rng().gen_range(0..10000);
                add_slice(
                    &add_slice(
                        &e_prime,
                        POOL[ridx * NUM_PRIMES..ridx * NUM_PRIMES + NUM_PRIMES]
                            .try_into()
                            .unwrap(),
                    ),
                    POOL[ridx2 * NUM_PRIMES..ridx2 * NUM_PRIMES + NUM_PRIMES]
                        .try_into()
                        .unwrap(),
                )
            };
            let t = dlw_reduce(shifted, pool_size);
            let norm_t = l1(&t);
            if norm_t < best_len {
                best_len = norm_t;
                e_prime = t;
            }
        }
        ReducedClassGroupElement::new(e_prime)
    }

    pub fn reduce_one_round(&self) -> ReducedClassGroupElement {
        let pool_size = 7500usize;
        let mut b: [Integer; NUM_PRIMES] = [Integer::ZERO; 74];
        let p_mp: Integer = Integer::from_digits(&self.limbs.map(|x| x.to_le()), Order::LsfLe);
        b[0].clone_from(&p_mp);
        for basis_idx in (0..NUM_PRIMES).rev() {
            let (c, _) = (dot(&b, basis_idx) / &ORTHO_NORMS[basis_idx])
                .to_integer_round(Round::Nearest)
                .unwrap();
            let slice = &BASIS[basis_idx * NUM_PRIMES..(basis_idx + 1) * NUM_PRIMES];
            for dim_idx in 0..b.len() {
                b[dim_idx] -= &c * slice[dim_idx];
            }
        }
        // ReducedClassGroupElement::new(b.map(|x| x.to_i8().unwrap()))
        ReducedClassGroupElement::new(dlw_reduce(b.map(|x| x.to_i8().unwrap()), pool_size))
    }
}

impl ReducedClassGroupElement {
    pub fn new(exponents: [i8; NUM_PRIMES]) -> ReducedClassGroupElement {
        ReducedClassGroupElement { exponents }
    }
}

#[cfg(test)]
mod tests {
    use rand::thread_rng;
    use crate::csifish::field_arithmetic::arithmetic::ModularArithmetic;
    use crate::csifish::variable_time::VariableTimeAction;
    use crate::csifish::constants::BASE_CURVE;
    use crate::csifish::field_arithmetic::classgroup::ClassGroupElement;

    use crate::csifish::lattice::ReducedClassGroupElement;

    #[test]
    fn reduce_basic() {
        let mut one = [0i8; 74];
        one[0] = 1;
        let mut el = ClassGroupElement::random(&mut thread_rng());
        let res1 = el.reduce_one_round();
        // for i in res1.exponents {
        //     print!("{i} ")
        // }
        println!(
            "{}",
            ReducedClassGroupElement::new(one)
                .variable_time_action(&res1.variable_time_action(&BASE_CURVE))
                .a
        );
        el.limbs[0] += 1;
        let res2 = el.reduce_one_round();
        // for i in res2.exponents {
        //     print!("{i} ")
        // }
        println!("{}", res2.variable_time_action(&BASE_CURVE).a);

        // const N: usize = 1_000;
        // let index_distr: Vec<usize> = (0..N).par_bridge().map(|_| ClassGroupElement::random().reduce_one_round()).collect::<Vec<ReducedClassGroupElement>>().iter().map(
        //     |j| {
        //         // j.exponents.iter().enumerate().for_each(|(i, x)| {
        //         //     let xabs = x.unsigned_abs() as u64;
        //         //     sum[i] += xabs;
        //         //     if xabs > max[i] { max[i] = xabs; }
        //         // });
        //         let mut idx: usize = 0;
        //         let mut best_val: u8 = 0;
        //         for (i, e) in j.exponents.iter().enumerate() {
        //             if e.unsigned_abs() > best_val {
        //                 idx = i;
        //                 best_val = e.unsigned_abs();
        //             }
        //         }
        //         idx
        //     }
        // ).collect();
        // for i in index_distr {
        //     print!("{i},");
        // }
        // let mut vars = [0f64; NUM_PRIMES];
        //
        // let mut sums = [0u64; NUM_PRIMES];
        // let distrs: Vec<[u64; NUM_PRIMES]> = (0..N).par_bridge().map(|_| ClassGroupElement::random().reduce_one_round()).collect::<Vec<ReducedClassGroupElement>>().iter().map(
        //     |j| {
        //         // j.exponents.iter().enumerate().for_each(|(i, x)| {
        //         //     let xabs = x.unsigned_abs() as u64;
        //         //     sum[i] += xabs;
        //         //     if xabs > max[i] { max[i] = xabs; }
        //         // });
        //         let mut j = j.exponents.map(|x| x.unsigned_abs() as u64);
        //         j.sort();
        //         j
        //     }
        // ).collect();
        // distrs.clone().into_iter().for_each(|x| {
        //     for i in 0..NUM_PRIMES {
        //         sums[i] += x[i];
        //     }
        // });
        // let means = sums.map(|x| x as f64 / N as f64);
        // distrs.clone().into_iter().for_each(|x| {
        //     // let mut q = x.clone();
        //     for i in 0..NUM_PRIMES {
        //         vars[i] += (x[i] as f64 - means[i]) * (x[i] as f64 - means[i]);
        //     }
        // });
        // let std_devs = vars.map(|x| (x / N as f64).sqrt());
        //
        // println!();
        // for i in means {
        //     print!("{i},");
        // }
        // println!();
        // for i in std_devs {
        //     print!("{i},");
        // }
        // let maxes = distrs.clone().into_iter().map(|x| x[NUM_PRIMES - 1]);
        // for max in maxes {
        //     print!("{max},");
        // }
        // println!();
        // println!("{}", sum.iter().sum::<u64>() as f64 / N as f64)
    }
}
