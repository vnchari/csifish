// #![allow(dead_code)]
//
// use criterion::{Criterion, criterion_group, criterion_main};
// use crypto_bigint::{impl_modulus, Random, U512};
// use crypto_bigint::modular::{ConstMontyForm, ConstMontyParams};
// use crypto_bigint::rand_core::RngCore;
// use rand::thread_rng;
//
// use csifish::csifish::classgroup::ReducedClassGroupElement;
// use csifish::csifish::montgomery::{MontgomeryCurve, Point};
// use csifish::csifish::signatures::KeyPair;
// use csifish::field::field::FieldElement;
// use csifish::precomputed::constants::BASE_CURVE;
//
// type Bencher = Criterion;
//
//
// fn bench_add(c: &mut Bencher) {
//     let p = FieldElement::random(&mut thread_rng());
//     let q = FieldElement::random(&mut thread_rng());
//     c.bench_function("bench_add", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p + q
//             );
//         });
//     });
// }
//
// fn bench_sub(c: &mut Bencher) {
//     let p = FieldElement::random(&mut thread_rng());
//     let q = FieldElement::random(&mut thread_rng());
//     c.bench_function("bench_sub", |b| {
//         b.iter(|| {
//             criterion::black_box(
//                 p - q
//             );
//         });
//     });
// }
//
//
// impl_modulus!(ModulusP, U512, "65b48e8f740f89bffc8ab0d15e3e4c4ab42d083aedc88c425afbfcc69322c9cda7aac6c567f35507516730cc1f0b4f25c2721bf457aca8351b81b90533c6c87b");
// pub type ModP = ConstMontyForm<ModulusP, { ModulusP::LIMBS }>;
//
// fn bench_mul(c: &mut Bencher) {
//     let p = ModP::random(&mut thread_rng());
//     let q = ModP::random(&mut thread_rng());
//     let p1 = FieldElement::from_montgomery(p.as_montgomery().clone().to_words());
//     let q1 = FieldElement::from_montgomery(q.as_montgomery().clone().to_words());
//     c.bench_function("bench_mul", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p1 * q1
//             );
//         });
//     });
// }
//
// fn bench_mul_bigint(c: &mut Bencher) {
//     let p = ModP::random(&mut thread_rng());
//     let q = ModP::random(&mut thread_rng());
//     c.bench_function("bench_mul_bigint", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p * q
//             );
//         });
//     });
// }
//
//
// fn bench_sqr_bigint(c: &mut Bencher) {
//     let p = ModP::random(&mut thread_rng());
//     c.bench_function("bench_sqr_bigint", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p.square()
//             );
//         });
//     });
// }
//
// fn bench_sqr(c: &mut Bencher) {
//     let p = ModP::random(&mut thread_rng());
//     let p1 = FieldElement::from_montgomery(p.as_montgomery().clone().to_words());
//     c.bench_function("bench_sqr", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p1.square()
//             );
//         });
//     });
// }
//
//
// fn bench_inv(c: &mut Bencher) {
//     let mut p = ModP::random(&mut thread_rng());
//     let p1 = FieldElement::from(p.retrieve().to_words());
//     c.bench_function("bench_inv", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p1.inv()
//             );
//         });
//     });
// }
//
// fn bench_inv_bigint(c: &mut Bencher) {
//     let p = ModP::random(&mut thread_rng());
//     c.bench_function("bench_inv_bigint", |b| {
//         b.iter(|| {
//             std::hint::black_box(
//                 p.inv()
//             );
//         });
//     });
// }
//
// fn bench_double_curve(c: &mut Bencher) {
//     let mut p = Point::random();
//     while !BASE_CURVE.on_curve(&p) {
//         p = Point::random();
//     }
//     c.bench_function("bench_double_curve", |b| {
//         b.iter(|| {
//             std::hint::black_box(BASE_CURVE.double(&p));
//         });
//     });
// }
//
// fn bench_isogeny(c: &mut Bencher) {
//     let a = FieldElement::from_be_hex("48912381B13014DF7E10F242424DFE6D43860ED48A2913843A45E75E15615849B2E2C8191E6CEF70A931E20883E8B59B87046926B8E534DCA88722A8E204496C");
//     let order: usize = 173usize;
//     let k = Point::from_x(FieldElement::from_be_hex("54C8CDA4F5B40B4DD5EF9011AFA313195A68106114B157B53270CB1005C8338E4CE00C826ECEE406027F383FA1D5037DBB81D92E4203B4092B9C3D20A32D49A8"));
//     let p = Point::from_x(FieldElement::from_be_hex("55E59A6BB770F1477F38E747D4C45F61CDA4D068736398DCB7C3A6B872208E6BA55FA42377A4B3EB25AEF4D0CE59C91A1D3A291B87700FFFE21805A7DEED199B"));
//     let ea = MontgomeryCurve::new(a);
//     c.bench_function("bench_isogeny", |b| {
//         b.iter(|| {
//             std::hint::black_box(ea.isogeny(&k, order, &p));
//         });
//     });
// }
//
// fn bench_action(c: &mut Bencher) {
//     let e = MontgomeryCurve::new(FieldElement::ZERO);
//     let exp: [i8; 74] = [
//         -5, 2, 0, -3, 4, -4, -5, 3, 5, -1, -2, -4, 0, -2, -3, 3, 1, -2, 5, 3, 4, 3, -4, 2, 2,
//         3, -1, 0, 1, -3, 0, 1, -5, -2, 0, 2, 0, 0, -5, 5, 4, 5, 0, -5, 0, -1, 0, 1, 5, 1, 1,
//         -3, 0, 5, 1, 2, -1, 1, -5, 0, 1, 5, 3, 2, -1, -5, 4, 2, 1, 2, -2, 0, 1, 5,
//     ];
//
//     let g1 = ReducedClassGroupElement::new(exp);
//     c.bench_function("bench_action", |b| {
//         b.iter(|| {
//             std::hint::black_box(g1.act_on(&e));
//         });
//     });
// }
//
//
// fn bench_keygen_256_13(c: &mut Criterion) {
//     c.bench_function("bench_keygen_256_13", |b| {
//         b.iter(|| {
//             std::hint::black_box({ KeyPair::<256, 13, 1>::generate() });
//         });
//     });
// }
//
// fn bench_sign_256_13(c: &mut Criterion) {
//     let keypair = KeyPair::<256, 13, 1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     c.bench_function("bench_sign_256_13", |b| {
//         b.iter(|| {
//             std::hint::black_box(keypair.sign(&msg));
//         });
//     });
// }
//
// fn bench_verify_256_13(c: &mut Criterion) {
//     let keypair = KeyPair::<256, 13, 1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     let signature = keypair.sign(&msg);
//     let pk = keypair.public_key();
//     c.bench_function("bench_verify_256_13", |b| {
//         b.iter(|| {
//             std::hint::black_box(pk.verify(&signature, &msg));
//         });
//     });
// }
//
// fn bench_keygen_1024_11(c: &mut Criterion) {
//     c.bench_function("bench_keygen_1024_11", |b| {
//         b.iter(|| {
//             std::hint::black_box({ KeyPair::<1024, 11,1>::generate() });
//         });
//     });
// }
//
// fn bench_sign_1024_11(c: &mut Criterion) {
//     let keypair = KeyPair::<1024, 11,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     c.bench_function("bench_sign_1024_11", |b| {
//         b.iter(|| {
//             std::hint::black_box({ keypair.sign(&msg) });
//         });
//     });
// }
//
// fn bench_verify_1024_11(c: &mut Criterion) {
//     let keypair = KeyPair::<1024, 11,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     let signature = keypair.sign(&msg);
//     let pk = keypair.public_key();
//     c.bench_function("bench_verify_1024_11", |b| {
//         b.iter(|| {
//             std::hint::black_box({ pk.verify(&signature, &msg) });
//         });
//     });
// }
//
// fn bench_keygen_4096_9(c: &mut Criterion) {
//     c.bench_function("bench_keygen_4096_9", |b| {
//         b.iter(|| {
//             std::hint::black_box({ KeyPair::<4096, 9,1>::generate() });
//         });
//     });
// }
//
// fn bench_sign_4096_9(c: &mut Criterion) {
//     let keypair = KeyPair::<4096, 9,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     c.bench_function("bench_sign_4096_9", |b| {
//         b.iter(|| {
//             std::hint::black_box({ keypair.sign(&msg) });
//         });
//     });
// }
//
// fn bench_verify_4096_9(c: &mut Criterion) {
//     let keypair = KeyPair::<4096, 9,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     let signature = keypair.sign(&msg);
//     let pk = keypair.public_key();
//     c.bench_function("bench_verify_4096_9", |b| {
//         b.iter(|| {
//             std::hint::black_box({ pk.verify(&signature, &msg) });
//         });
//     });
// }
//
// fn bench_keygen_32768_7(c: &mut Criterion) {
//     c.bench_function("bench_keygen_32768_7", |b| {
//         b.iter(|| {
//             std::hint::black_box({ KeyPair::<32768, 7,1>::generate() });
//         });
//     });
// }
//
// fn bench_sign_32768_7(c: &mut Criterion) {
//     let keypair = KeyPair::<32768, 7,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     c.bench_function("bench_sign_32768_7", |b| {
//         b.iter(|| {
//             std::hint::black_box({ keypair.sign(&msg) });
//         });
//     });
// }
//
// fn bench_verify_32768_7(c: &mut Criterion) {
//     let keypair = KeyPair::<32768, 7,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     let signature = keypair.sign(&msg);
//     let pk = keypair.public_key();
//     c.bench_function("bench_verify_32768_7", |b| {
//         b.iter(|| {
//             std::hint::black_box({ pk.verify(&signature, &msg) });
//         });
//     });
// }
//
// fn bench_keygen_262144_6(c: &mut Criterion) {
//     c.bench_function("bench_keygen_262144_6", |b| {
//         b.iter(|| {
//             std::hint::black_box({ KeyPair::<262144, 6,1>::generate() });
//         });
//     });
// }
//
// fn bench_sign_262144_6(c: &mut Criterion) {
//     let keypair = KeyPair::<262144, 6,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     c.bench_function("bench_sign_262144_6", |b| {
//         b.iter(|| {
//             std::hint::black_box({ keypair.sign(&msg) });
//         });
//     });
// }
//
// fn bench_verify_262144_6(c: &mut Criterion) {
//     let keypair = KeyPair::<262144, 6,1>::generate();
//     let mut msg = [0u8; 1024];
//     thread_rng().fill_bytes(&mut msg);
//     let signature = keypair.sign(&msg);
//     let pk = keypair.public_key();
//     c.bench_function("bench_verify_262144_6", |b| {
//         b.iter(|| {
//             std::hint::black_box({ pk.verify(&signature, &msg) });
//         });
//     });
// }
//
//
//
//
//
// criterion_group!(name = benches;
//     config = Criterion::default();//.sample_size(50);
//     targets =
//         // bench_double_curve,
//         // bench_isogeny,
//         // bench_keygen_256_13,
//         bench_sign_256_13,
//         bench_verify_256_13,
//         // bench_sign_1024_11,
//         // bench_verify_1024_11
//         // bench_action
//         // bench_add,
//         // bench_sub,
//         // bench_inv,
//         // bench_inv_bigint,
//         // bench_mul,
//         // bench_sqr,
//         // bench_mul_bigint,
//         // bench_sqr_bigint
//         // bench_mul_bigint
// );
// criterion_main!(benches);
