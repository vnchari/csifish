use crypto_bigint::rand_core::{CryptoRng, RngCore};
use rand::thread_rng;
use rug::integer::Order::LsfLe;
use rug::Integer;

use crate::csifish::montgomery::{MontgomeryCurve, Point};
use crate::csifish::variable_time::{VariableTimeAction, VariableTimeCurve};
use crate::csifish::field_arithmetic::arithmetic::{ModularArithmetic, MontgomeryArithmetic};
use crate::csifish::field_arithmetic::base_field::FieldElement;
use crate::csifish::constants::{NUM_PRIMES, P_GMP, PRIMES16};
use crate::csifish::lattice::ReducedClassGroupElement;

pub trait OneTimeCurve {
    fn elligator(&self, rng: &mut (impl CryptoRng + RngCore)) -> (Point, Point);
    fn two_point_isogeny(
        &self,
        k: &Point,
        ell: usize,
        p1: &Point,
        p2: &Point,
    ) -> (Point, Point, MontgomeryCurve);
}

pub trait OneTimeAction {
    //consumes self on purpose
    fn one_time_blinded_action(self, e: &MontgomeryCurve) -> MontgomeryCurve;
}

impl OneTimeCurve for MontgomeryCurve {
    fn elligator(&self, rng: &mut (impl CryptoRng + RngCore)) -> (Point, Point) {
        loop {
            let u = FieldElement::random_under_half(rng);
            let u2 = u.square();
            let mut d = u2 - FieldElement::ONE;
            if u.is_zero() || d.is_zero() {
                continue;
            }
            let mut m = self.a.x * u2;
            let mut t = self.a.x * m;
            let mut p = self.a.x.clone();
            let is_base = self.a.x.is_zero() as u64;
            p.conditional_move(is_base, &FieldElement::ONE);
            m.conditional_move(is_base, &FieldElement::ONE);
            t.conditional_move(is_base, &FieldElement::ONE);
            d *= self.a.z;
            t += d.square();
            t *= d;
            t *= p;

            let mneg = m.neg();
            let mut p_plus = Point { x: p, z: d };
            let mut p_minus = Point { x: mneg, z: d };
            //TODO: leaky
            let leg = (Integer::from_digits(&t.get_montgomery().map(|x| x.to_le()), LsfLe)
                .legendre(&P_GMP)
                >> 31) as u64
                & 1;
            p_plus.x.conditional_move(leg, &mneg);
            p_minus.x.conditional_move(leg, &p);
            return (p_plus, p_minus);
        }
    }

    fn two_point_isogeny(
        &self,
        k: &Point,
        ell: usize,
        p1: &Point,
        p2: &Point,
    ) -> (Point, Point, MontgomeryCurve) {
        //compute twisted Edwards curve coefficients
        let mut edwards_z = self.a.z + self.a.z;
        let mut edwards_x = self.a.x + edwards_z;
        edwards_z = self.a.x - edwards_z;

        let p1add: FieldElement = p1.x + p1.z;
        let p1sub: FieldElement = p1.x - p1.z;

        let p2add: FieldElement = p2.x + p2.z;
        let p2sub: FieldElement = p2.x - p2.z;

        let mut prod = Point {
            x: k.x - k.z,
            z: k.x + k.z,
        };

        let tmp1 = prod.x * p1add;
        let tmp0 = prod.z * p1sub;
        let mut q1 = Point {
            x: tmp0 + tmp1,
            z: tmp0 - tmp1,
        };
        let tmp1 = prod.x * p2add;
        let tmp0 = prod.z * p2sub;
        let mut q2 = Point {
            x: tmp0 + tmp1,
            z: tmp0 - tmp1,
        };

        // precompute a24.x = A.x+2*A.z, a24.z = 4*A.z
        let t = self.a.z + self.a.z;
        let a24 = Point {
            x: self.a.x + t,
            z: t + t,
        };

        let p_2 = Point {
            x: (k.x + k.z).square(),
            z: (k.x - k.z).square(),
        };
        let c = p_2.x - p_2.z;
        let b = a24.z * p_2.z;
        let a = (c * a24.x) + b;
        let mut kernel_buffer = [
            k.clone(),
            Point {
                x: p_2.x * b,
                z: a * c,
            },
            Point::zero(),
        ];

        for i in 1..(ell / 2) {
            let cur = kernel_buffer[i % 3];
            let tmp1 = cur.x - cur.z;
            let tmp0 = cur.x + cur.z;
            prod.x *= tmp1;
            prod.z *= tmp0;

            let t1 = tmp1 * p1add;
            let t0 = tmp0 * p1sub;
            q1.x *= t0 + t1;
            q1.z *= t0 - t1;

            let t1 = tmp1 * p2add;
            let t0 = tmp0 * p2sub;
            q2.x *= t0 + t1;
            q2.z *= t0 - t1;

            kernel_buffer[(i + 1) % 3] =
                self.differential_add(&cur, k, &kernel_buffer[(i - 1) % 3]);
        }
        q1.x = q1.x.square() * p1.x;
        q1.z = q1.z.square() * p1.z;

        q2.x = q2.x.square() * p2.x;
        q2.z = q2.z.square() * p2.z;

        edwards_x = edwards_x.constant_time_bounded_exp(&(ell as u64));
        edwards_z = edwards_z.constant_time_bounded_exp(&(ell as u64));
        edwards_x *= prod.z.square().square().square();
        edwards_z *= prod.x.square().square().square();
        let mut ax = edwards_x + edwards_z;
        ax += ax;
        let az = edwards_x - edwards_z;
        let codomain = MontgomeryCurve::projective(ax, az);
        (q1, q2, codomain)
    }
}

impl OneTimeAction for ReducedClassGroupElement {
    fn one_time_blinded_action(self, e: &MontgomeryCurve) -> MontgomeryCurve {
        const NUM_BATCHES: usize = 4;
        const MERGE_AFTER: usize = 2;
        const BLIND_MAX_EXP: u8 = 2;
        let mut blinded_exponents = self.exponents;
        let mut blinding = [0i8; NUM_PRIMES];
        // log_2(5^74) ~= 2^(178)
        for i in 0..NUM_PRIMES {
            let b = loop {
                let mut tmp = [0u8; 1];
                thread_rng().fill_bytes(&mut tmp);
                let val = tmp[0] >> 5;
                if val <= 2 * BLIND_MAX_EXP {
                    break val.wrapping_sub(BLIND_MAX_EXP) as i8;
                }
            };
            blinded_exponents[i] += b;
            blinding[i] -= b;
        }
        let mut e = ReducedClassGroupElement::new(blinded_exponents).variable_time_action(e);
        let mut isogeny_count = [2u8; NUM_PRIMES];
        let mut done: [bool; NUM_BATCHES] = [false; NUM_BATCHES];
        let mut batch_masks: [u128; NUM_BATCHES] = [0; NUM_BATCHES];
        for j in 0..NUM_PRIMES {
            batch_masks[j % NUM_BATCHES] |= 1 << j;
        }
        let mut cur_batch = 0;
        let mut i = 0;
        loop {
            assert!(cur_batch < NUM_BATCHES);
            let mut early_finish = 0;
            if i > MERGE_AFTER * NUM_BATCHES {
                cur_batch = 0;
                batch_masks[cur_batch] = 0;
                for i in 0..NUM_PRIMES {
                    if isogeny_count[i] != 0 {
                        batch_masks[cur_batch] |= 1 << i;
                        done[cur_batch] = false;
                    }
                }
                if done[cur_batch] {
                    return e;
                }
            } else {
                while done[cur_batch] {
                    if early_finish == NUM_BATCHES {
                        return e;
                    }
                    early_finish += 1;
                    cur_batch = (cur_batch + 1) % NUM_BATCHES;
                }
            }
            let (p_0, p_1) = e.elligator(&mut thread_rng());
            let mut p_0 = e.double(&e.variable_time_differential_addition_chain(
                &e.double(&p_0),
                &(!batch_masks[cur_batch]),
            ));
            let mut p_1 = e.double(&e.variable_time_differential_addition_chain(
                &e.double(&p_1),
                &(!batch_masks[cur_batch]),
            ));
            for prime_index in (0..NUM_PRIMES).rev() {
                if (batch_masks[cur_batch] >> prime_index) & 1 == 0 {
                    continue;
                }
                let cur_exp = blinding[prime_index];
                let sign_bit = ((cur_exp >> 7) & 1) as u64;
                let (mut p_s, mut p_1s) = (p_0, p_1);
                p_s.conditional_move(sign_bit, &p_1);
                p_1s.conditional_move(sign_bit, &p_0);

                let mask: u128 = 1 << prime_index;
                batch_masks[cur_batch] &= !mask;
                let k = e.variable_time_differential_addition_chain(&p_s, &batch_masks[cur_batch]);
                p_1s = e.variable_time_differential_addition_chain(&p_1s, &mask);

                //random information
                if !k.is_zero() {
                    let uexp = cur_exp as u8;
                    let is_non_zero = (uexp | !uexp.wrapping_sub(1)) >> 7;
                    let (p_s_tmp, p_1s_tmp, e_tmp) =
                        e.two_point_isogeny(&k, PRIMES16[prime_index] as usize, &p_s, &p_1s);
                    p_s = e.variable_time_differential_addition_chain(&p_s, &mask);
                    p_s.conditional_move(is_non_zero as u64, &p_s_tmp);
                    p_1s.conditional_move(is_non_zero as u64, &p_1s_tmp);
                    e.a.conditional_move(is_non_zero as u64, &e_tmp.a);
                    let update = 1 - (2 * sign_bit as i8);
                    blinding[prime_index] -= update * is_non_zero as i8;
                    isogeny_count[prime_index] -= 1;
                }
                (p_0, p_1) = (p_s, p_1s);
                p_0.conditional_move(sign_bit, &p_1s);
                p_1.conditional_move(sign_bit, &p_s);
            }
            assert_eq!(batch_masks[cur_batch], 0);
            for i in 0..NUM_PRIMES {
                if i % NUM_BATCHES == cur_batch {
                    batch_masks[cur_batch] |= ((isogeny_count[i] != 0) as u128) << i;
                }
            }
            done[cur_batch] = batch_masks[cur_batch] == 0;
            cur_batch = (cur_batch + 1) % NUM_BATCHES;
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::csifish::field_arithmetic::classgroup::ClassGroupElement;
    use super::*;
    use crate::csifish::constants::BASE_CURVE;

    #[test]
    fn ct_sdac() {
        for _ in 0..1 {
            let j = ClassGroupElement::random(&mut thread_rng()).reduce_one_round();
            let e1 = j.one_time_blinded_action(&BASE_CURVE);
            println!("{}", e1.a);
        }
        // let j = ClassGroupElement::random().reduce();
        // let e1 = j.one_time_blinded_action(&BASE_CURVE);
        // let e2 = j.variable_time_action(&BASE_CURVE);
        // println!("{}", e1.a);
        // println!("{}", e2.a);
    }

    #[test]
    fn constant_time_isogeny() {
        let a = FieldElement::ZERO;
        let b = FieldElement::from_be_hex("53BAA451F759835A01933C76BC58C0C203A9B6B02F7F086B30C3469A8452750AAECA8A4F7C26BFF43876F4510F405F4D2A006635D89A42D327D9A2E8C00BF340");
        let order: usize = 3;

        let k = Point::from_x(FieldElement::from_be_hex("22B668C942BF7D5F5DF869A215F7E9463A0A873CFE2953721F129EC98B8123A8E62DF0D1F100AA92F4C6C8552AD62C42C11DB1AE8540F46ADC16D8939808553A"));
        let p = Point::from_x(FieldElement::from_be_hex("0A3A72458C434F22FD1F2B441C3BAD38C0C069872F69372A43E818126CFF49DC3CA63E87BC5F0443201F9DA03EFE8DA618C4D207954D40F774A923CBC11F2CA7"));
        let im_p = Point::from_x(FieldElement::from_be_hex("1ED168610F98DC95AAB55E2B067E92B32AF0A436A73EF7142F31BC3CBE2A532F8D51061DA110C5EB01FEC1838C6D0AA3B643D90181AAA3184CF02ABB20ECFB2A"));

        let ea = MontgomeryCurve::new(a);

        let (xp, _, xb) = ea.two_point_isogeny(&k, order, &p, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = FieldElement::from_be_hex("48211766D23E629D22C38ED44B3D8A02622B7022E5CE2CE5CCF7CDD4F901213AE61B00371E74AD24C9F71C59C0B0269287B36EC9652F4ACC421B8975C8C9EE4F");
        let b = FieldElement::from_be_hex("20B68C844B20BBA2271497C8ECD471D2EB0E3640A3D238F142C13C3C86BDF9D2F758186586740B2A15F9709E18F93F7894704B23CCBB533AC8AD2F1031AE309B");
        let order: usize = 5;
        let k = Point::from_x(FieldElement::from_be_hex("409548FAF7B5117391A5AD4D1202CA9EE096D69F44188441796F2ACED23C0C21DA29C9286AD5A46636CE1E41F9F54CEF4F453F7EFFCD595E168CC519DD68EA51"));
        let p = Point::from_x(FieldElement::from_be_hex("3B2DC0FD5EE8C65F43DAD597D8C48C32138A9FE4A1008802D5CED33523731EB432469E2D7F2276625E3DF38566576180E559E1C13D5F9565696A6D0D83830FF4"));
        let im_p = Point::from_x(FieldElement::from_be_hex("3BF01DE995EB675B0C2303367BC0FC3F82AC3D7123F842DEC8DE1E34F6FE14FBFDC1BDF203914BF7F6C52A7AF66CA745D4C682A4D1C9F40D1D5CFB066BB46B2D"));
        let ea = MontgomeryCurve::new(a);

        let (xp, _, xb) = ea.two_point_isogeny(&k, order, &p, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = FieldElement::from_be_hex("48912381B13014DF7E10F242424DFE6D43860ED48A2913843A45E75E15615849B2E2C8191E6CEF70A931E20883E8B59B87046926B8E534DCA88722A8E204496C");
        let b = FieldElement::from_be_hex("4E7D723C463A2F779721CD1F53CB1F2F3F9ECEB60E1831A2DDC665687C1F7BD1B479670592F4967DFAD3F9675B6229ED2B4ECFC2AAE33258DAF6A1B5CCBF2B78");
        let order: usize = 173;
        let k = Point::from_x(FieldElement::from_be_hex("54C8CDA4F5B40B4DD5EF9011AFA313195A68106114B157B53270CB1005C8338E4CE00C826ECEE406027F383FA1D5037DBB81D92E4203B4092B9C3D20A32D49A8"));
        let p = Point::from_x(FieldElement::from_be_hex("55E59A6BB770F1477F38E747D4C45F61CDA4D068736398DCB7C3A6B872208E6BA55FA42377A4B3EB25AEF4D0CE59C91A1D3A291B87700FFFE21805A7DEED199B"));
        let im_p = Point::from_x(FieldElement::from_be_hex("5BBC4E76109DAF1BCC7A597C78DCA56C1645CA6C72859B3F316F972054BF200C0F8059E2CCD9B1886F7518230CC5A75E210A3A5C07D843FF79BD832B675E1BD3"));
        let ea = MontgomeryCurve::new(a);

        let (xp, _, xb) = ea.two_point_isogeny(&k, order, &p, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = FieldElement::from_be_hex("4F80CF43DC32028D21AF9A4596C7067352C764156B62056D1DBC2A528E367DAC0BC65E453F01864DF53E0A775E064988EDD71034EBE1D5B95C1235F11DC6ACFD");
        let b = FieldElement::from_be_hex("1AFAA394A786BAEB11895EE8A455AE6A6872C74C9D2F0F47773AAB2FD1481BACC7695E7B81177C643C054D3BD36268F5ABD7AE225EEBEFA531F9153F532E4BA5");
        let order: usize = 83;
        let k = Point::from_x(FieldElement::from_be_hex("425A0D7407BF49078B071367E138506CDF3CF5C5231384524F9C62C7E84BF1536C47B5AC7B981BD1D8B8A4FF5ED0A75471F0D80ED1515CE18C2D31780929FD58"));
        let p = Point::from_x(FieldElement::from_be_hex("1C33BAADF7E34ACA1AFE98CEF02BA3948B0A09A2996BC9BC2C28A4E33A4943E8FA63370DF59C0CDA3C1D943473E50B2D4334DEE8263F6CF6450619460BF09DC4"));
        let im_p = Point::from_x(FieldElement::from_be_hex("46C02021FCC07A06A4C41958C3C40EB31DF7D10947C1021FE2638A9DADB7C8792D8EC0271FE63DFDE0BF6E1B4D44E550A9606DC91541DC15263292469892BD8D"));
        let ea = MontgomeryCurve::new(a);

        let (xp, _, xb) = ea.two_point_isogeny(&k, order, &p, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);

        let a = FieldElement::from_be_hex("38C29C0734D6DD3278A7DCDB797A4D8D1E6C41ADAC0FA768EE3BF9EE8FBC3EECE077F09DE1644631226B822B3BA01868DE4A9E7603DADA2D024BB5B16E083020");
        let b = FieldElement::from_be_hex("2FDC10B37A9570DFA25AEE9482802ECDC60E7D5D47B6E06AC5C3114BA70DDF38C6D820DFAD5A126794DC0CCF78A3BB91283DEAE9D6B540EB45506934D5C145B7");
        let order: usize = 149;
        let k = Point::from_x(FieldElement::from_be_hex("5E9D1A6638D9610AE0568BC36A483E512E3AC582C45E79A5388D0C213F7315052B4B74784F468E9D5CE9EA882D9511AF4A7B92E7CDEC4D5AE22D32D8B9F805DE"));
        let p = Point::from_x(FieldElement::from_be_hex("35EE2441FB15C6C330837FA2950A9860C33A14E6847D78DE6EB62FF85291477CEB7E69CE825B88637283A87379AC17D3A1E319A2D95172CFFC31FE6380C54749"));
        let im_p = Point::from_x(FieldElement::from_be_hex("3A12C8F02C85892393291F5860DB7B8C86C198FE89B44B165A91E9C05F185896C036B64331A418347706C6D124B73AECE248925112F207E3E53114FEECE14545"));
        let ea = MontgomeryCurve::new(a);

        let (xp, _, xb) = ea.two_point_isogeny(&k, order, &p, &p);
        let xb = MontgomeryCurve::new(xb.a.normalize().x);
        assert_eq!(xb.a.x, b);
        assert_eq!(im_p.normalize().x, xp.normalize().x);
    }

    #[test]
    fn constant_time_action() {
        let e = MontgomeryCurve::new(FieldElement::ZERO);
        let exp: [i8; NUM_PRIMES] = [
            -5, 2, 0, -3, 4, -4, -5, 3, 5, -1, -2, -4, 0, -2, -3, 3, 1, -2, 5, 3, 4, 3, -4, 2, 2,
            3, -1, 0, 1, -3, 0, 1, -5, -2, 0, 2, 0, 0, -5, 5, 4, 5, 0, -5, 0, -1, 0, 1, 5, 1, 1,
            -3, 0, 5, 1, 2, -1, 1, -5, 0, 1, 5, 3, 2, -1, -5, 4, 2, 1, 2, -2, 0, 1, 5,
        ];

        let g1 = ReducedClassGroupElement::new(exp);
        let b = g1.clone().one_time_blinded_action(&e);

        let correct1 = FieldElement::from_be_hex("2D3F42F31F984ACE1F45E62D35F7C9936BA51863A204A7AF9562DF7822E01323EAECAB2D86BBA42CB9B1DAA7DAA565800BD5BF35A0297218E8CBDB0399618180");

        assert_eq!(b.normalize().a.x, correct1);

        let exp: [i8; NUM_PRIMES] = [
            1, -2, 5, 1, 2, 4, -1, 0, -2, -1, 2, 5, -3, 3, 3, -1, -2, -1, 0, -5, -1, -1, -5, 4, 2,
            -1, -1, -5, -4, -3, 4, 1, 4, -2, 4, -5, 3, -1, 1, 2, 0, 4, 1, -5, 4, 1, 4, -1, 0, -5,
            3, -2, -3, 0, -1, 4, 3, -2, -5, -5, 4, 3, 2, 1, -2, 3, 3, -2, -3, -5, 5, 3, -5, 2,
        ];

        let g2 = ReducedClassGroupElement::new(exp);
        let b = g2.clone().one_time_blinded_action(&e);

        let correct2 = FieldElement::from_be_hex("09EB001955B4E84ECFFE86806E0C8313800D0475CFF3519FAF30DC5F3A060E97AE258051DABED0245406DF3BD41B4A03F3C7756C2DE8DE4AD28AC8CD8D506695");

        assert_eq!(b.normalize().a.x, correct2);

        let e1 = MontgomeryCurve::new(correct1);
        let c = g2.clone().one_time_blinded_action(&e1);

        let correct3 = FieldElement::from_be_hex("2BA3EBCD76B29349F525D3B73BA841065926870C3A1F23902EF53652D880BCF6E8D2705B2F94E23551BBFE9F4FD9A4DA1EADF24EA62DC2A7F425A8EB901E31A6");

        assert_eq!(c.normalize().a.x, correct3);

        let e2 = MontgomeryCurve::new(correct2);
        let c = g1.clone().one_time_blinded_action(&e2);
        assert_eq!(c.normalize().a.x, correct3);

        let a = FieldElement::from_be_hex("5EB2AEEF49060ED93CC067CC83EDDA45D2494F1CF0EB19F41DA034D00A61CBFDFE7C05C0E2730E14EE51B1C0DD5F10CD4958FB9567E9125410860FADDE6D5306");
        let e3 = MontgomeryCurve::new(a);
        let exp: [i8; NUM_PRIMES] = [
            -5, 2, -5, -1, -4, -3, 5, 4, -2, 5, 3, -4, 4, 4, 5, 5, -5, -1, -2, -1, 2, -3, 1, -5,
            -2, 5, 5, 5, -2, 2, 3, 4, 2, -5, 4, 2, 1, 4, -3, 1, -3, 0, 5, 4, -4, 0, 0, -3, 1, 3, 0,
            -1, -4, -4, -5, -4, -5, 3, -3, 0, -4, 2, -1, 4, 5, 0, -3, 3, -4, -1, -2, -2, 2, -5,
        ];

        let g = ReducedClassGroupElement::new(exp);
        let b = g.clone().one_time_blinded_action(&e3);

        let correct = FieldElement::from_be_hex("244BFECEC58AB059E806D5E001BFA1230F5FD3735C2D78EA8F901E4C3FDE881D2FBE39781C948436C538EFA0E2C54B650E390D2B519BD6A3BA6026AFCC819DCB");
        assert_eq!(correct, b.normalize().a.x);

        let b = g.clone().one_time_blinded_action(&e);
        let correct = FieldElement::from_be_hex("0313B6847C6679D3E73A9DD53E2C48E7E1279BE4749D519B2CC13FF5F7D8B235944A1994761C0DFD8306A899567D1DE98ECE0F2431C907EAC61CD5E1F34E0E9E");
        assert_eq!(correct, b.normalize().a.x);
    }
}
