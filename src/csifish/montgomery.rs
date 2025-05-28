use core::fmt;

use rand::thread_rng;

use crate::csifish::field_arithmetic::base_field::FieldElement;
use crate::csifish::field_arithmetic::arithmetic::{ModularArithmetic, MontgomeryArithmetic};
use crate::csifish::constants::DeserializationError;

#[derive(Clone, Debug, PartialEq, Copy)]
pub struct Point {
    pub x: FieldElement,
    pub z: FieldElement,
}

impl Point {
    pub const fn from_x(x: FieldElement) -> Point {
        Point {
            x,
            z: FieldElement::ONE,
        }
    }

    pub const fn zero() -> Point {
        Point {
            x: FieldElement::ONE,
            z: FieldElement::ZERO,
        }
    }

    pub fn random() -> Point {
        Point::from_x(FieldElement::random(&mut thread_rng()))
    }

    pub fn is_zero(&self) -> bool {
        self.z.is_zero()
    }

    #[inline]
    pub fn normalize(&self) -> Point {
        let res = self.z.inv();
        if bool::from(res.is_some()) {
            Point {
                x: self.x * res.unwrap(),
                z: FieldElement::ONE,
            }
        } else {
            Point {
                x: FieldElement::ONE,
                z: FieldElement::ZERO,
            }
        }
    }

    pub fn to_be_bytes(&self) -> Vec<u8> {
        self.x.get_be_bytes().to_vec()
    }

    pub fn from_be_bytes(b: &[u8]) -> Result<Self, DeserializationError> {
        let x = FieldElement::from_be_bytes(b)?;
        Ok(Point::from_x(x))
    }

    pub fn conditional_move(&mut self, c: u64, b: &Point) {
        self.x.conditional_move(c, &b.x);
        self.z.conditional_move(c, &b.z)
    }

    pub fn conditional_swap(&mut self, c: u64, b: &mut Point) {
        self.x.conditional_swap(c, &mut b.x);
        self.z.conditional_swap(c, &mut b.z);
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.normalize();
        write!(f, "{} : {}", p.x, p.z)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct MontgomeryCurve {
    pub a: Point, // projective coefficient, A = a.x/a.z
}

impl MontgomeryCurve {
    pub const fn new(a: FieldElement) -> MontgomeryCurve {
        MontgomeryCurve {
            a: Point::from_x(a),
        }
    }

    #[inline(always)]
    pub(crate) const fn projective(ax: FieldElement, az: FieldElement) -> MontgomeryCurve {
        MontgomeryCurve {
            a: Point { x: ax, z: az },
        }
    }


    pub fn twist(&self) -> MontgomeryCurve {
        MontgomeryCurve::new(self.a.x.neg())
    }

    pub fn normalize(&self) -> MontgomeryCurve {
        let a = self.a.normalize();
        MontgomeryCurve { a }
    }

    // Given points P, Q, and either P+Q or P-Q, computes P-Q or P+Q
    // P, Q must be distinct points, neither of which are the origin or infinity (zero)
    // Algorithm 1 in Costello & Smith
    #[inline(always)]
    pub fn differential_add(&self, p: &Point, q: &Point, pq: &Point) -> Point {
        let v1 = (p.x + p.z) * (q.x - q.z);
        let v2 = (p.x - p.z) * (q.x + q.z);
        let v3 = (v1 + v2).square();
        let v4 = (v1 - v2).square();
        Point {
            x: pq.z * v3,
            z: pq.x * v4,
        }
    }

    // Given P not equal to origin or infinity, computes P+P
    // Projective version of Algorithm 2 in Costello & Smith,
    // obtained by multiplying through by 4a.z
    #[inline(always)]
    pub fn double(&self, p: &Point) -> Point {
        // V1 in Costello & Smith
        let vplus = (p.x + p.z).square(); // (x + z)^2

        // V2 in Costello & Smith
        let vminus = (p.x - p.z).square(); // (x - z)^2

        let vdelta = vplus - vminus; // V1 in Costello & Smith
        let va = vminus * self.a.z;
        let va = va + va + va + va;

        let x = vplus * va; // vminus * vplus * denominator

        let vb = (self.a.z + self.a.z + self.a.x) * vdelta + va;

        Point { x, z: vb * vdelta }
    }

    pub fn j(&self) -> FieldElement {
        let two_fifty_six: FieldElement = FieldElement::from_u16(256);
        let three: FieldElement = FieldElement::from_u8(3);
        let four: FieldElement = FieldElement::from_u8(4);

        let mut num = self.a.x.square() - three;
        num = two_fifty_six * num.square() * num;
        //not zero by construction
        let den = (self.a.x.square() - four).inv().unwrap();
        num * den
    }

    pub fn to_be_bytes(&self) -> Vec<u8> {
        self.a.to_be_bytes()
    }

    pub fn from_be_bytes(b: &[u8]) -> Result<Self, DeserializationError> {
        let x = FieldElement::from_be_bytes(b)?;
        Ok(MontgomeryCurve::new(x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn on_curve() {
    //     let a = FieldElement::from_be_hex("57164DAD2DAA6B17538CC28E418D0B93540024EC5F7038951049142A2FA46F030D7E5247B792A894FF526D7126DCB9CDEF42704493B6F8109CC5B127FD6F4888");
    //     let e = MontgomeryCurve::new(a);
    //
    //     let p = Point::from_x(FieldElement::from_be_hex("39F65EDE480BF8C5E5ACA9CE8EAC7EC98B8C02E3768B83444D77F06961B50B9CDCD8A3644D244624766ABECED59F881BA06B04033B6FA7652396BA8798A16CEA"));
    //     assert!(e.on_curve(&p));
    //
    //     let p = Point::from_x(FieldElement::from_be_hex("4ACCB8B58B35F4D9787314CB062D264A5CD43EC672B48CEAD6FE63FD49A94CC36912F751EAE0D262F1584DA663F2A3A18506EF9F4B444D7F40AAC0A7E5838869"));
    //     assert!(e.on_curve(&p));
    //
    //     let p = Point::from_x(FieldElement::from_be_hex("5727FDAA0E1070CFF054606C531F4BE7B6D55B9B2CC2F343C19306C76BFA01D247BBDC9F05E5AFD1D08453C532F1E733CC3419C868167AF4F03AC860A90E258F"));
    //     assert!(!e.on_curve(&p));
    //
    //     let p = Point::from_x(FieldElement::from_be_hex("63576655EC3AA520BF7A2B022D5253C18E8676C03EF81FA05030B2A4509F2E4C11A37DDDB06606A6CA94DAD30DD0876E300FD9AEAEF2B075DD5CE6AD74255B13"));
    //     assert!(!e.on_curve(&p));
    // }

    #[test]
    fn add3() {
        let m = FieldElement::from_be_hex("54C8A0DADC3C7B204FDF48616AF757968326DC25866F018424FCD27D45C809CAC0D3F58553D6CB42704819843C67406977C51CD790BE78350FADAB6CB72AFA8D");
        let mut e = MontgomeryCurve::new(FieldElement::ZERO);
        e.a.x *= m;
        e.a.z *= m;
        let p = Point::from_x(FieldElement::from_be_hex("2C269236F8C147E177F0A3B5C95EA91A4C04DE1BBA1FFC84622CE367805C8551A7D628FFEA33900BA7F6F90F65AA0EAFEA189CEC225A1730B7F1E28FF2C253B4"));
        let q = Point::from_x(FieldElement::from_be_hex("258DD34D9A793CE8396185303C65012683D12744E672B2195B38D7A2DF26242024C3E8740CCAFD3071388C488266E2BBCBA9F5B3F2D4252CBCDB4EF69DD58EDB"));
        let pq_sum = Point::from_x(FieldElement::from_be_hex("55DA2197E667E35145692AB34A06D9C62A379533C8430B98295F3F0A6828682BC7046F3DF84760F6C22CBBEA55B67EF83C53AB1CF103DE0CF0774719D9038C46"));
        let pq_diff = Point::from_x(FieldElement::from_be_hex("5A54CDD887887905BB57A2FB2C5D108A361052CAF9B80E93812C69E5BF5EA6C9869F43748FB12937C5A91C999AAF4A5D0CBA6A0B9B67A8339F15DA786626F7F3"));
        assert_eq!(e.differential_add(&p, &q, &pq_sum).normalize(), pq_diff);
        assert_eq!(e.differential_add(&p, &q, &pq_diff).normalize(), pq_sum);

        let a = FieldElement::from_be_hex("47D112C8D0BBF39D1983F677BE0CD423445C8BACA91B516EB3350F1CB95FFB454F4B0C18CE2EA540CE7B0932B951B365511CDBB82458DCA4D0ABBA04DB00D84D");
        let m = FieldElement::from_be_hex("52C521933A01AD67352ABAEE2BB6FDB4025BA653A1B6C5C8B939B5647EF56A8111640A7717FEEB38967FE1F7B653384D1E4C5E4DAA76686F97CA4DC9E405C543");
        let mut e = MontgomeryCurve::new(a);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(FieldElement::from_be_hex("55D5F791250BF6F55EF6DCDE2F1249DB59DAEA4FDAB940A53FC796293B137D58FC0888CF466C523C21EED52CDA1554705E3ECE421FDB0653E052DE3009F60C29"));
        let q = Point::from_x(FieldElement::from_be_hex("38AC0C3673F3DE0E5C049EA68C97CE5261EEE5EC6BB02F568E0D3B70A40E139F22C64995A40B1B5F90D3AFF943557C0E14BB8227C477FC7E822C5787450198DB"));
        let pq_sum = Point::from_x(FieldElement::from_be_hex("5C58A818F02D09EF5DC1CAABEC133D07D3D9EFB98984B1DD9B7636D4A4FBC6B47FBF51FCDB71570CE00C723F4083CF4E70F1467077FC99B0AB746004ECC445BD"));
        let pq_diff = Point::from_x(FieldElement::from_be_hex("174B7233B392E09C38F7DB5DA5D352585FBE4D16CBD626404D894CA3407AB399400FD9FFE0FDF1C76ADB7F11955F88F45805EBDD727271F5465BEB39D96E6782"));

        assert_eq!(e.differential_add(&p, &q, &pq_sum).normalize(), pq_diff);
        assert_eq!(e.differential_add(&p, &q, &pq_diff).normalize(), pq_sum);
    }

    #[test]
    fn double() {
        let m = FieldElement::from_be_hex("09F4B54FF7BCF319672BEDF8D865F241C27D52BE5D70F8E4BD18806AF9E5BFF3FC85635574EC6D8513679E8F1BBE290672AC5AB4A0106D05C10DD74B758F2589");
        let mut e = MontgomeryCurve::new(FieldElement::ZERO);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(FieldElement::from_be_hex("24FECD3460B2F53C30EC39BE2E4CE435EAC8EE1E5E57CC2A022D9C885F941346D78DFA27AF0BE7CBFE913E461310BB5713A85D7AC02AAD5B8E802D2F7D045199"));
        let q = Point::from_x(FieldElement::from_be_hex("50809AE8A1F22F00D2D655E4E6E7580A1C868E59AE4BE4E98A7AF3FD3BB0A26C24724B6E0C98EFA529C914E7F5DFA89D1CBBCEF8E2DB859AB9538A7787B10EFB"));
        let two_p = Point::from_x(FieldElement::from_be_hex("08D21FCB978BE16F2A6A9A5753418098BE7CE928B417E4672B448F99029BF4F9D771C6AE5356889A0A7F0DD30AC467343CEDFD25BDE221C8B82B9D208CB3F2C1"));
        let two_q = Point::from_x(FieldElement::from_be_hex("02C862DCE53D5BE8C3FDDBC31152CAF717D52579033EF877E8DCAD1BF333C3B590D045023BD3B365D0FB0EDDAFB2D74871BA618EAB37C285B91CD83CF395F848"));

        assert_eq!(e.double(&p).normalize(), two_p);
        assert_eq!(e.double(&q).normalize(), two_q);

        let a = FieldElement::from_be_hex("20A24B2BE41B432B63B7C485E95F781E344744881BB3C9C36A10995DFB62592FFFDF6EAEA75F632894CC2D30F0C273AF4002D580E272FC2B4EB2C46730A6D7CD");
        let m = FieldElement::from_be_hex("3DB6207324441ECC3884C3865E972993CA0C239CAABDC1E4BFEB70E73CA776A83C8AF9E383B93CEDE3E6BA05503D422EEE7FAA5484F4711E69939FA3FDF7DD97");
        let mut e = MontgomeryCurve::new(a);
        e.a.x *= m;
        e.a.z *= m;

        let p = Point::from_x(FieldElement::from_be_hex("228B3C1B541F1C50E6E46545E9E3E26F747EC6D18D436CD02DA23B8AF7A37FD81F1132833AC5F96EBD3AF861F2CD83EAF429017B859EB55D6A7BAF514849B1BB"));
        let q = Point::from_x(FieldElement::from_be_hex("3F53508C3824E5B18630C63E35FD312E39B36345E7BA855D1AE5DA044682327A8FBE7627C46E707CD4E5F87BD0EB191D366DC0A7D7F2ECDFC7E6F9ECFC6A55D2"));
        let two_p = Point::from_x(FieldElement::from_be_hex("29306C62A1CFDDDEFE22DB90F5C3BF81DDA85BF6F25B45F37D6ECA62D633066CE47D707B9698EED38F1C2FD9D60E13EF33B712603FEE555A2EEFF68177B50A47"));
        let two_q = Point::from_x(FieldElement::from_be_hex("0E417B11BA938C20C957E22F787757DD7B138C633DF90F529DD0EC490A7FA49E847644AE5728E981DDCFC24D37B93EF169502C66FDBDE16D6D0A986DB0640738"));

        assert_eq!(e.double(&p).normalize(), two_p);
        assert_eq!(e.double(&q).normalize(), two_q);
    }
}
