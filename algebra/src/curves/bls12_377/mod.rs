//! [Zexe's BLS12-377](https://eprint.iacr.org/2018/962.pdf), the outer curve of edwards_bls12 and
//! inner curve of the sw6.
//! Its a BLS curve with embedding degree 12 and base field size is 377 bits.
//!
//! Its security against modern version of the number field sieve (STNFS) is estimated as 125 bit
//! by [HG 2020](https://eprint.iacr.org/2020/351.pdf).

// The HG 2020 result states an unexpected negligible loss in security compared to the recommendations
// from [Guillevic 2019](https://eprint.iacr.org/2019/1371) which establishes a base field size
// of 446 Bit to accomodate that security level.
// I hope that HG will give us some hints on that.

use crate::{
    curves::{
        bls12::{
            Bls12, Bls12Parameters, G1Affine as Bls12G1Affine, G1Prepared,
            G1Projective as Bls12G1Projective, G2Affine as Bls12G2Affine, G2Prepared,
            G2Projective as Bls12G2Projective, TwistType,
        },
        PairingCurve, PairingEngine,
    },
    fields::bls12_377::{Fq, Fq12, Fq12Parameters, Fq2Parameters, Fq6Parameters},
};

pub mod g1;
pub mod g2;
#[cfg(test)]
mod tests;

use self::{g1::Bls12_377G1Parameters, g2::Bls12_377G2Parameters};

pub struct Bls12_377Parameters;

impl Bls12Parameters for Bls12_377Parameters {
    /// BLS parameter x = 9586122913090633729
    const X: &'static [u64] = &[0x8508c00000000001];
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::D;
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = Bls12_377G1Parameters;
    type G2Parameters = Bls12_377G2Parameters;
}

pub type Bls12_377 = Bls12<Bls12_377Parameters>;

pub type G1Affine = Bls12G1Affine<Bls12_377Parameters>;
pub type G1Projective = Bls12G1Projective<Bls12_377Parameters>;
pub type G2Affine = Bls12G2Affine<Bls12_377Parameters>;
pub type G2Projective = Bls12G2Projective<Bls12_377Parameters>;

impl PairingCurve for G1Affine {
    type Engine = Bls12_377;
    type Prepared = G1Prepared<Bls12_377Parameters>;
    type PairWith = G2Affine;
    type PairingResult = Fq12;

    fn prepare(&self) -> Self::Prepared {
        Self::Prepared::from_affine(*self)
    }

    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult {
        Bls12_377::pairing(*self, *other)
    }
}

impl PairingCurve for G2Affine {
    type Engine = Bls12_377;
    type Prepared = G2Prepared<Bls12_377Parameters>;
    type PairWith = G1Affine;
    type PairingResult = Fq12;

    fn prepare(&self) -> Self::Prepared {
        Self::Prepared::from_affine(*self)
    }

    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult {
        Bls12_377::pairing(*other, *self)
    }
}
