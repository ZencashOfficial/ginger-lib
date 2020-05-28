//! Model for BLS12 curves and their Ate pairing.
//!
//! Supports the following standard efficiency measures for pairings:
//!
//! - G2 is represented by a sextic twist over a quadratic extension of the base field,
//! - the Frobenius operator is applied to reduce the cost of the final exponentiation, and
//! - we do pre-computations of (essentially) the line coefficients of the Miller loop.

use crate::{
    curves::{
        models::{ModelParameters, SWModelParameters},
        PairingCurve, PairingEngine,
    },
    fields::{
        fp12_2over3over2::{Fp12, Fp12Parameters},
        fp2::Fp2Parameters,
        fp6_3over2::Fp6Parameters,
        BitIterator, Field, Fp2, PrimeField, SquareRootField,
    },
};

use std::marker::PhantomData;

pub enum TwistType {
    M, //  E': y^2 = x^3 + b/i^5 which is isomorphic to y^2=x^3 + i*b
    D, // E': y^2 = x^3 + b/i
}

/// Parameters of an BLS12 curve E: y^2 = x^3 + b mod p as needed by the Ate pairing.
///
/// As in [ABLR 2013](https://eprint.iacr.org/2013/722.pdf), the embedding field F12 is regarded
/// as towered extension
///
///     F12 = F6[Z]/(Z^2-Y),
///     F6 = F2[Y]/(Y^3-U),
///     F2 = Fp[X]/(X^2-alpha),
///
/// using a non-square alpha from Fp, and a non-cube U from F2, respectively.
/// We apply standard efficiency measures: G2 is represented by a subgroup
/// of prime order r=ord(G1) of the sextic twist over F2,
///
///     E6: y^2 = x^3 + b/i,
/// or
///
///     E6: y^2 = x^3 + b/i^5,
///
/// depending on the twist type D and M, respectively.
/// Here, the twist element i = U, the non-residue for building F6 on top of F2.
pub trait Bls12Parameters: 'static {
    // Ate loop count, equals abs(Frobenius trace - 1)
    const X: &'static [u64];
    // sign of the Ate loop count
    const X_IS_NEGATIVE: bool;
    // type of the sextic twist
    const TWIST_TYPE: TwistType;
    // Base field
    type Fp: PrimeField + SquareRootField + Into<<Self::Fp as PrimeField>::BigInt>;
    // extension field for the sextic twist
    type Fp2Params: Fp2Parameters<Fp = Self::Fp>;
    // intermediate field as cubic extension of Fp2
    type Fp6Params: Fp6Parameters<Fp2Params = Self::Fp2Params>;
    // embedding field
    type Fp12Params: Fp12Parameters<Fp6Params = Self::Fp6Params>;
    type G1Parameters: SWModelParameters<BaseField = Self::Fp>;
    type G2Parameters: SWModelParameters<
        BaseField = Fp2<Self::Fp2Params>,
        ScalarField = <Self::G1Parameters as ModelParameters>::ScalarField,
    >;
}

pub mod g1;
pub mod g2;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Bls12<P: Bls12Parameters>(PhantomData<fn() -> P>);

impl<P: Bls12Parameters> Bls12<P> {
    // Evaluate the line function at point p.
    fn ell(
        f: &mut Fp12<P::Fp12Params>,
        coeffs: &(Fp2<P::Fp2Params>, Fp2<P::Fp2Params>, Fp2<P::Fp2Params>),
        p: &G1Affine<P>,
    ) {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;
        let mut c2 = coeffs.2;

        match P::TWIST_TYPE {
            TwistType::M => {
                c2.mul_by_fp(&p.y);
                c1.mul_by_fp(&p.x);
                f.mul_by_014(&c0, &c1, &c2);
            },
            TwistType::D => {
                c0.mul_by_fp(&p.y);
                c1.mul_by_fp(&p.x);
                f.mul_by_034(&c0, &c1, &c2);
            },
        }
    }

    fn exp_by_x(mut f: Fp12<P::Fp12Params>) -> Fp12<P::Fp12Params> {
        f = f.cyclotomic_exp(P::X);
        if P::X_IS_NEGATIVE {
            f.conjugate();
        }
        f
    }
}

impl<P: Bls12Parameters> PairingEngine for Bls12<P>
where
    G1Affine<P>: PairingCurve<
        BaseField = <P::G1Parameters as ModelParameters>::BaseField,
        ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
        Projective = G1Projective<P>,
        PairWith = G2Affine<P>,
        Prepared = G1Prepared<P>,
        PairingResult = Fp12<P::Fp12Params>,
    >,
    G2Affine<P>: PairingCurve<
        BaseField = <P::G2Parameters as ModelParameters>::BaseField,
        ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
        Projective = G2Projective<P>,
        PairWith = G1Affine<P>,
        Prepared = G2Prepared<P>,
        PairingResult = Fp12<P::Fp12Params>,
    >,
{
    type Fr = <P::G1Parameters as ModelParameters>::ScalarField;
    type G1Projective = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G2Projective = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type Fq = P::Fp;
    type Fqe = Fp2<P::Fp2Params>;
    type Fqk = Fp12<P::Fp12Params>;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<
            Item = &'a (
                &'a <Self::G1Affine as PairingCurve>::Prepared,
                &'a <Self::G2Affine as PairingCurve>::Prepared,
            ),
        >,
    {
        let mut pairs = vec![];
        for &(p, q) in i {
            if !p.is_zero() && !q.is_zero() {
                pairs.push((p, q.ell_coeffs.iter()));
            }
        }

        let mut f = Self::Fqk::one();

        for i in BitIterator::new(P::X).skip(1) {
            f.square_in_place();

            for &mut (p, ref mut coeffs) in &mut pairs {
                Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
            }

            if i {
                for &mut (p, ref mut coeffs) in &mut pairs {
                    Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
                }
            }
        }

        if P::X_IS_NEGATIVE {
            f.conjugate();
        }

        f
    }

    fn final_exponentiation(f: &Self::Fqk) -> Option<Self::Fqk> {
        // Computing the final exponentation following
        // https://eprint.iacr.org/2016/130.pdf.
        // We don't use their "faster" formula because it is difficult to make
        // it work for curves with odd `P::X`.
        // Hence we implement the algorithm from Table 1 below.

        // f1 = r.conjugate() = f^(p^6)
        let mut f1 = *f;
        f1.frobenius_map(6);

        match f.inverse() {
            Some(mut f2) => {
                // f2 = f^(-1);
                // r = f^(p^6 - 1)
                let mut r = f1 * &f2;

                // f2 = f^(p^6 - 1)
                f2 = r;
                // r = f^((p^6 - 1)(p^2))
                r.frobenius_map(2);

                // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
                // r = f^((p^6 - 1)(p^2 + 1))
                r *= &f2;

                // Hard part of the final exponentation is below:
                // From https://eprint.iacr.org/2016/130.pdf, Table 1
                let mut y0 = r.cyclotomic_square();
                y0.conjugate();

                let mut y5 = Self::exp_by_x(r);

                let mut y1 = y5.cyclotomic_square();
                let mut y3 = y0 * &y5;
                y0 = Self::exp_by_x(y3);
                let y2 = Self::exp_by_x(y0);
                let mut y4 = Self::exp_by_x(y2);
                y4 *= &y1;
                y1 = Self::exp_by_x(y4);
                y3.conjugate();
                y1 *= &y3;
                y1 *= &r;
                y3 = r;
                y3.conjugate();
                y0 *= &r;
                y0.frobenius_map(3);
                y4 *= &y3;
                y4.frobenius_map(1);
                y5 *= &y2;
                y5.frobenius_map(2);
                y5 *= &y0;
                y5 *= &y4;
                y5 *= &y1;
                Some(y5)
            },
            None => None,
        }
    }
}
