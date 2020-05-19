/*
MNT4 group gadgets G1Gadget and G2Gadget as well as G1PreparedGadget and G2PreparedGadget for
the Miller loop, including their serialization (toBytesGadget)
    Maybe better to move all the preperatory gadgets to pairing-related code
*/

use algebra::Field;

use crate::{fields::{
    FieldGadget, fp::FpGadget, fp2::Fp2Gadget,
}, groups::curves::short_weierstrass::short_weierstrass_projective::AffineGadget,
    bits::uint8::UInt8, Assignment,
    alloc::AllocGadget,
    ToBytesGadget,
};

use r1cs_core::{ConstraintSystem, SynthesisError};
use algebra::curves::models::mnt4::MNT4Parameters;
use std::fmt::Debug;
use std::ops::Mul;
use crate::bits::boolean::Boolean;

pub mod mnt4753;

pub type G1Gadget<P> = AffineGadget<
    <P as MNT4Parameters>::G1Parameters,
    <P as MNT4Parameters>::Fp,
    FpG<P>
>;
pub type G2Gadget<P> = AffineGadget<
    <P as MNT4Parameters>::G2Parameters,
    <P as MNT4Parameters>::Fp,
    Fp2G<P>
>;

type FpG<P> = FpGadget<<P as MNT4Parameters>::Fp>;
type Fp2G<P> = Fp2Gadget<<P as MNT4Parameters>::Fp2Params, <P as MNT4Parameters>::Fp>;

#[derive(Derivative)]
#[derivative(
Clone(bound = "FpGadget<P::Fp>: Clone"),
Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
Debug(bound = "FpGadget<P::Fp>: Debug"),
Debug(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Debug"),
)]
pub struct G1PreparedGadget<P: MNT4Parameters> {
    pub p:                   G1Gadget<P>,
    pub p_y_twist_squared:   Fp2G<P>,
}

impl<P: MNT4Parameters> G1PreparedGadget<P> {
    pub fn from_affine<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        value: &G1Gadget<P>,
    ) -> Result<Self, SynthesisError> {

        let p = value.clone();
        let twist_squared = P::TWIST.square();
        let c0 = p.y.mul_by_constant(cs.ns(||"p.y * twist_squared.c0"), &twist_squared.c0)?;
        let c1 = p.y.mul_by_constant(cs.ns(||"p.y * twist_squared.c1"), &twist_squared.c1)?;
        let p_y_twist_squared = Fp2G::<P>::new(c0, c1);
        Ok(G1PreparedGadget{p, p_y_twist_squared})
    }
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for G1PreparedGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut p = self.p.to_bytes(&mut cs.ns(|| "p to bytes"))?;
        p.extend_from_slice(&self.p_y_twist_squared.to_bytes(&mut cs.ns(|| "p_y_twist_squared to bytes"))?);
        Ok(p)
    }

    fn to_bytes_strict<CS: ConstraintSystem<P::Fp>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut p = self.p.to_bytes_strict(&mut cs.ns(|| "p to bytes"))?;
        p.extend_from_slice(&self.p_y_twist_squared.to_bytes_strict(&mut cs.ns(|| "p_y_twist_squared to bytes"))?);
        Ok(p)
    }
}

#[derive(Derivative)]
#[derivative(
Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
Debug(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Debug")
)]
pub struct G2CoefficientsGadget<P: MNT4Parameters> {
    pub(crate) r_y:            Fp2G<P>,
    pub(crate) gamma:          Fp2G<P>,
    pub(crate) gamma_x:        Fp2G<P>,
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for G2CoefficientsGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.r_y.to_bytes(&mut cs.ns(|| "r_y to bytes"))?;
        x.extend_from_slice(&self.gamma.to_bytes(&mut cs.ns(|| "gamma to bytes"))?);
        x.extend_from_slice(&self.gamma_x.to_bytes(&mut cs.ns(|| "gamma_x to bytes"))?);
        Ok(x)
    }

    fn to_bytes_strict<CS: ConstraintSystem<P::Fp>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.r_y.to_bytes_strict(&mut cs.ns(|| "r_y to bytes"))?;
        x.extend_from_slice(&self.gamma.to_bytes_strict(&mut cs.ns(|| "gamma to bytes"))?);
        x.extend_from_slice(&self.gamma_x.to_bytes_strict(&mut cs.ns(|| "gamma_x to bytes"))?);
        Ok(x)
    }
}

#[derive(Derivative)]
#[derivative(
Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
Debug(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Debug")
)]
pub struct G2PreparedGadget<P: MNT4Parameters>{
    pub q:      G2Gadget<P>,
    pub coeffs: Vec<G2CoefficientsGadget<P>>
}

impl<P: MNT4Parameters>G2PreparedGadget<P> {
    /* Takes as input a (non-zero) G2Gadget Q, and outputs the
     (P-independent) pre-computable coefficients for all operations of the Miller loop.
     These are exactly the same as in the pairing primitive
         s.y = the y-coordinate of internal state S,
         gamma = the F2-slope of the tangent/P-chord at S,
         gamma_x = the F2-slope times the x-coordinate s.x of S.
    */
    pub fn from_affine<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        value: &G2Gadget<P>,
    ) -> Result<Self, SynthesisError> {

        let mut s = value.clone();

        let mut g2p = G2PreparedGadget{
            q: s.clone(),
            coeffs: vec![]
        };

        for (i, &n) in P::WNAF.iter().rev().enumerate(){

            let mut cs = cs.ns(|| format!("Iteration {}", i));

            let (s2, c) = Self::doubling_step_for_flipped_miller_loop(cs.ns(|| "double"), &s.clone())?;
            g2p.coeffs.push(c);
            s = s2;
            if n != 0 {
                let (s2, c) = Self::mixed_addition_step_for_flipped_miller_loop(cs.ns(|| "add"), &value.x, &value.y, &s.clone(), n)?;
                g2p.coeffs.push(c);
                s = s2;
            }
        }
        Ok(g2p)
    }

    /* computes the preparation coefficients at S and its subsequent doubling
    */
    fn doubling_step_for_flipped_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        s: &G2Gadget<P>,
    ) -> Result<(G2Gadget<P>, G2CoefficientsGadget<P>), SynthesisError>
    {
        //Compute gamma, i.e. the F2-slope of the tangent at S
        let s_x_squared = s.x.square(cs.ns(||"s_x^2"))?;
        let three_sx_squared_plus_a = s_x_squared
            .double(cs.ns(|| "2s_x^2"))?
            .add(cs.ns(|| "3s_x^2"), &s_x_squared)?
            .add_constant(cs.ns(|| "3s_x^2 + a"), &P::TWIST_COEFF_A)?;

        let two_sy = s.y.double(cs.ns(||"2s_y"))?;

        let gamma = Fp2G::<P>::alloc(cs.ns(|| "gamma"), || {
            Ok(three_sx_squared_plus_a.get_value().get()?.mul(&two_sy.get_value().get()?.inverse().get()?))
        })?;

        //Check gamma (gamma*2s_y = sx^2 + 3a)
        gamma.mul_equals(cs.ns(|| "Check gamma"), &two_sy, &three_sx_squared_plus_a)?;

        //Compute and check gamma_x
        let gamma_x = gamma.mul(cs.ns(|| "Compute gamma_x"), &s.x)?;


        /* as we already computed the slope of the tangent, we re-implement the Weierstrass
        doubling formulas to save constraints.
        */
        //Compute and check new_sx
        let two_sx = s.x.double(cs.ns(|| "2s_x"))?;
        let new_sx = gamma.square(cs.ns(|| "gamma^2"))?
            .sub(cs.ns(|| "gamma^2 - 2s_x"), &two_sx)?;

        //Compute and check new_sy
        let new_sy = s.x
            .sub(cs.ns(|| "s_x - new_s_x"), &new_sx)?
            .mul(cs.ns(|| "gamma * (s_x - new_s_x)"), &gamma)?
            .sub(cs.ns(|| "gamma * (s_x - new_s_x) - s_y"), &s.y)?;

        let c = G2CoefficientsGadget{r_y: s.y.clone(), gamma, gamma_x};
        let s2 = G2Gadget::<P>::new(new_sx, new_sy, Boolean::constant(false));

        Ok((s2, c))
    }

    /* computes the preparation coefficients at S and its subsequent state by adding/subtracting
      Q=(x,y) depending on the naf "bit"
    */
    fn mixed_addition_step_for_flipped_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        x: &Fp2G<P>,
        y: &Fp2G<P>,
        s: &G2Gadget<P>,
        naf_i: i32,
    ) -> Result<(G2Gadget<P>, G2CoefficientsGadget<P>), SynthesisError>
    {
        //Compute gamma, i.e. the F2-slope of the chord between Q and S
        let sx_minus_x = s.x
            .sub(cs.ns(|| "s_x - x"), &x)?;

        let sy_plus_y = s.y.add(cs.ns(||"(s_y + y)"), &y)?;
        let sy_minus_y = s.y.sub(cs.ns(|| "(s_y - y)"), &y)?;
        let numerator = if naf_i > 0 {sy_minus_y} else {sy_plus_y};

        let gamma = Fp2G::<P>::alloc(cs.ns(|| "Compute gamma"), || {
            let sx_minus_x_inv = sx_minus_x.get_value().get()?.inverse().get()?;
            Ok(numerator.get_value().get()?.mul(&sx_minus_x_inv))
        })?;

        //Check gamma
        gamma.mul_equals(cs.ns(||"Check gamma"), &sx_minus_x, &numerator)?;

        //Compute and check gamma_x
        let gamma_x = gamma.mul(cs.ns(|| "Compute gamma_x"), &x)?;

        /* as we already computed the slope gamma, we re-implement Weierstrass addition to save constraints
        */
        //Compute and check new_sx
        let new_sx = gamma.square(cs.ns(|| "gamma^2"))?
            .sub(cs.ns(|| "gamma^2 - s_x"), &s.x)?
            .sub(cs.ns(|| "gamma^2 - s_x - x"), &x)?;

        //Compute and check new_sy
        let new_sy = s.x
            .sub(cs.ns(|| "s_x - new_s_x"), &new_sx)?
            .mul(cs.ns(|| "gamma * (s_x - new_s_x)"), &gamma)?
            .sub(cs.ns(|| "gamma * (s_x - new_s_x) - s_y"), &s.y)?;

        let c = G2CoefficientsGadget{r_y: s.y.clone(), gamma, gamma_x};
        let s2 = G2Gadget::<P>::new(new_sx, new_sy, Boolean::constant(false));

        Ok((s2, c))
    }
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for G2PreparedGadget<P>
{
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError>
    {
        let mut x = self.q.to_bytes(&mut cs.ns(|| "q to bytes"))?;

        for (i, c) in self.coeffs.iter().enumerate() {
            x.extend_from_slice(&c.to_bytes(&mut cs.ns(|| format!("coefficients {} to bytes", i)))?);
        }

        Ok(x)
    }

    fn to_bytes_strict<CS: ConstraintSystem<P::Fp>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.q.to_bytes_strict(&mut cs.ns(|| "q to bytes"))?;

        for (i, c) in self.coeffs.iter().enumerate() {
            x.extend_from_slice(&c.to_bytes_strict(&mut cs.ns(|| format!("coefficients {} to bytes", i)))?);
        }

        Ok(x)
    }
}