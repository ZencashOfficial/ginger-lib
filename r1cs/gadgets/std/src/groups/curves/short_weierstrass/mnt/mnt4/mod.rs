use algebra::{Field, AffineCurve};
use algebra::curves::{
    models::mnt4::{MNT4Parameters, g1::G1Prepared, g2::{G2PreparedCoefficients, G2Prepared}},
    ProjectiveCurve,
};

use crate::{fields::{
    FieldGadget, fp::FpGadget, fp2::Fp2Gadget,
}, groups::{
    curves::short_weierstrass::short_weierstrass_projective::AffineGadget,
    GroupGadget,
},
            bits::{
                boolean::Boolean,
                uint8::UInt8,
            }, Assignment,
            alloc::{AllocGadget, ConstantGadget},
            select::CondSelectGadget,
            ToBytesGadget,
};

use r1cs_core::{ConstraintSystem, SynthesisError};
use std::fmt::Debug;
use std::ops::Mul;

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
Clone(bound = "G1Gadget<P>: Clone"),
Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
Debug(bound = "G1Gadget<P>: Debug"),
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

impl<P: MNT4Parameters> ConstantGadget<G1Prepared<P>, P::Fp> for G1PreparedGadget<P> {

    #[inline]
    fn from_value<CS: ConstraintSystem<P::Fp>>(mut cs: CS, value: &G1Prepared<P>) -> Self {
        let p = G1Gadget::<P>::from_value(cs.ns(|| "hardcode p"), &value.p.into_projective());
        let p_y_twist_squared = Fp2G::<P>::from_value(cs.ns(|| "hardcode p_y_twist_squared"), &value.py_twist_squared);

        Self {
            p,
            p_y_twist_squared,
        }
    }

    fn get_constant(&self) -> G1Prepared<P> {
        G1Prepared::<P>{
            p: self.p.get_value().unwrap().into_affine(),
            py_twist_squared: self.p_y_twist_squared.get_value().unwrap(),
        }
    }
}


#[derive(Derivative)]
#[derivative(
Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
Debug(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Debug"),
Eq(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Eq"),
PartialEq(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: PartialEq"),
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

impl<P: MNT4Parameters> ConstantGadget<G2PreparedCoefficients<P>, P::Fp> for G2CoefficientsGadget<P> {

    #[inline]
    fn from_value<CS: ConstraintSystem<P::Fp>>(mut cs: CS, value: &G2PreparedCoefficients<P>) -> Self {
        let r_y = Fp2G::<P>::from_value(cs.ns(|| "hardcode r_y"), &value.r_y);
        let gamma = Fp2G::<P>::from_value(cs.ns(|| "hardcode gamma"), &value.gamma);
        let gamma_x = Fp2G::<P>::from_value(cs.ns(|| "hardcode gamma_x"), &value.gamma_x);

        Self {
            r_y,
            gamma,
            gamma_x,
        }
    }

    fn get_constant(&self) -> G2PreparedCoefficients<P> {
        G2PreparedCoefficients::<P>{
            r_y: self.r_y.get_value().unwrap(),
            gamma: self.gamma.get_value().unwrap(),
            gamma_x: self.gamma_x.get_value().unwrap(),
        }
    }
}

impl<P: MNT4Parameters> CondSelectGadget<P::Fp> for G2CoefficientsGadget<P> {
    fn conditionally_select<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self
    ) -> Result<Self, SynthesisError>
    {
        let r_y = Fp2G::<P>::conditionally_select(
            cs.ns(|| "r_y"),
            cond,
            &first.r_y,
            &second.r_y
        )?;

        let gamma = Fp2G::<P>::conditionally_select(
            cs.ns(|| "gamma"),
            cond,
            &first.gamma,
            &second.gamma
        )?;

        let gamma_x = Fp2G::<P>::conditionally_select(
            cs.ns(|| "gamma_x"),
            cond,
            &first.gamma_x,
            &second.gamma_x
        )?;

        Ok(Self {
            r_y,
            gamma,
            gamma_x,
        })
    }

    fn cost() -> usize {
        3 * <Fp2G<P> as CondSelectGadget<P::Fp>>::cost()
    }
}

#[derive(Derivative)]
#[derivative(
Clone(bound = "G2Gadget<P>: Clone"),
Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
Debug(bound = "G2Gadget<P>: Debug"),
Debug(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Debug"),
Eq(bound = "G2Gadget<P>: Eq"),
Eq(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Eq"),
PartialEq(bound = "G2Gadget<P>: PartialEq"),
PartialEq(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: PartialEq"),
)]
pub struct G2PreparedGadget<P: MNT4Parameters>{
    pub q:      G2Gadget<P>,
    pub coeffs: Vec<G2CoefficientsGadget<P>>
}

impl<P: MNT4Parameters>G2PreparedGadget<P> {
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

    fn doubling_step_for_flipped_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        s: &G2Gadget<P>,
    ) -> Result<(G2Gadget<P>, G2CoefficientsGadget<P>), SynthesisError>
    {
        //Compute gamma
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

    fn mixed_addition_step_for_flipped_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        x: &Fp2G<P>,
        y: &Fp2G<P>,
        s: &G2Gadget<P>,
        naf_i: i32,
    ) -> Result<(G2Gadget<P>, G2CoefficientsGadget<P>), SynthesisError>
    {
        //Compute gamma
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

impl<P: MNT4Parameters> ConstantGadget<G2Prepared<P>, P::Fp> for G2PreparedGadget<P> {

    #[inline]
    fn from_value<CS: ConstraintSystem<P::Fp>>(mut cs: CS, value: &G2Prepared<P>) -> Self {
        let q = G2Gadget::<P>::from_value(cs.ns(|| "hardcode q"), &value.q.into_projective());
        let coeffs = value.coeffs
            .iter()
            .enumerate()
            .map(|(i, query_i)| {
                G2CoefficientsGadget::<P>::from_value(cs.ns(|| format!("coeff_{}", i)), query_i)
            })
            .collect::<Vec<_>>();
        Self {
            q,
            coeffs,
        }
    }

    fn get_constant(&self) -> G2Prepared<P> {
        G2Prepared::<P>{
            q: self.q.get_value().unwrap().into_affine(),
            coeffs: self.coeffs.iter().map(|coeff| coeff.get_constant()).collect::<Vec<_>>()
        }
    }
}

impl<P: MNT4Parameters> CondSelectGadget<P::Fp> for G2PreparedGadget<P> {
    fn conditionally_select<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self
    ) -> Result<Self, SynthesisError>
    {
        let q = G2Gadget::<P>::conditionally_select(
            cs.ns(|| "q"),
            cond,
            &first.q,
            &second.q
        )?;

        let mut coeffs = Vec::new();
        for (i, (first, second)) in
            first.coeffs.iter().zip(second.coeffs.iter()).enumerate() {
            let val = G2CoefficientsGadget::<P>::conditionally_select(
                cs.ns(|| format!("coeff_{}", i)),
                cond,
                &first,
                &second
            )?;
            coeffs.push(val);
        }

        Ok(Self {
            q,
            coeffs,
        })
    }

    fn cost() -> usize {
        let total_steps = P::WNAF.len() +
            P::WNAF.iter().filter(|&val| *val != 0).collect::<Vec<_>>().len();

        <G2Gadget<P> as CondSelectGadget<P::Fp>>::cost() +
            (total_steps * <G2CoefficientsGadget<P> as CondSelectGadget<P::Fp>>::cost())
    }
}