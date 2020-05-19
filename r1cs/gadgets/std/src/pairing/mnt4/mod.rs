/*
Pairing gadget for MNT4 curves:
    - MNT4PairingGadget, and
    - the implementation of the PairingGadget in alignment to its primitive
      (Ate pairing with flipped Miller loop, using pre-computations).

To do: generic treatment of the sign of the trace  using  ATE_IS_LOOP_COUNT_NEG as well as
the sign of the last chunk component using FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG,
see below.
*/

use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{fields::{fp4::Fp4Gadget, FieldGadget}, groups::curves::short_weierstrass::mnt::mnt4::{G1Gadget, G2Gadget, G1PreparedGadget, G2PreparedGadget},
};

use algebra::{ModelParameters, Fp4, PairingCurve};
use crate::pairing::PairingGadget;
use algebra::curves::models::mnt4::{MNT4p, MNT4Parameters,
                                    G1Affine, G1Projective, G1Prepared,
                                    G2Affine, G2Projective, G2Prepared,
};
use std::marker::PhantomData;

pub mod mnt4753;

pub struct MNT4PairingGadget<P: MNT4Parameters>(PhantomData<P>);

impl<P: MNT4Parameters> PairingGadget<MNT4p<P>, P::Fp> for MNT4PairingGadget<P>
    where
        G1Affine<P>: PairingCurve<
            BaseField = <P::G1Parameters as ModelParameters>::BaseField,
            ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
            Projective = G1Projective<P>,
            PairWith = G2Affine<P>,
            Prepared = G1Prepared<P>,
            PairingResult = Fp4<P::Fp4Params>,
        >,
        G2Affine<P>: PairingCurve<
            BaseField = <P::G2Parameters as ModelParameters>::BaseField,
            ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
            Projective = G2Projective<P>,
            PairWith = G1Affine<P>,
            Prepared = G2Prepared<P>,
            PairingResult = Fp4<P::Fp4Params>,
        >,
{
    type G1Gadget = G1Gadget<P>;
    type G2Gadget = G2Gadget<P>;
    type G1PreparedGadget = G1PreparedGadget<P>;
    type G2PreparedGadget = G2PreparedGadget<P>;
    type GTGadget = Fp4Gadget<P::Fp4Params, P::Fp>;

    /* Miller loop for the computing pairing product of a slice of (P,Q)-values
    (using the precomputed G1- and G2PreparedGadgets)
    */
    fn miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        p: &[Self::G1PreparedGadget],
        q: &[Self::G2PreparedGadget],
    ) -> Result<Self::GTGadget, SynthesisError>
    {
        let mut result = Self::GTGadget::one(cs.ns(|| "one"))?;
        let it = p.iter().zip(q.iter());

        for (i, (ps, qs)) in it.into_iter().enumerate() {

            let mut cs = cs.ns(|| format!("Pair_{}", i));

            let mut f = Self::GTGadget::one(cs.ns(|| "f"))?;

            let mut idx: usize = 0;

            for (j, &n) in P::WNAF.iter().rev().enumerate() {

                let mut cs = cs.ns(|| format!("Iteration_{}", j));

                let c = &qs.coeffs[idx];
                idx += 1;

                /* evaluate the tangent line g_{R,R} at P in F4 (scaled by twist^2) using the
                pre-computed data:
                    g_{R,R}(P) = (y_P - lambda*x_p - d) * twist^2,
                where
                     lambda = gamma * Y/twist,
                     d = (y'-gamma * x')* Y/twist^2,
                with (x',y') being the twist coordinates of R. Thus
                    g_{R,R}(P) = y_p*twist^2 + (gamma*x'- gamma*twist*x_p - y') *Y.
                The scale factor twist^2 from F2 is cancelled out by the final exponentiation.
                */
                let g_rr_at_p_c0 = ps.clone().p_y_twist_squared;

                let mut t = c.gamma.mul_by_constant(cs.ns(|| "double compute gamma_twist"), &P::TWIST)?;
                t.mul_assign_by_fp_gadget(cs.ns(|| "double gamma_twist * ps.p.x"), &ps.p.x)?;
                let g_rr_at_p_c1 = c.gamma_x
                    .sub(cs.ns(|| "gamma_x - r_y"), &c.r_y)?
                    .sub(cs.ns(|| "gamma_x - r_y - t"), &t)?;

                //Compute g_rr_at_p
                let g_rr_at_p = Self::GTGadget::new(g_rr_at_p_c0.clone(), g_rr_at_p_c1);

                //Compute new_f
                f = f.square(cs.ns(|| "f^2"))?.mul_by_023(cs.ns(||"double compute f"), &g_rr_at_p)?;

                if n != 0 {
                    // Addition/substraction step
                    // evaluate chord g_{RQ}(P) in F4 using pre-computed data as above
                    let c = &qs.coeffs[idx];
                    idx += 1;

                    let g_rq_at_p_c0 = ps.clone().p_y_twist_squared;

                    //Compute g_rq_at_p_c1
                    let neg_q_y = qs.q.y.negate(cs.ns(|| "- q.y"))?;
                    let q_y = if n > 0 {qs.clone().q.y} else {neg_q_y};

                    let mut t = c.gamma.mul_by_constant(cs.ns(|| "add compute gamma_twist"), &P::TWIST)?;
                    t.mul_assign_by_fp_gadget(cs.ns(|| "add gamma_twist * ps.p.x"), &ps.p.x)?;
                    let g_rq_at_p_c1 = c.gamma_x
                        .sub(cs.ns(|| "gamma_x - q_y"), &q_y)?
                        .sub(cs.ns(|| "gamma_x - q_y - t"), &t)?;

                    //Compute g_rq_at_p
                    let g_rq_at_p = Self::GTGadget::new(g_rq_at_p_c0, g_rq_at_p_c1);

                    //Compute and check f
                    f = f.mul_by_023(cs.ns(||"add compute f"), &g_rq_at_p)?;
                }
            }

            /*
            CAUTION, if clause missing!
            as in the pairing primitive, unitary inverse if and only if P::ATE_IS_LOOP_COUNT_NEG == TRUE
            */
            f = f.unitary_inverse(cs.ns(|| "f unitary inverse"))?;
            result.mul_in_place(cs.ns(|| format!("mul_assign_{}", i)), &f)?;
        }
        Ok(result)
    }

    fn final_exponentiation<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        value: &Self::GTGadget,
    ) -> Result<Self::GTGadget, SynthesisError>{
        // Final exp first chunk, the "easy" part,
        // using the Frobenius map and value_inv= elt^{-1} to compute value^(q^2-1)
        let value_inv = value.inverse(cs.ns(|| "value_inverse"))?;
        // elt = value^(q^2 - 1)
        let elt = value.clone()
            .frobenius_map(cs.ns(|| "value_frobenius_2"), 2)?
            .mul(cs.ns(|| "value_frobenius_2_div_value"), &value_inv)?;
        // and its inverse for later purpose
        let elt_inv = value_inv.clone()
            .frobenius_map(cs.ns(|| "value_inv_frobenius_2"), 2)?
            .mul(cs.ns(|| "value_inv_frobenius_2_div_value"), &value)?;

        // Final exp last chunk, the "hard part", i.e. the
        // remaining exponentiaton by m_1*q + m_0, m_0 can be signed.
        let elt_q = elt.clone()
            .frobenius_map(cs.ns(|| "elt_q_frobenius_1"), 1)?;

        // exponentiation by m_1 and m_0 using optimized exponentiation for r-th roots of unity
        // elt^{q*m_1}
        let w1_part = elt_q
            .cyclotomic_exp(cs.ns(|| "compute w1"), P::FINAL_EXPONENT_LAST_CHUNK_1)?;

        /* CAUTION, code not generic here
        as in the pairing primitive, depending on P::FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG we have to
        choose either elt or elt_inv to compute w0
        */

        //elt^{m_0}
        let w0_part = elt_inv.clone()
            .cyclotomic_exp(cs.ns(|| "compute w0"),P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)?;

        w1_part.mul(cs.ns(|| "w0 * w1"), &w0_part)

    }

    fn prepare_g1<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        q: &Self::G1Gadget,
    ) -> Result<Self::G1PreparedGadget, SynthesisError>
    {
        Self::G1PreparedGadget::from_affine(cs, q)
    }

    fn prepare_g2<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        q: &Self::G2Gadget,
    ) -> Result<Self::G2PreparedGadget, SynthesisError>
    {
        Self::G2PreparedGadget::from_affine(cs, q)
    }
}