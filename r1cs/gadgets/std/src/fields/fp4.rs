/*
Definition of the degree 4 extension field gadget Fp4Gadget, and implementation of the
following gadgets for it:
    - FieldGadget:
        mul, inverse and mul_equ gadget use Karatsuba multiplication
        square gadget can be improved by one constraint as in Fp2Gadget
        NEqGadget has to be checked if it meets it's purpose by demanding all three components to
        be different.
    - cyclotomic operations gadgets as used by the Ate pairing gadget
    - AllocGadget, CloneGadget, ConstantGadget,
    - PartialEqGadget, ConditionalEqGadget, NEqGadget,
    - CondSelectGadget, TwoBitLookupGadget, ThreeBitNegLookupGadget,
    - ToBitsGadget, FromBitsGadget, ToBytesGadget
*/

use algebra::{fields::{
    fp4::{Fp4, Fp4Parameters},
    Field, Fp2Parameters,
}, PrimeField, Fp2, BigInteger, SquareRootField};
use r1cs_core::{ConstraintSystem, ConstraintVar, SynthesisError};
use std::{borrow::Borrow, marker::PhantomData};

use crate::{prelude::*, Assignment};

type Fp2Gadget<P, ConstraintF> = super::fp2::Fp2Gadget<<P as Fp4Parameters>::Fp2Params, ConstraintF>;

#[derive(Derivative)]
#[derivative(Debug(bound = "P: Fp4Parameters, ConstraintF: PrimeField + SquareRootField"))]
#[must_use]
pub struct Fp4Gadget<P, ConstraintF: PrimeField + SquareRootField>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    pub c0: Fp2Gadget<P, ConstraintF>,
    pub c1: Fp2Gadget<P, ConstraintF>,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P, ConstraintF: PrimeField + SquareRootField> Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    pub fn new(c0: Fp2Gadget<P, ConstraintF>, c1: Fp2Gadget<P, ConstraintF>) -> Self {
        Self {
            c0,
            c1,
            _params: PhantomData,
        }
    }

    /// Multiply a Fp2Gadget by the quadratic nonresidue P::NONRESIDUE which defines the extension
    /// field arithmetics
    #[inline]
    pub fn mul_fp2_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &Fp2Gadget<P, ConstraintF>,
    ) -> Result<Fp2Gadget<P, ConstraintF>, SynthesisError> {
        let new_c0 = Fp2Gadget::<P, ConstraintF>::mul_fp_gadget_by_nonresidue(cs, &fe.c1)?;
        let new_c1 = fe.c0.clone();
        Ok(Fp2Gadget::<P, ConstraintF>::new(new_c0, new_c1))
    }

    /* gadgets for the cyclotomic operations (used in the Ate pairing evaluation)
    */
    #[inline]
    pub fn unitary_inverse<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        let new_c0 = self.c0.clone();
        let new_c1 = self.c1.clone().negate(cs.ns(|| "c1 negation"))?;
        Ok(Self::new(new_c0, new_c1))
    }

    pub fn cyclotomic_square<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        let c1_squared = self.c1.square(cs.ns(||"c1^2"))?;
        let c1_squared_nr = Self::mul_fp2_gadget_by_nonresidue(cs.ns(||"nr * c1^2"), &c1_squared)?;
        let one = Fp2::<P::Fp2Params>::one();

        let c0 = {
            let c1_squared_nr_doubled = c1_squared_nr.double(cs.ns(|| "2(nr*c1^2)"))?;
            c1_squared_nr_doubled.add_constant(cs.ns(|| "2(nr*c1^2) + 1"),  &one)?
        };

        let c1 = {
            let c1_plus_c0 = self.c0.add(cs.ns(||"c1 + c0"), &self.c1)?;
            let c1_plus_c0_squared = c1_plus_c0.square(cs.ns(|| "(c1 + c0)^2"))?;
            c1_plus_c0_squared
                .sub(cs.ns(|| "(c1 + c0)^2 - nr*c1^2"), &c1_squared_nr)?
                .sub(cs.ns(|| "(c1 + c0)^2 - nr*c1^2 - c1^2"), &c1_squared)?
                .sub_constant(cs.ns(|| "(c1 + c0)^2 - nr*c1^2 - c1^2 - 1"), &one)?
        };
        Ok(Self::new(c0, c1))
    }

    #[inline]
    pub fn cyclotomic_exp<CS: ConstraintSystem<ConstraintF>, B: BigInteger>(
        &self,
        mut cs: CS,
        exp: B,
    ) -> Result<Self, SynthesisError> {
        let mut res = Self::one(cs.ns(|| "one"))?;
        let self_inverse = self.unitary_inverse(cs.ns(|| "unitary inverse"))?;
        let mut found_nonzero = false;
        let naf = exp.find_wnaf();

        for (j, &value) in naf.iter().rev().enumerate() {
            if found_nonzero {
                res = res.cyclotomic_square(cs.ns(||format!("res_square_{:?}", j)))?;
            }
            if value != 0 {
                found_nonzero = true;

                if value > 0 {
                    res.mul_in_place(cs.ns(|| format!("res_mul_{:?}", j)), self)?;
                } else {
                    res.mul_in_place(cs.ns(|| format!("res_mul_inverse_{:?}", j)), &self_inverse)?;
                }
            }
        }

        Ok(res)
    }

    /* Optimized Karatsuba multiplication of Fp4 gadgets a = a0 + Y*a1 by Fp2 gadgets of the form
            b = b0 + Y*b1  with b0 = (b0.0, 0),
       saves one multiplication. Here, Fp4=Fp2[Y]/(Y^2-beta).
    */
    #[inline]
    pub fn mul_by_023<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        // optimized v0 = (b0.0*a0.0, b0.0*a0.1)
        let v0 =
        {
            let v0_c0 = self.c0.c0.mul(cs.ns(|| "self.c0.c0 * other.c0.c0"), &other.c0.c0)?;
            let v0_c1 = self.c0.c1.mul(cs.ns(|| "self.c0.c1 * other.c0.c0"), &other.c0.c0)?;
            Fp2Gadget::<P, ConstraintF>::new(v0_c0, v0_c1)
        };
        // v1 = a1*b1
        let v1 = self.c1.mul(cs.ns(|| "self.c1 * other.c1"), &other.c1)?;
        // c0 = v0 + beta * v1
        let c0 = {
            let non_residue_times_v1 =
                Self::mul_fp2_gadget_by_nonresidue(cs.ns(|| "v1 mul_by_nr"), &v1)?;
            v0.add(cs.ns(|| "v0 + beta * v1"), &non_residue_times_v1)?
        };
        // c1 = (a0 + a1) * (b0 + b1) - v0 -v1.
        let c1 = {
            let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
            let b0_plus_b1 = other.c0.add(cs.ns(|| "b0 + b1"), &other.c1)?;
            let a0_plus_a1_times_b0_plus_b1 =
                a0_plus_a1.mul(&mut cs.ns(|| "(a0 + a1) * (b0 + b1)"), &b0_plus_b1)?;
            a0_plus_a1_times_b0_plus_b1
                .sub(cs.ns(|| "res - v0"), &v0)?
                .sub(cs.ns(|| "res - v0 - v1"), &v1)?
        };

        Ok(Self::new(c0, c1))
    }
}

/* FieldGadget implementation for Fp4Gadget as quadratic extension of Fp2
*/
impl<P, ConstraintF: PrimeField + SquareRootField> FieldGadget<Fp4<P>, ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    type Variable = (
        (ConstraintVar<ConstraintF>, ConstraintVar<ConstraintF>),
        (ConstraintVar<ConstraintF>, ConstraintVar<ConstraintF>),
    );

    #[inline]
    fn get_value(&self) -> Option<Fp4<P>> {
        match (
            self.c0.get_value(),
            self.c1.get_value(),
        ) {
            (Some(c0), Some(c1)) => Some(Fp4::new(c0, c1)),
            (..) => None,
        }
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        (
            self.c0.get_variable(),
            self.c1.get_variable(),
        )
    }

    #[inline]
    fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = Fp2Gadget::<P, ConstraintF>::zero(cs.ns(|| "c0"))?;
        let c1 = Fp2Gadget::<P, ConstraintF>::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = Fp2Gadget::<P, ConstraintF>::one(cs.ns(|| "c0"))?;
        let c1 = Fp2Gadget::<P, ConstraintF>::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    /*
    addition gadgets
    */

    #[inline]
    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.add(&mut cs.ns(|| "add c0"), &other.c0)?;
        let c1 = self.c1.add(&mut cs.ns(|| "add c1"), &other.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Fp4<P>,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.add_constant(cs.ns(|| "c0"), &other.c0)?;
        let c1 = self.c1.add_constant(cs.ns(|| "c1"), &other.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        other: &Fp4<P>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.add_constant_in_place(cs.ns(|| "c0"), &other.c0)?;
        self.c1.add_constant_in_place(cs.ns(|| "c1"), &other.c1)?;
        Ok(self)
    }

    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bit: &Boolean,
        coeff: Fp4<P>,
    ) -> Result<Self, SynthesisError> {
        let c0 = self
            .c0
            .conditionally_add_constant(cs.ns(|| "c0"), bit, coeff.c0)?;
        let c1 = self
            .c1
            .conditionally_add_constant(cs.ns(|| "c1"), bit, coeff.c1)?;
        Ok(Self::new(c0, c1))
    }

    /*
    substraction gadgets
    */

    #[inline]
    fn sub<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.sub(&mut cs.ns(|| "sub c0"), &other.c0)?;
        let c1 = self.c1.sub(&mut cs.ns(|| "sub c1"), &other.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn negate<CS: ConstraintSystem<ConstraintF>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = self.c0.negate(&mut cs.ns(|| "negate c0"))?;
        let c1 = self.c1.negate(&mut cs.ns(|| "negate c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn negate_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.negate_in_place(&mut cs.ns(|| "negate c0"))?;
        self.c1.negate_in_place(&mut cs.ns(|| "negate c1"))?;
        Ok(self)
    }

    /*
    multiplication gadgets
    */

    #[inline]
    fn mul<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {

        /*
          Karatsuba multiplication for Fp4 as a quadratic extension of Fp2:
          v0 = A.c0 * B.c0
          v1 = A.c1 * B.c1
          result.c0 = v0 + non_residue * v1
          result.c1 = (A.c0 + A.c1) * (B.c0 + B.c1) - v0 - v1
          where "non_residue * elem" := (non_residue * elt.c1, elt.c0)

          Enforced with 3 Fp2_mul_gadget's that ensure that:
          A.c1 * B.c1 = v1
          A.c0 * B.c0 = v0
          (A.c0+A.c1)*(B.c0+B.c1) = result.c1 + v0 + v1

          Reference:
          "Multiplication and Squaring on Pairing-Friendly Fields"
          Devegili, OhEigeartaigh, Scott, Dahab
        */

        let v0 = self.c0.mul(cs.ns(|| "v0"), &other.c0)?;
        let v1 = self.c1.mul(cs.ns(|| "v1"), &other.c1)?;
        let c0 = {
            let non_residue_times_v1 =
                Self::mul_fp2_gadget_by_nonresidue(cs.ns(|| "first mul_by_nr"), &v1)?;
            v0.add(cs.ns(|| "v0 + beta * v1"), &non_residue_times_v1)?
        };

        let c1 = {
            let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
            let b0_plus_b1 = other.c0.add(cs.ns(|| "b0 + b1"), &other.c1)?;
            let a0_plus_a1_times_b0_plus_b1 =
                a0_plus_a1.mul(&mut cs.ns(|| "(a0 + a1) * (b0 + b1)"), &b0_plus_b1)?;
            a0_plus_a1_times_b0_plus_b1
                .sub(cs.ns(|| "res - v0"), &v0)?
                .sub(cs.ns(|| "res - v0 - v1"), &v1)?
        };

        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Fp4<P>,
    ) -> Result<Self, SynthesisError> {
        /* ordinary complex multiplication
            c0 + c1 *X = (a0 + a1 * X) * (b0 + b1 * X)
           of a field gadget self =  (a0 + a1 * X) by a constant element (b0 + b1 * X).
            c0 = a0*b0 + non_residue *a1*b1
            c1 = a0*b1 + a1*b0
        Doesn't need any constraints; returns linear combinations of a0, a1.
        */
        let (a0, a1) = (&self.c0, &self.c1);
        let (b0, b1) = (other.c0, other.c1);
        let mut v0 = a0.mul_by_constant(&mut cs.ns(|| "v0"), &b0)?;
        // misleading names
        let mut v1 = Self::mul_fp2_gadget_by_nonresidue(&mut cs.ns(|| "beta * a1"), a1)?;
        let beta_v1 = v1.mul_by_constant_in_place(&mut cs.ns(|| "beta * a1*b1"), &b1)?;

        v0.add_in_place(&mut cs.ns(|| "c0"), beta_v1)?;
        let c0 = v0;

        let mut a0b1 = a0.mul_by_constant(&mut cs.ns(|| "a0b1"), &b1)?;
        let a1b0 = a1.mul_by_constant(&mut cs.ns(|| "a1b0"), &b0)?;
        a0b1.add_in_place(&mut cs.ns(|| "c1"), &a1b0)?;
        let c1 = a0b1;
        Ok(Self::new(c0, c1))
    }

    // can be improved by one constraint over Fp2 by using the same trick as for the Fp2Gadget
    #[inline]
    fn square<CS: ConstraintSystem<ConstraintF>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
        /*
        Karatsuba squaring for Fp4 as a quadratic extension of Fp2:
        v0 = A.c0^2
        v1 = A.c1^2
        result.c0 = v0 + non_residue * v1
        result.c1 = (A.c0 + A.c1)^2 - v0 - v1

        Enforced with 3 Fp2_sqr_gadget's that ensure that:
        A.c1^2 = v1
        A.c0^2 = v0
        (A.c0+A.c1)^2 = result.c1 + v0 + v1

        Reference:
        "Multiplication and Squaring on Pairing-Friendly Fields"
        Devegili, OhEigeartaigh, Scott, Dahab
        */
        let v0 = self.c0.square(cs.ns(|| "v0 = a0^2"))?;
        let v1 = self.c1.square(cs.ns(|| "v1 = a1^2"))?;
        let c0 = {
            let non_residue_times_v1 =
                Self::mul_fp2_gadget_by_nonresidue(cs.ns(|| "first mul_by_nr"), &v1)?;
            v0.add(cs.ns(|| "v0 + beta * v1"), &non_residue_times_v1)?
        };

        let c1 = {
            let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
            let a0_plus_a1_squared = a0_plus_a1.square(cs.ns(|| "(a0 + a1)^2"))?;
            a0_plus_a1_squared
                .sub(cs.ns(|| "res - v0"), &v0)?
                .sub(cs.ns(|| "res - v0 - v1"), &v1)?
        };

        Ok(Self::new(c0, c1))
    }

    // now that the patched code does not save any constraints:
    // why not replace Karatsuba code by simple call to mul_equals?
    #[inline]
    fn inverse<CS: ConstraintSystem<ConstraintF>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
        let inverse = Self::alloc(&mut cs.ns(|| "alloc inverse"), || {
            self.get_value().and_then(|val| val.inverse()).get()
        })?;

        // Karatsuba multiplication for Fp2 with the inverse:
        //     v0 = A.c0 * B.c0,
        //     v1 = A.c1 * B.c1,
        //      1 = v0 + non_residue * v1,
        //      0 = result.c1 = (A.c0 + A.c1) * (B.c0 + B.c1) - v0 - v1.
        // Enforced with 3 constraints (substituting v0 by v1)
        //    (1)  A.c1 * B.c1 = v1,
        //    (2) (A.c0 + A.c1) * (B.c0 + B.c1) =  1 + (1 - non_residue) * v1
        //                                      = 1 - non_residue * v1 + v1
        //    (3)  A.c0 * B.c0 = 1 - non_residue * v1,
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab

        // Constraint 1
        let v1 = self.c1.mul(cs.ns(|| "inv_constraint_1"), &inverse.c1)?;

        // Constraint 2
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
        let b0_plus_b1 = inverse.c0.add(cs.ns(|| "b0 + b1"), &inverse.c1)?;

        let one = Fp2::<P::Fp2Params>::one();
        let rhs = Self::mul_fp2_gadget_by_nonresidue(cs.ns(|| "nr * v1"), &v1)?
            .sub(cs.ns(|| "sub v1"), &v1)?
            .negate(cs.ns(|| "negate it"))?
            .add_constant(cs.ns(|| "add one"), &one)?;
        a0_plus_a1.mul_equals(cs.ns(|| "inv_constraint_2"), &b0_plus_b1, &rhs)?;

        // Constraint 3
        let rhs2 = rhs.sub(cs.ns(|| " 1 - nonresidue * v1"), &v1)?;
        self.c0.mul_equals(cs.ns(||"inv_constraint_3"),&inverse.c0, &rhs2)?;

        Ok(inverse)
    }

    // does not save any constraints over default implementation?
    fn mul_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        // Karatsuba multiplication for Fp2:
        //     v0 = A.c0 * B.c0
        //     v1 = A.c1 * B.c1
        //     result.c0 = v0 + non_residue * v1
        //     result.c1 = (A.c0 + A.c1) * (B.c0 + B.c1) - v0 - v1
        // Enforced with 3 constraints:
        //     A.c1 * B.c1 = v1
        //     A.c0 * B.c0 = result.c0 - non_residue * v1
        //     (A.c0+A.c1)*(B.c0+B.c1) = result.c1 + result.c0 + (1 - non_residue) * v1
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab
        let mul_cs = &mut cs.ns(|| "mul");

        // Compute v1
        let v1 = self.c1.mul(mul_cs.ns(|| "v1"), &other.c1)?;

        // Perform second check
        let non_residue_times_v1 = Self::mul_fp2_gadget_by_nonresidue(mul_cs.ns(|| "nr * v1"), &v1)?;
        let rhs = result
            .c0
            .sub(mul_cs.ns(|| "sub from result.c0"), &non_residue_times_v1)?;
        self.c0
            .mul_equals(mul_cs.ns(|| "second check"), &other.c0, &rhs)?;

        // Last check
        let a0_plus_a1 = self.c0.add(mul_cs.ns(|| "a0 + a1"), &self.c1)?;
        let b0_plus_b1 = other.c0.add(mul_cs.ns(|| "b0 + b1"), &other.c1)?;
        let one_minus_non_residue_v1 =
            v1.sub(mul_cs.ns(|| "sub from v1"), &non_residue_times_v1)?;

        let result_c1_plus_result_c0_plus_one_minus_non_residue_v1 = result
            .c1
            .add(mul_cs.ns(|| "c1 + c0"), &result.c0)?
            .add(mul_cs.ns(|| "rest of stuff"), &one_minus_non_residue_v1)?;

        a0_plus_a1.mul_equals(
            mul_cs.ns(|| "third check"),
            &b0_plus_b1,
            &result_c1_plus_result_c0_plus_one_minus_non_residue_v1,
        )?;

        Ok(())
    }

    fn frobenius_map<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        power: usize,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.frobenius_map_in_place(cs, power)?;
        Ok(result)
    }

    fn frobenius_map_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.frobenius_map_in_place(&mut cs.ns(|| "c0"), power)?;
        self.c1.frobenius_map_in_place(&mut cs.ns(|| "c1"), power)?;

        self.c1.c0.mul_by_constant_in_place(
            cs.ns(|| "c1_c0_power"),
            &P::FROBENIUS_COEFF_FP4_C1[power % 4],
        )?;
        self.c1.c1.mul_by_constant_in_place(
            cs.ns(|| "c1_c1_power"),
            &P::FROBENIUS_COEFF_FP4_C1[power % 4],
        )?;

        Ok(self)
    }

    fn cost_of_mul() -> usize {
        3 * Fp2Gadget::<P, ConstraintF>::cost_of_mul()
    }

    fn cost_of_mul_equals() -> usize {
        3 * Fp2Gadget::<P, ConstraintF>::cost_of_mul_equals()
    }

    fn cost_of_inv() -> usize {
        1 * Fp2Gadget::<P,ConstraintF>::cost_of_mul()
            + 2 * Fp2Gadget::<P, ConstraintF>::cost_of_mul_equals()
    }
}


/*
Alloc-, Clone and ConstantGadget for the Fp2Gadget
*/

impl<P, ConstraintF: PrimeField + SquareRootField> AllocGadget<Fp4<P>, ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<Fp4<P>>,
    {
        let (c0, c1) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1))
            },
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = Fp2Gadget::<P, ConstraintF>::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = Fp2Gadget::<P, ConstraintF>::alloc(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<Fp4<P>>,
    {
        let (c0, c1) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1))
            },
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = Fp2Gadget::<P, ConstraintF>::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = Fp2Gadget::<P, ConstraintF>::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> Clone for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    fn clone(&self) -> Self {
        Self::new(self.c0.clone(), self.c1.clone())
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> ConstantGadget<Fp4<P>, ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn from_value<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value: &Fp4<P>,
    ) -> Self
    {
        let c0 = Fp2Gadget::<P, ConstraintF>::from_value(&mut cs.ns(|| "c0"), &value.c0);
        let c1 = Fp2Gadget::<P, ConstraintF>::from_value(&mut cs.ns(|| "c1"), &value.c1);
        Self::new(c0, c1)
    }

    #[inline]
    fn get_constant(&self) -> Fp4<P> {
        self.get_value().unwrap()
    }
}

/*
relational and conditional gadgets (incl. lookup tables) for Fp4Gadgets
*/

impl<P, ConstraintF: PrimeField + SquareRootField> PartialEq for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> Eq for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
}

impl<P, ConstraintF: PrimeField + SquareRootField> EqGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
}

impl<P, ConstraintF: PrimeField + SquareRootField> ConditionalEqGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.c0
            .conditional_enforce_equal(&mut cs.ns(|| "c0"), &other.c0, condition)?;
        self.c1
            .conditional_enforce_equal(&mut cs.ns(|| "c1"), &other.c1, condition)?;
        Ok(())
    }


    fn cost() -> usize {
        2 * <Fp2Gadget<P, ConstraintF> as ConditionalEqGadget<ConstraintF>>::cost()
    }
}

/* enforces that all two components of two Fp4Gadgets are not equal.
Note: This is not the canonical notion of inequality for field elements, we need to check
whether the implementation depicts the need of futher application
*/
impl<P, ConstraintF: PrimeField + SquareRootField> NEqGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.c0.enforce_not_equal(&mut cs.ns(|| "c0"), &other.c0)?;
        self.c1.enforce_not_equal(&mut cs.ns(|| "c1"), &other.c1)?;
        Ok(())
    }

    fn cost() -> usize {
        2 * <Fp2Gadget<P, ConstraintF> as NEqGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> CondSelectGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = Fp2Gadget::<P, ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c0"),
            cond,
            &first.c0,
            &second.c0,
        )?;
        let c1 = Fp2Gadget::<P, ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c1"),
            cond,
            &first.c1,
            &second.c1,
        )?;

        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <Fp2Gadget<P, ConstraintF> as CondSelectGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> TwoBitLookupGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    type TableConstant = Fp4<P>;
    /* Two bit lookup circuit, chooses c[i], i=0,..,3, according to i = b[1]*2 + b[0] given
    by two Booleans b[0], b[1].
    */
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = Fp2Gadget::<P, ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = Fp2Gadget::<P, ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    /* low level version of two bit lookup table, outsources the quadratic term b0 * b1 in
       (c0-c1-c2+c3) * b0 * b1 = result - [c0 + (c1-c0)*b0 + (c2-c0)*b1]
    to an input Boolean 'precomp'.
    */
    fn two_bit_lookup_lc<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        precomp: &Boolean,
        b: &[Boolean],
        c: &[Self::TableConstant])
        -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = Fp2Gadget::<P, ConstraintF>::two_bit_lookup_lc(cs.ns(|| "Lookup c0"), precomp, b, &c0s)?;
        let c1 = Fp2Gadget::<P, ConstraintF>::two_bit_lookup_lc(cs.ns(|| "Lookup c1"), precomp, b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <Fp2Gadget<P, ConstraintF> as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> ThreeBitCondNegLookupGadget<ConstraintF>
for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    type TableConstant = Fp4<P>;

    /* Three bit look up circuit combining a two bit table and conditional negation.
    Outputs
        c = (-1)^b[2] * c[i] = (-1)^b[2] * c[b[1]*2 + b[0]]
    according to the bits b[0],b[1],b[2].
    As two_bit_lookup_lc, the quadratic term b[0]*[1] is outsourced to the input Boolean b0b1.
    */
    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = Fp2Gadget::<P, ConstraintF>::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c0"), b, b0b1, &c0s)?;
        let c1 = Fp2Gadget::<P, ConstraintF>::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c1"), b, b0b1, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <Fp2Gadget<P, ConstraintF> as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

/*
Packing and unpacking gadgets
*/
impl<P, ConstraintF: PrimeField + SquareRootField> ToBitsGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, mut cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits(&mut cs)?;
        let mut c1 = self.c1.to_bits(&mut cs)?;
        c0.append(&mut c1);

        Ok(c0)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits_strict(&mut cs)?;
        let mut c1 = self.c1.to_bits_strict(&mut cs)?;

        c0.append(&mut c1);

        Ok(c0)
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> ToBytesGadget<ConstraintF> for Fp4Gadget<P, ConstraintF>
    where
        P: Fp4Parameters,
        P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_bytes(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bytes(cs.ns(|| "c1"))?;

        c0.append(&mut c1);

        Ok(c0)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_bytes_strict(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bytes_strict(cs.ns(|| "c1"))?;

        c0.append(&mut c1);

        Ok(c0)
    }
}

