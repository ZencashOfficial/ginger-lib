/*
Definition of the quadratic extension field gadget Fp2Gadget and implementation of the
following traits for it:
    - FieldGadget:
        mul and related gadgets using Karatsuba multiplication,
        extra implementations for square and square_in_place saving one constraint,
        NEqGadget has to be checked if it meets it's purpose by demanding all two components to be
        different.
    - AllocGadget, CloneGadget, ConstantGadget,
    - PartialEqGadget, ConditionalEqGadget, NEqGadget,
    - CondSelectGadget, TwoBitLookupGadget, ThreeBitNegLookupGadget,
    - ToBitsGadget, FromBitsGadget, ToBytesGadget
*/

use algebra::{
    fields::{Fp2, Fp2Parameters},
    Field, PrimeField, SquareRootField,
};
use r1cs_core::{ConstraintSystem, ConstraintVar, SynthesisError};
use std::{borrow::Borrow, marker::PhantomData};

use crate::{fields::fp::FpGadget, prelude::*, Assignment};

#[derive(Derivative)]
#[derivative(Debug(bound = "P: Fp2Parameters, ConstraintF: PrimeField + SquareRootField"))]
#[must_use]
pub struct Fp2Gadget<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> {
    pub c0: FpGadget<ConstraintF>,
    pub c1: FpGadget<ConstraintF>,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

/* extra extension field gadgets for the Fp2Gadget
*/
impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> Fp2Gadget<P, ConstraintF> {
    pub fn new(c0: FpGadget<ConstraintF>, c1: FpGadget<ConstraintF>) -> Self {
        Self {
            c0,
            c1,
            _params: PhantomData,
        }
    }

    /// Multiply a FpGadget by quadratic nonresidue P::NONRESIDUE which defines the Fp2 arithemtics
    #[inline]
    pub fn mul_fp_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &FpGadget<ConstraintF>,
    ) -> Result<FpGadget<ConstraintF>, SynthesisError> {
        fe.mul_by_constant(cs, &P::NONRESIDUE)
    }

    /// Multiply a Fp2Gadget by a FpGadget.
    #[inline]
    pub fn mul_assign_by_fp_gadget<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &FpGadget<ConstraintF>,
    ) -> Result<&mut Self, SynthesisError>
    {
        self.c0.mul_in_place(cs.ns(||"compute new_c0"), &fe)?;
        self.c1.mul_in_place(cs.ns(||"compute new_c1"), &fe)?;
        Ok(self)
    }

    /// Multiply a Fp2Gadget by an element of fp.
    #[inline]
    pub fn mul_by_fp_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &P::Fp,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.mul_by_constant_in_place(cs.ns(|| "c0"), fe)?;
        self.c1.mul_by_constant_in_place(cs.ns(|| "c1"), fe)?;
        Ok(self)
    }

    /// Multiply a Fp2Gadget by an element of fp.
    #[inline]
    pub fn mul_by_fp_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        fe: &P::Fp,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.mul_by_fp_constant_in_place(cs, fe)?;
        Ok(result)
    }
}

/* FieldGadget implementation for the Fp2Gadget
*/
impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> FieldGadget<Fp2<P>, ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    type Variable = (ConstraintVar<ConstraintF>, ConstraintVar<ConstraintF>);

    #[inline]
    fn get_value(&self) -> Option<Fp2<P>> {
        match (self.c0.value, self.c1.value) {
            (Some(c0), Some(c1)) => Some(Fp2::new(c0, c1)),
            (..) => None,
        }
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        (
            self.c0.get_variable().clone(),
            self.c1.get_variable().clone(),
        )
    }

    #[inline]
    fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::zero(cs.ns(|| "c0"))?;
        let c1 = FpGadget::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::one(cs.ns(|| "c0"))?;
        let c1 = FpGadget::zero(cs.ns(|| "c1"))?;
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
        cs: CS,
        other: &Fp2<P>,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        let _ = result.add_constant_in_place(cs, other)?;
        Ok(result)
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        other: &Fp2<P>,
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
        coeff: Fp2<P>,
    ) -> Result<Self, SynthesisError> {
        let c0 = self
            .c0
            .conditionally_add_constant(cs.ns(|| "c0"), bit, coeff.c0)?;
        let c1 = self
            .c1
            .conditionally_add_constant(cs.ns(|| "c1"), bit, coeff.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn double<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.double_in_place(cs)?;
        Ok(result)
    }

    #[inline]
    fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.double_in_place(&mut cs.ns(|| "double c0"))?;
        self.c1.double_in_place(&mut cs.ns(|| "double c1"))?;
        Ok(self)
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
    fn negate<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.negate_in_place(cs)?;
        Ok(result)
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
        /* Karatsuba multiplication for Fp2:
             v0 = A.c0 * B.c0
             v1 = A.c1 * B.c1
             result.c0 = v0 + non_residue * v1
             result.c1 = (A.c0 + A.c1) * (B.c0 + B.c1) - v0 - v1.
        Reference:
        "Multiplication and Squaring on Pairing-Friendly Fields"
         Devegili, OhEigeartaigh, Scott, Dahab
        Can be enforced with 3 constraints (but not done in the code below, which uses four constr.)
             A.c1 * B.c1 = v1
             A.c0 * B.c0 = result.c0 - non_residue * v1
             (A.c0+A.c1)*(B.c0+B.c1) = result.c1 + result.c0 + (1 - non_residue) * v1
        */
        let mul_cs = &mut cs.ns(|| "mul");

        let v0 = self.c0.mul(mul_cs.ns(|| "v0"), &other.c0)?;
        let v1 = self.c1.mul(mul_cs.ns(|| "v1"), &other.c1)?;
        let c0 = {
            let non_residue_times_v1 =
                v1.mul_by_constant(mul_cs.ns(|| "non_residue * v1"), &P::NONRESIDUE)?;
            v0.add(mul_cs.ns(|| "v0 + non_residue * v1"), &non_residue_times_v1)?
        };
        let c1 = {
            let a0_plus_a1 = self.c0.add(mul_cs.ns(|| "a0 + a1"), &self.c1)?;
            let b0_plus_b1 = other.c0.add(mul_cs.ns(|| "b0 + b1"), &other.c1)?;
            let a0_plus_a1_times_b0_plus_b1 =
                a0_plus_a1.mul(&mut mul_cs.ns(|| "(a0 + a1) * (b0 + b1)"), &b0_plus_b1)?;
            a0_plus_a1_times_b0_plus_b1
                .sub(mul_cs.ns(|| "(a0 + a1) * (b0 + b1) - v0"), &v0)?
                .sub(mul_cs.ns(|| "(a0 + a1) * (b0 + b1) - v0 - v1"), &v1)?
        };
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        fe: &Fp2<P>,
    ) -> Result<Self, SynthesisError> {
        /* ordinary complex multiplication
            c0 + c1 *X = (a0 + a1 * X) * (b0 + b1 * X)
           of a field gadget self =  (a0 + a1 * X) by a constant element (b0 + b1 * X).
            c0 = a0*b0 + non_residue *a1*b1
            c1 = a0*b1 + a1*b0
        Doesn't need any constraints; returns linear combinations of a0, a1.
        */
        let (a0, a1) = (&self.c0, &self.c1);
        let (b0, b1) = (fe.c0, fe.c1);
        let mut v0 = a0.mul_by_constant(&mut cs.ns(|| "v0"), &b0)?;
        let beta_v1 = a1.mul_by_constant(&mut cs.ns(|| "non_residue * v1"), &(b1 * &P::NONRESIDUE))?;

        v0.add_in_place(&mut cs.ns(|| "c0 = v0 + non_residue * v1"), &beta_v1)?;
        let c0 = v0;

        let mut a0b1 = a0.mul_by_constant(&mut cs.ns(|| "a0b1"), &b1)?;
        let a1b0 = a1.mul_by_constant(&mut cs.ns(|| "a1b0"), &b0)?;
        a0b1.add_in_place(&mut cs.ns(|| "c1 = a0b1 + a1b0"), &a1b0)?;
        let c1 = a0b1;
        Ok(Self::new(c0, c1))
    }

    /* improves default implementation by 1 constraint.
    */
    #[inline]
    fn square<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        /* From Libsnark/fp2_gadget.tcc
         Complex squaring for Fp2,
             result.c0 = A.c0^2 + A.c1^2 * non_residue,
             result.c1 = 2 * A.c0 * A.c1
         Using the Karatsuba trick, enforced by two rank 1 constraints,
             A.c0 * A.c1 = v0,
             (A.c0 + A.c1) * (A.c0 + non_residue * A.c1) - v0 * (1 + non_residue) = result.c0,
         and  2*v0 = result.c1 operating directly on LC level.

        Reference:
        "Multiplication and Squaring on Pairing-Friendly Fields"
        Devegili, OhEigeartaigh, Scott, Dahab
        */
        // constraint 1
        let mut v0 = self.c0.mul(cs.ns(|| "v0"), &self.c1)?;
        // prepare for constraint 2
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
        let non_residue_c1 = self
            .c1
            .mul_by_constant(cs.ns(|| "non_residue * a1"), &P::NONRESIDUE)?;
        let a0_plus_non_residue_c1 = self
            .c0
            .add(cs.ns(|| "a0 + non_residue * a1"), &non_residue_c1)?;
        let one_plus_non_residue_v0 = v0.mul_by_constant(
            cs.ns(|| "(1 + non_residue) * v0"),
            &(P::Fp::one() + &P::NONRESIDUE),
        )?;
        // constraint 2
        let c0 = a0_plus_a1
            .mul(
                cs.ns(|| "(a0 + a1) * (a0 + non_residue * a1)"),
                &a0_plus_non_residue_c1,
            )?
            .sub(cs.ns(|| "- (1 + non_residue) v0"), &one_plus_non_residue_v0)?;

        // no extra constraint, works directly works on LC:
        v0.double_in_place(cs.ns(|| "2v0"))?;
        let c1 = v0;

        Ok(Self::new(c0, c1))
    }

    /* no improvement in number of constraints,
    improves default implementation by one private variable?
    */
    #[inline]
    fn square_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        /* From Libsnark/fp2_gadget.tcc
         Complex squaring for Fp2,
             result.c0 = A.c0^2 + A.c1^2 * non_residue,
             result.c1 = 2 * A.c0 * A.c1
         Using the Karatsuba trick, enforced by two rank 1 constraints,
             A.c0 * A.c1 = v0,
             (A.c0 + A.c1) * (A.c0 + non_residue * A.c1) - v0 * (1 + non_residue) = result.c0,
         and  2*v0 = result.c1 operating directly on LC level.
        Reference:
        "Multiplication and Squaring on Pairing-Friendly Fields"
        Devegili, OhEigeartaigh, Scott, Dahab
        */
        // constraint 1
        let mut v0 = self.c0.mul(cs.ns(|| "v0"), &self.c1)?;
        // prepare for constraint 2
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;

        let _ = self
            .c1
            .mul_by_constant_in_place(cs.ns(|| "non_residue * a1"), &P::NONRESIDUE)?;
        let a0_plus_non_residue_c1 = self
            .c0
            .add(cs.ns(|| "a0 + non_residue * a1"), &self.c1)?;
        let one_plus_non_residue_v0 = v0.mul_by_constant(
            cs.ns(|| "(1 + non_residue) * v0"),
            &(P::Fp::one() + &P::NONRESIDUE),
        )?;
        // constraint 2
        self.c0 = a0_plus_a1
            .mul(
                cs.ns(|| "(a0 + a1) * (a0 + non_residue * a1)"),
                &a0_plus_non_residue_c1,
            )?
            .sub(cs.ns(|| "- (1 + non_residue) * v0"), &one_plus_non_residue_v0)?;

        // no extra constraint, works directly on LC
        v0.double_in_place(cs.ns(|| "2v0"))?;
        self.c1 = v0;

        Ok(self)
    }


    // now that the patched code does not save any constraints:
    // why not replace Karatsuba code by simple call to mul_equals?
    #[inline]
    fn inverse<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
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
        let v1 = self.c1.mul(cs.ns(|| "inv_constraint_1"),
                             &inverse.c1)?;

        // Constraint 2
        let one = P::Fp::one();
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
        let b0_plus_b1 = inverse.c0.add(cs.ns(|| "b0 + b1"), &inverse.c1)?;

        let rhs =
            v1.mul_by_constant(
                cs.ns(|| "(1 - nonresidue) * v1"),
                &(one - &P::NONRESIDUE)
            )?
                .add_constant(
                    cs.ns(|| "add one"),
                    &one)?;

        a0_plus_a1.mul_equals(cs.ns(|| "inv_constraint_2"), &b0_plus_b1, &rhs)?;

        // Constraint 3
        let rhs2 = rhs.sub(cs.ns(|| " 1 - nonresidue * v1"), &v1)?;
        self.c0.mul_equals(cs.ns(||"inv_constraint_3"),&inverse.c0, &rhs2)?;

        Ok(inverse)
    }

    // does not save any constraint over default implementation?
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
        let mut v1 = self.c1.mul(mul_cs.ns(|| "v1"), &other.c1)?;

        // Perform second check
        let non_residue_times_v1 =
            v1.mul_by_constant(mul_cs.ns(|| "non_residue * v0"), &P::NONRESIDUE)?;
        let rhs = result
            .c0
            .sub(mul_cs.ns(|| "sub from result.c0"), &non_residue_times_v1)?;
        self.c0
            .mul_equals(mul_cs.ns(|| "second check"), &other.c0, &rhs)?;

        // Last check
        let a0_plus_a1 = self.c0.add(mul_cs.ns(|| "a0 + a1"), &self.c1)?;
        let b0_plus_b1 = other.c0.add(mul_cs.ns(|| "b0 + b1"), &other.c1)?;
        let one_minus_non_residue_v1 =
            v1.sub_in_place(mul_cs.ns(|| "sub from v1"), &non_residue_times_v1)?;

        let result_c1_plus_result_c0_plus_one_minus_non_residue_v1 = result
            .c1
            .add(mul_cs.ns(|| "c1 + c0"), &result.c0)?
            .add(mul_cs.ns(|| "rest of stuff"), one_minus_non_residue_v1)?;

        a0_plus_a1.mul_equals(
            mul_cs.ns(|| "third check"),
            &b0_plus_b1,
            &result_c1_plus_result_c0_plus_one_minus_non_residue_v1,
        )?;

        Ok(())
    }

    #[inline]
    fn frobenius_map<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        power: usize,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        let _ = result.frobenius_map_in_place(cs, power)?;
        Ok(result)
    }

    /* k-th power of the Frobenius map x->x^p in Fp2
       pi^k(c0 + c1*X) = c0 + c1 *pi^k(X) = c0 + c1*FROBENIUS_COEFF[k%2]* X
    */
    #[inline]
    fn frobenius_map_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        self.c1
            .mul_by_constant_in_place(cs, &P::FROBENIUS_COEFF_FP2_C1[power % 2])?;
        Ok(self)
    }

    fn cost_of_mul() -> usize {
        3
    }

    fn cost_of_mul_equals() -> usize {
        3
    }

    fn cost_of_inv() -> usize {
        3
    }
}

/*
Alloc-, Clone and ConstantGadget for the Fp2Gadget
*/

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> AllocGadget<Fp2<P>, ConstraintF>
for Fp2Gadget<P, ConstraintF>
{
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<Fp2<P>>,
    {
        let (c0, c1) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1))
            },
            Err(_) => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = FpGadget::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = FpGadget::alloc(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<Fp2<P>>,
    {
        let (c0, c1) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1))
            },
            Err(_) => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = FpGadget::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = FpGadget::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> Clone
for Fp2Gadget<P, ConstraintF>
{
    fn clone(&self) -> Self {
        Self {
            c0:      self.c0.clone(),
            c1:      self.c1.clone(),
            _params: PhantomData,
        }
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
ConstantGadget<Fp2<P>, ConstraintF> for Fp2Gadget<P, ConstraintF>
{
    #[inline]
    fn from_value<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value: &Fp2<P>,
    ) -> Self
    {
        let c0 = FpGadget::<ConstraintF>::from_value(&mut cs.ns(|| "c0"), &value.c0);
        let c1 = FpGadget::<ConstraintF>::from_value(&mut cs.ns(|| "c1"), &value.c1);
        Self::new(c0, c1)
    }

    #[inline]
    fn get_constant(&self) -> Fp2<P> {
        self.get_value().unwrap()
    }
}

/*
relational and conditional gadgets (incl. lookup tables) for the Fp2Gadget
*/

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> PartialEq
    for Fp2Gadget<P, ConstraintF>
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> Eq for Fp2Gadget<P, ConstraintF> {}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> EqGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> ConditionalEqGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    /* enforces equality of two Fp2Gadgets x,y if the condition Boolean is true, otherwise it does
    not enforce anything.
    */
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
        2
    }
}

/* enforces that both  "real" and "imaginary" components of two Fp2Gadgets are not equal.
Note: This is not the canonical notion of inequality for field elements, we need to check
whether the implementation depicts the need of futher application
*/
impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> NEqGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
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
        2
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> CondSelectGadget<ConstraintF>
for Fp2Gadget<P, ConstraintF>
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c0"),
            cond,
            &first.c0,
            &second.c0,
        )?;
        let c1 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c1"),
            cond,
            &first.c1,
            &second.c1,
        )?;

        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> TwoBitLookupGadget<ConstraintF>
for Fp2Gadget<P, ConstraintF>
{
    type TableConstant = Fp2<P>;
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
        let c0 = FpGadget::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = FpGadget::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
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
        let c0 = FpGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c0"), precomp, b, &c0s)?;
        let c1 = FpGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c1"), precomp, b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <FpGadget<ConstraintF> as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
ThreeBitCondNegLookupGadget<ConstraintF> for Fp2Gadget<P, ConstraintF>
{
    type TableConstant = Fp2<P>;

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
        let c0 = FpGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c0"), b, b0b1, &c0s)?;
        let c1 = FpGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c1"), b, b0b1, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <FpGadget<ConstraintF> as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

/*
Packing and unpacking gadgets for Fp2Gadget
*/

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> ToBitsGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits(&mut cs)?;
        let mut c1 = self.c1.to_bits(cs)?;
        c0.append(&mut c1);
        Ok(c0)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits_strict(&mut cs)?;
        let mut c1 = self.c1.to_bits_strict(cs)?;
        c0.append(&mut c1);
        Ok(c0)
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> ToBytesGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
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
