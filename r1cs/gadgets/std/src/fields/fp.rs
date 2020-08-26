/*
Definition of the FpGadget and implementation of the following gadgets for it:
    - FieldGadget,
    - AllocGadget, CloneGadget, ConstantGadget,
    - PartialEqGadget, ConditionalEqGadget, NEqGadget,
    - CondSelectGadget, TwoBitLookupGadget, ThreeBitNegLookupGadget,
    - ToBitsGadget, FromBitsGadget, ToBytesGadget
*/

use algebra::{bytes::ToBytes, FpParameters, PrimeField};
use r1cs_core::{
    ConstraintSystem,
    ConstraintVar::{self, *},
    LinearCombination, SynthesisError,
};

use std::borrow::Borrow;

use crate::{boolean::AllocatedBit, prelude::*, Assignment};

#[derive(Debug)]
pub struct FpGadget<F: PrimeField> {
    pub value:    Option<F>,
    pub variable: ConstraintVar<F>,
}

// extra functions for the FpGadget
impl<F: PrimeField> FpGadget<F> {

    #[inline]
    pub fn from<CS: ConstraintSystem<F>>(mut cs: CS, value: &F) -> Self {
        Self::alloc(cs.ns(|| "from"), || Ok(*value)).unwrap()
    }

    /* the odd test for FpGadgets, relys on a uses a rather expensive comparison in value by
    using to_bits_strict
    */
    #[inline]
    pub fn is_odd<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
    ) -> Result<Boolean, SynthesisError> {
        let bits = self.to_bits_strict(cs.ns(|| "to bits strict"))?;
        Ok(bits[bits.len() - 1])
    }

    /* An inexpensive alternative for secure unpacking of FpGadgets of bit length strictly
     smaller than the bit length of the modulus.

     Assumes that skip_leading_bits > 0.

     DANGER: do not use with skip_leading_bits = 0. In this case, unpacking is not secure!
     Use to_bits_strict instead.
     */
    #[inline]
    pub fn to_bits_with_length_restriction<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        skip_leading_bits: usize,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        // in order to avoid skip_leading_bits = 0 we could select num_bits as F::Params::CAPACITY,
        // and renaming skip_leading_bits into skip_further_bits or skip_further_leading_bits
        let num_bits = F::Params::MODULUS_BITS;

        let bit_values = match self.value {
            Some(value) => {
                value.write_bits().iter().map(|b| Some(*b)).collect::<Vec<_>>()
            },
            None => vec![None; num_bits as usize],
        };

        let mut bits = vec![];
        for (i, b) in bit_values.iter().skip(skip_leading_bits).enumerate() {
            bits.push(AllocatedBit::alloc(cs.ns(|| format!("bit {}", i)), || {
                b.get()
            })?);
        }

        let mut lc = LinearCombination::zero();
        let mut coeff = F::one();

        // construct sum_{i=0}^{bitlen} b_i * 2^i as linear combination of the powers of 2.
        for bit in bits.iter().rev() {
            lc = lc + (coeff, bit.get_variable());
            coeff.double_in_place();
        }

        lc = &self.variable - lc;

        cs.enforce(|| "unpacking_constraint", |lc| lc, |lc| lc, |_| lc);

        Ok(bits.into_iter().map(Boolean::from).collect())
    }

    /* An inexpensive but secure alternative for unpacking FpGadgets of byte length strictly
    smaller than the byte length of the modulus.

    Assumes that to_skip > 0.

    DANGER: do not use with to_skip = 0. In this case, unpacking is not secure!
    Use to_bytes_strict instead.
    */
    #[inline]
    pub fn to_bytes_with_length_restriction<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        to_skip: usize
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut byte_values = match self.value {
            Some(value) => to_bytes![&value.into_repr()]?
                .into_iter()
                .map(Some)
                .collect::<Vec<_>>(),
            None => {
                let default = F::default();
                let default_len = to_bytes![&default].unwrap().len();
                vec![None; default_len]
            },
        };

        for _ in 0..to_skip {byte_values.pop();}

        let bytes = UInt8::alloc_vec(cs.ns(|| "Alloc bytes"), &byte_values)?;

        let mut lc = LinearCombination::zero();
        let mut coeff = F::one();

        for bit in bytes
            .iter()
            .flat_map(|byte_gadget| byte_gadget.bits.clone())
            {
                match bit {
                    Boolean::Is(bit) => {
                        lc = lc + (coeff, bit.get_variable());
                        coeff.double_in_place();
                    },
                    Boolean::Constant(_) | Boolean::Not(_) => unreachable!(),
                }
            }

        lc = &self.variable - lc;

        cs.enforce(|| "unpacking_constraint", |lc| lc, |lc| lc, |_| lc);

        Ok(bytes)
    }

}


/* FieldGadget implementation for the FpGadget
*/
impl<F: PrimeField> FieldGadget<F, F> for FpGadget<F> {
    type Variable = ConstraintVar<F>;

    #[inline]
    fn get_value(&self) -> Option<F> {
        self.value
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        self.variable.clone()
    }

    #[inline]
    fn zero<CS: ConstraintSystem<F>>(_cs: CS) -> Result<Self, SynthesisError> {
        let value = Some(F::zero());
        Ok(FpGadget {
            value,
            variable: ConstraintVar::zero(),
        })
    }

    #[inline]
    fn one<CS: ConstraintSystem<F>>(_cs: CS) -> Result<Self, SynthesisError> {
        let value = Some(F::one());
        Ok(FpGadget {
            value,
            variable: CS::one().into(),
        })
    }

    /*
    addition gadgets
    */

    #[inline]
    fn add<CS: ConstraintSystem<F>>(
        &self,
        mut _cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let value = match (self.value, other.value) {
            (Some(val1), Some(val2)) => Some(val1 + &val2),
            (..) => None,
        };

        Ok(FpGadget {
            value,
            variable: &self.variable + &other.variable,
        })
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<F>>(
        &self,
        _cs: CS,
        other: &F,
    ) -> Result<Self, SynthesisError> {
        let value = self.value.map(|val| val + other);
        Ok(FpGadget {
            value,
            variable: self.variable.clone() + (*other, CS::one()),
        })
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<F>>(
        &mut self,
        _cs: CS,
        other: &F,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| *val += other);
        self.variable += (*other, CS::one());
        Ok(self)
    }

    #[inline]
    fn double<CS: ConstraintSystem<F>>(&self, _cs: CS) -> Result<Self, SynthesisError> {
        let value = self.value.map(|val| val.double());
        let mut variable = self.variable.clone();
        variable.double_in_place();
        Ok(FpGadget { value, variable })
    }

    #[inline]
    fn double_in_place<CS: ConstraintSystem<F>>(
        &mut self,
        _cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| val.double_in_place());
        self.variable.double_in_place();
        Ok(self)
    }

    /* conditional adding of a constant
        result = self + bit * constant,
    where bit is Boolean.
    */
    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<F>>(
        &self,
        mut _cs: CS,
        bit: &Boolean,
        coeff: F,
    ) -> Result<Self, SynthesisError> {
        let value = match (self.value, bit.get_value()) {
            (Some(v), Some(b)) => Some(if b { v + &coeff } else { v }),
            (..) => None,
        };
        Ok(FpGadget {
            value,
            variable: LC(bit.lc(CS::one(), coeff)) + &self.variable,
        })
    }

    /*
    substraction gadgets
    */

    #[inline]
    fn sub<CS: ConstraintSystem<F>>(
        &self,
        mut _cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let value = match (self.value, other.value) {
            (Some(val1), Some(val2)) => Some(val1 - &val2),
            (..) => None,
        };

        Ok(FpGadget {
            value,
            variable: &self.variable - &other.variable,
        })
    }

    // sub_in_place not implemented

    #[inline]
    fn negate<CS: ConstraintSystem<F>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.negate_in_place(cs)?;
        Ok(result)
    }

    #[inline]
    fn negate_in_place<CS: ConstraintSystem<F>>(
        &mut self,
        _cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| *val = -(*val));
        self.variable.negate_in_place();
        Ok(self)
    }

    /*
    multiplication gadgets
    */

    #[inline]
    fn mul<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let product = Self::alloc(cs.ns(|| "mul"), || {
            Ok(self.value.get()? * &other.value.get()?)
        })?;
        cs.enforce(
            || "mul_constraint",
            |lc| &self.variable + lc,
            |lc| &other.variable + lc,
            |lc| &product.variable + lc,
        );
        Ok(product)
    }

    #[inline]
    fn mul_by_constant<CS: ConstraintSystem<F>>(
        &self,
        cs: CS,
        other: &F,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.mul_by_constant_in_place(cs, other)?;
        Ok(result)
    }

    #[inline]
    fn mul_by_constant_in_place<CS: ConstraintSystem<F>>(
        &mut self,
        mut _cs: CS,
        other: &F,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| *val *= other);
        self.variable *= *other;
        Ok(self)
    }

    #[inline]
    fn inverse<CS: ConstraintSystem<F>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
        let inverse = Self::alloc(cs.ns(|| "inverse"), || {
            let result = self.value.get()?;
            let inv = result.inverse().expect("Inverse doesn't exist!");
            Ok(inv)
        })?;

        let one = CS::one();
        cs.enforce(
            || "inv_constraint",
            |lc| &self.variable + lc,
            |lc| &inverse.variable + lc,
            |lc| lc + one,
        );
        Ok(inverse)
    }

    #[inline]
    fn mul_equals<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        cs.enforce(
            || "mul_constraint",
            |lc| &self.variable + lc,
            |lc| &other.variable + lc,
            |lc| &result.variable + lc,
        );
        Ok(())
    }

    #[inline]
    fn square_equals<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        cs.enforce(
            || "sqr_constraint",
            |lc| &self.variable + lc,
            |lc| &self.variable + lc,
            |lc| &result.variable + lc,
        );
        Ok(())
    }

    #[inline]
    fn frobenius_map<CS: ConstraintSystem<F>>(
        &self,
        _: CS,
        _: usize,
    ) -> Result<Self, SynthesisError> {
        // as Fp is the base field of its own, the Frobenius map is trivial
        Ok(self.clone())
    }

    #[inline]
    fn frobenius_map_in_place<CS: ConstraintSystem<F>>(
        &mut self,
        _: CS,
        _: usize,
    ) -> Result<&mut Self, SynthesisError> {
        Ok(self)
    }

    // constraint counts
    fn cost_of_mul() -> usize { 1 }

    fn cost_of_mul_equals() -> usize { 1 }

    fn cost_of_inv() -> usize { 1 }
}

/*
Alloc-, Clone and ConstantGadget for the FpGadget
*/

impl<F: PrimeField> AllocGadget<F, F> for FpGadget<F> {
    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<F>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<F>,
    {
        let mut value = None;
        let variable = cs.alloc(
            || "alloc",
            || {
                let tmp = *value_gen()?.borrow();
                value = Some(tmp);
                Ok(tmp)
            },
        )?;
        Ok(FpGadget {
            value,
            variable: Var(variable),
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<F>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<F>,
    {
        let mut value = None;
        let variable = cs.alloc_input(
            || "alloc",
            || {
                let tmp = *value_gen()?.borrow();
                value = Some(tmp);
                Ok(tmp)
            },
        )?;

        Ok(FpGadget {
            value,
            variable: Var(variable),
        })
    }
}

impl<F: PrimeField> Clone for FpGadget<F> {
    fn clone(&self) -> Self {
        Self {
            value:    self.value.clone(),
            variable: self.variable.clone(),
        }
    }
}

impl<F: PrimeField> ConstantGadget<F, F> for FpGadget<F> {

    #[inline]
    fn from_value<CS: ConstraintSystem<F>>(
        _cs: CS,
        value: &F,
    ) -> Self
    {
        let value = *value;
        FpGadget {
            value: Some(value),
            variable: ConstraintVar::<F>::from((value, CS::one())),
        }
    }

    #[inline]
    fn get_constant(&self) -> F {
        self.get_value().unwrap()
    }
}

/*
relational and conditional gadgets (incl. lookup tables) for the FpGadget
*/

impl<F: PrimeField> PartialEq for FpGadget<F> {
    fn eq(&self, other: &Self) -> bool {
        !self.value.is_none() && !other.value.is_none() && self.value == other.value
    }
}

impl<F: PrimeField> Eq for FpGadget<F> {}

impl<F: PrimeField> EqGadget<F> for FpGadget<F> {}

impl<F: PrimeField> ConditionalEqGadget<F> for FpGadget<F> {
    /* enforces equality of two FpGadgets x,y if the condition Boolean is true, otherwise it does
    not enforce anything.
    */
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        let difference = self.sub(cs.ns(|| "difference"), other)?;
        let one = CS::one();
        let one_const = F::one();
        // (self - other) * condition = 0
        cs.enforce(
            || "conditional_equals",
            |lc| &difference.variable + lc,
            |lc| lc + &condition.lc(one, one_const),
            |lc| lc,
        );
        Ok(())
    }

    fn cost() -> usize {
        1
    }
}

impl<F: PrimeField> NEqGadget<F> for FpGadget<F> {
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        let a_minus_b = self.sub(cs.ns(|| "A - B"), other)?;
        a_minus_b.inverse(cs.ns(|| "Enforce inverse exists"))?;
        Ok(())
    }

    fn cost() -> usize {
        1
    }
}

impl<F: PrimeField> EquVerdictGadget<F> for FpGadget<F> {
    /* outputs a Boolean verdict v which is true/false depending on whether the values x,y
    of two FpGadgets are equal or not.
    */
    fn enforce_verdict<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        // The Boolean 'verdict' we want to constrain.
        let v = Boolean::alloc(cs.ns(|| "alloc verdict"), || {
            let self_val = self.get_value().get()?;
            let other_val = other.get_value().get()?;
            Ok(self_val == other_val)
        })?;
        /* The constraints on v,x,y are
                0 = v * (x - y),
            1 - v = c * (x - y).
        If v=1 the first equation enforces x=y (otherwise it enforces nothing). The auxiliary
        variable c is arbitrary in this case.
        If v=0 the second equation enforces x!=y (and otherwise nothing). In this case, c is the
        mult. inverse of x-y.
        */
        let c = Self::alloc(cs.ns(|| "alloc c"), || {
            let v_val = v.get_value().get()?;
            if v_val {
               Ok(F::one()) //Just an arbitrary value
            }
            else {
                let self_val = self.get_value().get()?;
                let other_val = other.get_value().get()?;
                Ok((self_val - &other_val).inverse().get()?)
            }
        })?;

        self.conditional_enforce_equal(cs.ns(|| "0 = v * (x - y)"), other, &v)?;

        let one = CS::one();
        cs.enforce(
            || "1 - v = c * (x - y)",
            |lc| (&self.variable - &other.variable) + lc,
            |lc| &c.variable + lc,
            |_| v.not().lc(one, F::one()),
        );

        Ok(v)
    }
}

impl<F: PrimeField> CondSelectGadget<F> for FpGadget<F> {
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<F>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        if let Boolean::Constant(cond) = *cond {
            if cond {
                Ok(first.clone())
            } else {
                Ok(second.clone())
            }
        } else {
            let result = Self::alloc(cs.ns(|| ""), || {
                cond.get_value()
                    .and_then(|cond| if cond { first } else { second }.get_value())
                    .get()
            })?;
            // a = self; b = other; c = cond, r = result;
            //
            // r = c * a + (1  - c) * b, or as rank 1 constraint
            // c * (a - b) = r - b
            let one = CS::one();
            cs.enforce(
                || "conditionally_select",
                |_| cond.lc(one, F::one()),
                |lc| (&first.variable - &second.variable) + lc,
                |lc| (&result.variable - &second.variable) + lc,
            );

            Ok(result)
        }
    }

    fn cost() -> usize {
        1
    }
}

/* Two bit lookup circuit, chooses c[i], i=0,..,3, according to i = b[1]*2 + b[0] given
by two Booleans b[0], b[1].
*/
impl<F: PrimeField> TwoBitLookupGadget<F> for FpGadget<F> {
    type TableConstant = F;
    fn two_bit_lookup<CS: ConstraintSystem<F>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        debug_assert!(b.len() == 2);
        debug_assert!(c.len() == 4);

        let result = Self::alloc(cs.ns(|| "Allocate lookup result"), || {
            match (b[0].get_value().get()?, b[1].get_value().get()?) {
                (false, false) => Ok(c[0]),
                (false, true) => Ok(c[2]),
                (true, false) => Ok(c[1]),
                (true, true) => Ok(c[3]),
            }
        })?;
        /* explicit lookup formula for table c = c[0],..,c[3]
            result = (1-b1) * [(1-b0)*c0 + b0*c1] + b1 * [(1-b0)*c2 + b0*c3],
        or as single rank 1 constraint
           (c0-c1-c2+c3) * b0 * b1 = result - [c0 + (c1-c0)*b0 + (c2-c0)*b1]
        */
        let one = CS::one();
        cs.enforce(
            || "Enforce lookup",
            |lc| lc + b[1].lc(one, c[3] - &c[2] - &c[1] + &c[0]) + (c[1] - &c[0], one),
            |lc| lc + b[0].lc(one, F::one()),
            |lc| result.get_variable() + lc + (-c[0], one) + b[1].lc(one, c[0] - &c[2]),
        );

        Ok(result)
    }

    /* low level version of two bit lookup table, outsources the quadratic term b0 * b1 in
        (c0-c1-c2+c3) * b0 * b1 = result - [c0 + (c1-c0)*b0 + (c2-c0)*b1]
    to an input Boolean 'precomp'.
    */
    fn two_bit_lookup_lc<CS: ConstraintSystem<F>>
    (   mut cs: CS,
        precomp: &Boolean,
        b: &[Boolean],
        c: &[Self::TableConstant]
    ) -> Result<Self, SynthesisError> {

        let result = Self::zero(cs.ns(|| "alloc result"))?
            .conditionally_add_constant(cs.ns(|| "add constant"),
                                        &Boolean::constant(true),
                                        c[0])?
            .conditionally_add_constant(cs.ns(|| "add b0"),
                                        &b[0],
                                        c[1] - &c[0])?
            .conditionally_add_constant(cs.ns(|| "add b1"),
                                        &b[1],
                                        c[2] - &c[0])?
            .conditionally_add_constant(cs.ns(|| "add b0 AND b1"),
                                        &precomp,
                                        c[3] + &c[0] - &c[1] - &c[2])?;

        Ok(result)
    }

    fn cost() -> usize {
        1
    }
}

/* Three bit look up circuit combining a two bit table and conditional negation.
Outputs
    c = (-1)^b[2] * c[i] = (-1)^b[2] * c[b[1]*2 + b[0]]
according to the bits b[0],b[1],b[2].
As two_bit_lookup_lc, the quadratic term b[0]*[1] is outsourced to the input Boolean b0b1.
*/
impl<F: PrimeField> ThreeBitCondNegLookupGadget<F> for FpGadget<F> {
    type TableConstant = F;

    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<F>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        debug_assert!(b.len() == 3);
        debug_assert!(c.len() == 4);

        let result = Self::alloc(cs.ns(|| "Allocate lookup result"), || {
            let y = match (b[0].get_value().get()?, b[1].get_value().get()?) {
                (false, false) => c[0],
                (false, true) => c[2],
                (true, false) => c[1],
                (true, true) => c[3],
            };
            if b[2].get_value().get()? {
                Ok(-y)
            } else {
                Ok(y)
            }
        })?;

        let one = CS::one();
        // y_lc = [c0 + (c1-c0)*b0 + (c2-c0)*b1] + (c0-c1-c2+c3) * b0b1
        let y_lc = b0b1.lc(one, c[3] - &c[2] - &c[1] + &c[0])
            + b[0].lc(one, c[1] - &c[0])
            + b[1].lc(one, c[2] - &c[0])
            + (c[0], one);
        // result = (1-b[2]) * y_lc  - b[2] * y_lc,
        // or as single rank 1 constraint
        // 2*b[2]*y_lc = y_lc - result
        cs.enforce(
            || "Enforce lookup",
            |_| y_lc.clone() + y_lc.clone(),
            |lc| lc + b[2].lc(one, F::one()),
            |_| -result.get_variable() + y_lc.clone(),
        );

        Ok(result)
    }

    fn cost() -> usize {
        2 // shouldn't be that 1 as quadratic term b0b1 is outsourced?
    }
}

/*
Packing and unpacking gadgets for FpGadget
*/

impl<F: PrimeField> ToBitsGadget<F> for FpGadget<F> {
    /// Outputs the binary representation of the value in `self` in *big-endian*
    /// form.

    /* Insecure unpacking, potentially allows that the vector of Booleans is not unique (they either
    represent the field element or the field element plus the modulus).
    DANGER: only use this when you really have thought about it!
    */
    fn to_bits<CS: ConstraintSystem<F>>(&self, mut cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        self.to_bits_with_length_restriction(&mut cs, 0)
    }

    /* Secure unpacking, enforces the Booleans to be the integer bits of the field element < modulus
    (and not the field element plus the modulus) involving an expensive comparison by value.
    */
    fn to_bits_strict<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let bits = self.to_bits(&mut cs)?;
        Boolean::enforce_in_field::<_, _, F>(&mut cs, &bits)?;

        Ok(bits)
    }
}

impl<F: PrimeField> FromBitsGadget<F> for FpGadget<F> {
    /* Secure packing of a Boolean vector into an FpGadget, assumes that the length of the vector
    is strictly smaller than the length of the field modulus.

    DANGER: for Boolean vector length >= modulus length this circuit is NOT secure, i.e. it does not
    enforce that the field element has the input bits as its integer representation.
    */
    fn from_bits<CS: ConstraintSystem<F>>(mut cs: CS, bits: &[Boolean]) -> Result<Self, SynthesisError> {
        // bits is assumed to be in big endian order
        let bits = bits.chunks(F::Params::CAPACITY as usize).next().unwrap();

        let mut num = Self::zero(cs.ns(|| "alloc_lc_{}"))?;
        let mut coeff = F::one();

        // Need to reverse in order to reconstruct the field element, because we
        // assume having a *big_endian* bit representation of `Self`.
        for (j, bit) in bits.iter().rev().enumerate() {

            // Use a support FpGadget to hold the linear combination (needed because
            // the allocated bit won't have a value until proving time.
            num = num.conditionally_add_constant(
                cs.ns(|| format!("add_bit_{}", j)),
                bit,
                coeff,
            )?;

            coeff.double_in_place();
        }

        // Alloc the field gadget with the value resulting from bit linear combination
        let variable = Self::alloc(
            cs.ns(|| "variable"),
            || {
                let value = num.get_value().get()?;
                Ok(value)
            }
        )?;

        // num * 1 = variable
        cs.enforce(
            || "packing constraint",
            |lc| lc,
            |lc| lc,
            |lc| &variable.variable - &num.variable + lc,
        );
        Ok(variable)
    }
}

impl<F: PrimeField> ToBytesGadget<F> for FpGadget<F> {
    // the insecure unpacking uses to_bytes_with_length_restriction with to_skip = 0
    fn to_bytes<CS: ConstraintSystem<F>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes_with_length_restriction(&mut cs, 0)
    }

    // the secure unpacking function using the rather expensive Boolean::enforce_in_field
    // circuit.
    fn to_bytes_strict<CS: ConstraintSystem<F>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let bytes = self.to_bytes(&mut cs)?;
        Boolean::enforce_in_field::<_, _, F>(
            &mut cs,
            &bytes.iter()
                .flat_map(|byte_gadget| byte_gadget.into_bits_le())
                // This reverse maps the bits into big-endian form, as required by `enforce_in_field`.
                .rev()
                .collect::<Vec<_>>(),
        )?;

        Ok(bytes)
    }
}


