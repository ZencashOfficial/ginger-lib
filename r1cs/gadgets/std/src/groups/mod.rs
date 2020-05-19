/*
GroupGadget:
    - interfaces for group arithmetics and
    - default implementations for variable and fixed base exponentiation, as well as
    - the interface for 3-bit lookup table fixed base exponentiation.
There is also a default implementation for multi-scalar multi-base exponentiation as used by the
Pedersen CRH/commitment scheme -> better move it to there (see comments below).
Can be improved by providing a generic implementation for 3-bit (signed) lookup table exponentiation.
*/

use crate::prelude::*;
use algebra::{Field, Group};
use r1cs_core::{ConstraintSystem, SynthesisError};

use std::{borrow::Borrow, fmt::Debug};

pub mod curves;

pub use self::curves::{
    short_weierstrass::bls12,
    twisted_edwards::{edwards_sw6, jubjub},
};

pub trait GroupGadget<G: Group, ConstraintF: Field>:
    Sized
    + ToBytesGadget<ConstraintF>
    + NEqGadget<ConstraintF>
    + EqGadget<ConstraintF>
    + ToBitsGadget<ConstraintF>
    + CondSelectGadget<ConstraintF>
    + AllocGadget<G, ConstraintF>
    + ConstantGadget<G, ConstraintF>
    + Clone
    + Debug
{
    type Value: Debug;
    type Variable;

    fn get_value(&self) -> Option<Self::Value>;

    fn get_variable(&self) -> Self::Variable;

    fn zero<CS: ConstraintSystem<ConstraintF>>(cs: CS) -> Result<Self, SynthesisError>;

    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError>;

    fn sub<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let neg_other = other.negate(cs.ns(|| "Negate other"))?;
        self.add(cs.ns(|| "Self - other"), &neg_other)
    }

    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &G,
    ) -> Result<Self, SynthesisError>;

    fn sub_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &G,
    ) -> Result<Self, SynthesisError> {
        let neg_other = -(*other);
        self.add_constant(cs.ns(|| "Self - other"), &neg_other)
    }

    fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
    ) -> Result<(), SynthesisError>;

    fn negate<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError>;


    /* Exponentiation of an element P via linear combination of 2^i-th powers of P,
    given the scalar as little endian vector of Booleans b_0, b_1, ..., b_{l-1}:
        result = b_0 * P + b_1 * (2*P) + b_2 * (4*P) + ... + b_{l-1} * (2^{2^l}*P)
    Uses add, which in this generic implementation is assumed to be complete addition.
    WARNING: If add is incomplete then one has to have control over exceptional cases!
    */
    fn mul_bits<'a, CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        result: &Self,
        bits: impl Iterator<Item = &'a Boolean>,
    ) -> Result<Self, SynthesisError> {
        let mut power = self.clone();
        let mut result = result.clone();
        for (i, bit) in bits.enumerate() {
            // as this circuit is generic, one has to compute add in every single loop anyway,
            // even if not taken into account
            let new_encoded = result.add(&mut cs.ns(|| format!("Add {}-th power", i)), &power)?;
            result = Self::conditionally_select(
                &mut cs.ns(|| format!("Select {}", i)),
                bit.borrow(),
                &new_encoded,
                &result,
            )?;
            power.double_in_place(&mut cs.ns(|| format!("{}-th Doubling", i)))?;
        }
        Ok(result)
    }

    /* Fixed base exponentiation via linear combination of precomputed 2^i-th powers of B,
    given a scalar as little endian vector of Booleans b_0, b_1, ..., b_{l-1}:
        result = b_0 * B + b_1 * (2*B) + b_2 * (4*B) + ... + b_{l-1} * (2^{2^l}*B)
    Uses add, which in this generic implementation is assumed to be complete addition.
    WARNING: If add is incomplete then one has to have control over exceptional cases!
    */
    fn precomputed_base_scalar_mul<'a, CS, I, B>(
        &mut self,
        mut cs: CS,
        scalar_bits_with_base_powers: I,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Iterator<Item = (B, &'a G)>,
        B: Borrow<Boolean>,
        G: 'a,
    {
        for (i, (bit, base_power)) in scalar_bits_with_base_powers.enumerate() {
            let new_encoded = self.add_constant(
                &mut cs.ns(|| format!("Add {}-th base power", i)),
                &base_power,
            )?;
            *self = Self::conditionally_select(
                &mut cs.ns(|| format!("Conditional Select {}", i)),
                bit.borrow(),
                &new_encoded,
                &self,
            )?;
        }
        Ok(())
    }

    // what is this good for? does not use any pre-computations, just ordinary
    // variable-base scalar multiplication.
    fn mul_bits_fixed_base<'a, CS: ConstraintSystem<ConstraintF>>(
        base: &'a G,
        mut cs: CS,
        result: &Self,
        bits: &[Boolean],
    ) -> Result<Self, SynthesisError> {
        let base_g = Self::from_value(cs.ns(|| "hardcode base"), base);
        base_g.mul_bits(cs, result, bits.into_iter())
    }

    // why no default implementation for 3 bit lookup tables?
    fn precomputed_base_3_bit_signed_digit_scalar_mul<'a, CS, I, J, B>(
        _: CS,
        _: &[B],
        _: &[J],
    ) -> Result<Self, SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Borrow<[Boolean]>,
        J: Borrow<[I]>,
        B: Borrow<[G]>,
    {
        Err(SynthesisError::AssignmentMissing)
    }

    /* Multi-scalar multi-base exponentiation
        result = sum_j scalar_j * B_j
    using precomputed 2^{2^i}-th of the base(s).
    Why here? This type of exponentiation is specific to the Pedersen CRH/Commitment scheme.
    Although it can be used for single-base multi-scalar exponentiation, it does not provide
    the natural interface for it.

    WARNING: in many applications 'a::to_bits need to be secure unpacking.
    */
    fn precomputed_base_multiscalar_mul<'a, CS, T, I, B>(
        mut cs: CS,
        bases: &[B], // a vector of a vector of the base multiples of the bases
        scalars: I, // typically a Vec of Vec of bits
    ) -> Result<Self, SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        T: 'a + ToBitsGadget<ConstraintF> + ?Sized,
        I: Iterator<Item = &'a T>,
        B: Borrow<[G]>,
    {
        let mut result = Self::zero(&mut cs.ns(|| "Declare Result"))?;
        // Compute ∏(h_i^{m_i}) for all i.
        for (i, (bits, base_powers)) in scalars.zip(bases).enumerate() {
            let base_powers = base_powers.borrow();
            let bits = bits.to_bits(&mut cs.ns(|| format!("Convert Scalar {} to bits", i)))?;
            result.precomputed_base_scalar_mul(
                cs.ns(|| format!("Window {}", i)),
                bits.iter().zip(base_powers),
            )?;
        }
        Ok(result)
    }

    fn cost_of_add() -> usize;

    fn cost_of_double() -> usize;
}

#[cfg(test)]
mod test {
    use algebra::{Field, ProjectiveCurve, ToCompressedBits};
    use r1cs_core::ConstraintSystem;

    use crate::{prelude::*, test_constraint_system::TestConstraintSystem, ToCompressedBitsGadget};
    use algebra::groups::Group;
    use rand;
    use crate::groups::curves::short_weierstrass::short_weierstrass_projective::CompressAffinePointGadget;

    pub(crate) fn group_test<
        ConstraintF: Field,
        G: Group,
        GG: GroupGadget<G, ConstraintF>,
        CS: ConstraintSystem<ConstraintF>,
    >(
        cs: &mut CS,
        a: GG,
        b: GG,
    ) {
        let zero = GG::zero(cs.ns(|| "Zero")).unwrap();
        assert_eq!(zero, zero);

        // a == a
        assert_eq!(a, a);

        // a + 0 = a
        assert_eq!(a.add(cs.ns(|| "a_plus_zero"), &zero).unwrap(), a);
        // a - 0 = a
        assert_eq!(a.sub(cs.ns(|| "a_minus_zero"), &zero).unwrap(), a);
        // a - a = 0
        assert_eq!(a.sub(cs.ns(|| "a_minus_a"), &a).unwrap(), zero);

        // a + b = b + a
        let a_b = a.add(cs.ns(|| "a_plus_b"), &b).unwrap();
        let b_a = b.add(cs.ns(|| "b_plus_a"), &a).unwrap();
        assert_eq!(a_b, b_a);
        // (a + b) + a = a + (b + a)
        let ab_a = a_b.add(&mut cs.ns(|| "a_b_plus_a"), &a).unwrap();
        let a_ba = a.add(&mut cs.ns(|| "a_plus_b_a"), &b_a).unwrap();
        assert_eq!(ab_a, a_ba);

        // a.double() = a + a
        let a_a = a.add(cs.ns(|| "a + a"), &a).unwrap();
        let mut a2 = a.clone();
        a2.double_in_place(cs.ns(|| "2a")).unwrap();
        assert_eq!(a2, a_a);
        // b.double() = b + b
        let mut b2 = b.clone();
        b2.double_in_place(cs.ns(|| "2b")).unwrap();
        let b_b = b.add(cs.ns(|| "b + b"), &b).unwrap();
        assert_eq!(b2, b_b);

        let _ = a.to_bytes(&mut cs.ns(|| "ToBytes")).unwrap();
        let _ = a.to_bytes_strict(&mut cs.ns(|| "ToBytes Strict")).unwrap();

        let _ = b.to_bytes(&mut cs.ns(|| "b ToBytes")).unwrap();
        let _ = b
            .to_bytes_strict(&mut cs.ns(|| "b ToBytes Strict"))
            .unwrap();
    }

    pub(crate) fn group_test_with_unsafe_add<
        ConstraintF: Field,
        G: Group,
        GG: GroupGadget<G, ConstraintF>,
        CS: ConstraintSystem<ConstraintF>,
    >(
        cs: &mut CS,
        a: GG,
        b: GG,
    ) {
        let zero = GG::zero(cs.ns(|| "Zero")).unwrap();
        assert_eq!(zero, zero);

        // a == a
        assert_eq!(a, a);

        // a + b = b + a
        let a_b = a.add(cs.ns(|| "a_plus_b"), &b).unwrap();
        let b_a = b.add(cs.ns(|| "b_plus_a"), &a).unwrap();
        assert_eq!(a_b, b_a);

        // (a + b) + a = a + (b + a)
        let ab_a = a_b.add(&mut cs.ns(|| "a_b_plus_a"), &a).unwrap();
        let a_ba = a.add(&mut cs.ns(|| "a_plus_b_a"), &b_a).unwrap();
        assert_eq!(ab_a, a_ba);

        // a.double() + b = (a + b) + a: Testing double() using b as shift
        let mut a2 = a.clone();
        a2.double_in_place(cs.ns(|| "2a")).unwrap();
        let a2_b = a2.add(cs.ns(|| "2a + b"), &b).unwrap();

        let a_b_a = a.add(cs.ns(|| "a + b"), &b).unwrap()
            .add(cs.ns(|| "a + b + a"), &a).unwrap();
        assert_eq!(a2_b, a_b_a);

        // (b.double() + a) = (b + a) + b: Testing double() using a as shift
        let mut b2 = b.clone();
        b2.double_in_place(cs.ns(|| "2b")).unwrap();
        let b2_a = b2.add(cs.ns(|| "2b + a"), &a).unwrap();

        let b_a_b = b.add(cs.ns(|| "b + a"), &a).unwrap()
            .add(cs.ns(|| "b + a + b"), &b).unwrap();
        assert_eq!(b2_a, b_a_b);

        let _ = a.to_bytes(&mut cs.ns(|| "ToBytes")).unwrap();
        let _ = a.to_bytes_strict(&mut cs.ns(|| "ToBytes Strict")).unwrap();

        let _ = b.to_bytes(&mut cs.ns(|| "b ToBytes")).unwrap();
        let _ = b
            .to_bytes_strict(&mut cs.ns(|| "b ToBytes Strict"))
            .unwrap();
    }

    #[test]
    fn jubjub_group_gadgets_test() {
        use crate::groups::jubjub::JubJubGadget;
        use algebra::{curves::jubjub::JubJubProjective, fields::jubjub::fq::Fq};

        let mut cs = TestConstraintSystem::<Fq>::new();

        let a: JubJubProjective = rand::random();
        let b: JubJubProjective = rand::random();

        let a = JubJubGadget::alloc(&mut cs.ns(|| "generate_a"), || Ok(a)).unwrap();
        let b = JubJubGadget::alloc(&mut cs.ns(|| "generate_b"), || Ok(b)).unwrap();
        group_test::<_, JubJubProjective, _, _>(&mut cs.ns(|| "GroupTest(a, b)"), a, b);
    }

    #[test]
    fn mnt4_group_gadgets_test() {
        use crate::groups::curves::short_weierstrass::mnt::mnt4::mnt4753::MNT4G1Gadget;
        use algebra::{
            curves::mnt4753::G1Projective as MNT4G1Projective,
            fields::mnt4753::Fq,
        };

        let mut cs = TestConstraintSystem::<Fq>::new();

        //Test G1
        let a: MNT4G1Projective = rand::random();
        let b: MNT4G1Projective = rand::random();

        let a = MNT4G1Gadget::alloc(&mut cs.ns(|| "generate_a_g1"), || Ok(a)).unwrap();
        let b = MNT4G1Gadget::alloc(&mut cs.ns(|| "generate_b_g1"), || Ok(b)).unwrap();
        group_test_with_unsafe_add::<_, MNT4G1Projective, _, _>(&mut cs.ns(|| "GroupTest(a, b)_g1"), a, b);

        let p1: MNT4G1Projective = rand::random();
        let p1_compressed = p1.into_affine().compress();
        let p1_gadget = MNT4G1Gadget::alloc(&mut cs.ns(|| "generate_p1"), || Ok(p1)).unwrap();
        let p1_compression_gadget = CompressAffinePointGadget::<Fq>::new(p1_gadget.x, p1_gadget.y, p1_gadget.infinity);
        let p1_compressed_by_gadget = p1_compression_gadget.to_compressed(cs.ns(|| "p1 compressed g1")).unwrap();
        let mut p1_compressed_by_gadget_conv = vec![];
        for b in p1_compressed_by_gadget {
            p1_compressed_by_gadget_conv.push(b.get_value().unwrap());
        }
        assert_eq!(p1_compressed_by_gadget_conv, p1_compressed);

    }

    #[test]
    fn mnt6_group_gadgets_test() {
        use crate::groups::curves::short_weierstrass::mnt::mnt6::mnt6753::MNT6G1Gadget;
        use algebra::{
            curves::mnt6753::G1Projective as MNT6G1Projective,
            fields::mnt6753::Fq,
        };

        let mut cs = TestConstraintSystem::<Fq>::new();

        //Test G1
        let a: MNT6G1Projective = rand::random();
        let b: MNT6G1Projective = rand::random();

        let a = MNT6G1Gadget::alloc(&mut cs.ns(|| "generate_a_g1"), || Ok(a)).unwrap();
        let b = MNT6G1Gadget::alloc(&mut cs.ns(|| "generate_b_g1"), || Ok(b)).unwrap();
        group_test_with_unsafe_add::<_, MNT6G1Projective, _, _>(&mut cs.ns(|| "GroupTest(a, b)_g1"), a, b);

        let p1: MNT6G1Projective = rand::random();
        let p1_compressed = p1.into_affine().compress();
        let p1_gadget = MNT6G1Gadget::alloc(&mut cs.ns(|| "generate_p1"), || Ok(p1)).unwrap();
        let p1_compression_gadget = CompressAffinePointGadget::<Fq>::new(p1_gadget.x, p1_gadget.y, p1_gadget.infinity);
        let p1_compressed_by_gadget = p1_compression_gadget.to_compressed(cs.ns(|| "p1 compressed g1")).unwrap();
        let mut p1_compressed_by_gadget_conv = vec![];
        for b in p1_compressed_by_gadget {
            p1_compressed_by_gadget_conv.push(b.get_value().unwrap());
        }
        assert_eq!(p1_compressed_by_gadget_conv, p1_compressed);

    }
}