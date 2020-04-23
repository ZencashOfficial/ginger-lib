use algebra::{
    Field, bytes::ToBytes
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

pub mod gm17;
pub mod groth16;


pub trait NIZK {
    type Circuit;
    type AssignedCircuit;
    type VerifierInput: ?Sized;
    type ProvingParameters: Clone;
    type VerificationParameters: Clone + Default;
    type PreparedVerificationParameters: Clone + Default + From<Self::VerificationParameters>;
    type Proof: ToBytes + Clone + Default;
}

pub trait NIZKVerifierGadget<N: NIZK, ConstraintF: Field> {
    type VerificationKeyGadget: AllocGadget<N::VerificationParameters, ConstraintF> +
                                ToBytesGadget<ConstraintF>;
    type PreparedVerificationKeyGadget: FromGadget<Self::VerificationKeyGadget, ConstraintF>;
    type ProofGadget: AllocGadget<N::Proof, ConstraintF>;

    fn check_verify<'a, CS, I, T>(
        cs: CS,
        pvk: &Self::PreparedVerificationKeyGadget,
        input: I,
        proof: &Self::ProofGadget,
    ) -> Result<(), SynthesisError>
        where
            CS: ConstraintSystem<ConstraintF>,
            I: Iterator<Item = &'a T>,
            T: 'a + ToBitsGadget<ConstraintF> + ?Sized;
}

pub trait NIZKVerifierConstantGadget<N: NIZK, ConstraintF: Field, NG: NIZKVerifierGadget<N, ConstraintF>> {

    type PreparedVerificationKeyGadget:         FromGadget<NG::VerificationKeyGadget, ConstraintF>;
    type PreparedVerificationKeyConstantGadget: ConstantGadget<N::PreparedVerificationParameters, ConstraintF> +
                                                CondSelectGadget<ConstraintF> +
                                                Into<NG::PreparedVerificationKeyGadget>;

    fn check_verify_with_constant_pvk<'a, CS, I, T>(
        cs: CS,
        pvk: &'a N::PreparedVerificationParameters,
        input: I,
        proof: &NG::ProofGadget,
    ) -> Result<(), SynthesisError>
        where
            CS: ConstraintSystem<ConstraintF>,
            I: Iterator<Item = &'a T>,
            T: 'a + ToBitsGadget<ConstraintF> + ?Sized;

    /// If `cond` is true enforce verification of `proof` with `pvk_1`, otherwise
    /// with `pvk_2`.
    fn check_verify_with_constant_pvk_selection<'a, CS, I, T>(
        mut cs: CS,
        pvk_1: &'a N::PreparedVerificationParameters,
        pvk_2: &'a N::PreparedVerificationParameters,
        cond:  &Boolean,
        input: I,
        proof: &NG::ProofGadget,
    ) -> Result<(), SynthesisError>
        where
            CS: ConstraintSystem<ConstraintF>,
            I: Iterator<Item = &'a T>,
            T: 'a + ToBitsGadget<ConstraintF> + ?Sized
    {
        let constant_pvk_1 = Self::PreparedVerificationKeyConstantGadget::from_value(
            cs.ns(|| "Hardcode pvk1"),
            pvk_1
        );

        let constant_pvk_2 = Self::PreparedVerificationKeyConstantGadget::from_value(
            cs.ns(|| "Hardcode pvk2"),
            pvk_2
        );

        let selected_constant_pvk = Self::PreparedVerificationKeyConstantGadget::conditionally_select(
            cs.ns(|| "select pvk"),
            cond,
            &constant_pvk_1,
            &constant_pvk_2,
        )?;

        let selected_pvk =
            Self::PreparedVerificationKeyConstantGadget::into(selected_constant_pvk);

        NG::check_verify(
            cs.ns(|| "verify"),
            &selected_pvk,
            input,
            proof
        )?;

        Ok(())
    }
}