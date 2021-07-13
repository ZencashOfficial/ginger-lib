use algebra::PrimeField;
use primitives::crh::{
    SBox, SpongeMode, AlgebraicSponge,
    poseidon::{
        PoseidonHash, PoseidonSponge, PoseidonParameters
    }
};
use crate::crh::{
    SBoxGadget, FieldBasedHashGadget, AlgebraicSpongeGadget,
};
use r1cs_std::{
    fields::{
        FieldGadget, fp::FpGadget
    },
    alloc::ConstantGadget,
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::marker::PhantomData;

#[cfg(feature = "mnt4_753")]
pub mod mnt4753;
#[cfg(feature = "mnt4_753")]
pub use self::mnt4753::*;

#[cfg(feature = "mnt6_753")]
pub mod mnt6753;
#[cfg(feature = "mnt6_753")]
pub use self::mnt6753::*;

#[cfg(feature = "tweedle")]
pub mod tweedle;
#[cfg(feature = "tweedle")]
pub use self::tweedle::*;

#[cfg(feature = "bn_382")]
pub mod bn382;
#[cfg(feature = "bn_382")]
pub use self::bn382::*;


pub struct PoseidonHashGadget
<
    ConstraintF: PrimeField,
    P:           PoseidonParameters<Fr = ConstraintF>,
    SB:          SBox<Field = ConstraintF, Parameters = P>,
    SBG:         SBoxGadget<ConstraintF, SB>,
>
{
    _field:             PhantomData<ConstraintF>,
    _parameters:        PhantomData<P>,
    _sbox:              PhantomData<SB>,
    _sbox_gadget:       PhantomData<SBG>,
}

impl<
    ConstraintF: PrimeField,
    P:   PoseidonParameters<Fr = ConstraintF>,
    SB:  SBox<Field = ConstraintF, Parameters = P>,
    SBG: SBoxGadget<ConstraintF, SB>
> PoseidonHashGadget<ConstraintF, P, SB, SBG>
{

    fn poseidon_perm<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        state: &mut [FpGadget<ConstraintF>],
    ) -> Result<(), SynthesisError>
    {
        // index that goes over the round constants
        let mut round_cst_idx = 0;

        // First full rounds
        for i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                // Temporary workaround: hardcoding the round constant and using it
                // in the following add() constraint, instead of using add_constant(),
                // helps reducing the R1CS density a little.
                let rc = FpGadget::<ConstraintF>::from_value(
                    cs.ns(|| format!("hardcode round constant {}", round_cst_idx)),
                    &P::ROUND_CST[round_cst_idx]
                );
                *d = rc.add(cs.ns(|| format!("add_constant_{}", round_cst_idx)), d)?;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(||format!("S-Box_1_{}_{}",i, j)), d)?;
            }

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_first_full_round_{}", i)), state)?;

        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                // Temporary workaround: hardcoding the round constant and using it
                // in the following add() constraint, instead of using add_constant(),
                // helps reducing the R1CS density a little.
                let rc = FpGadget::<ConstraintF>::from_value(
                    cs.ns(|| format!("hardcode round constant {}", round_cst_idx)),
                    &P::ROUND_CST[round_cst_idx]
                );
                *d = rc.add(cs.ns(|| format!("add_constant_{}", round_cst_idx)), d)?;
                round_cst_idx += 1;
            }

            // Apply S-Box only to the first element of the state vector
            SBG::apply(
                cs.ns(||format!("S-Box_2_{}_{}",_i, 0)),
                &mut state[0]
            )?;

            // Perform the matrix mix
            Self::matrix_mix (cs.ns(|| format!("poseidon_mix_matrix_partial_round_{}", _i)), state)?;
        }

        // Second full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                // Temporary workaround: hardcoding the round constant and using it
                // in the following add() constraint, instead of using add_constant(),
                // helps reducing the R1CS density a little.
                let rc = FpGadget::<ConstraintF>::from_value(
                    cs.ns(|| format!("hardcode round constant {}", round_cst_idx)),
                    &P::ROUND_CST[round_cst_idx]
                );
                *d = rc.add(cs.ns(|| format!("add_constant_{}", round_cst_idx)), d)?;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            for (j, d) in state.iter_mut().enumerate() {
                SBG::apply(cs.ns(|| format!("S-Box_3_{}_{}", _i, j)), d)?;
            }

            // Perform the matrix mix
            Self::matrix_mix(cs.ns(|| format!("poseidon_mix_matrix_second_full_round_{}", _i)), state)?;
        }
        Ok(())
    }

    // Function that does the dot product for the mix matrix
    fn dot_prod<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        res: &mut FpGadget<ConstraintF>,
        state: &mut [FpGadget<ConstraintF>],
        mut start_idx_cst: usize,
    ) -> Result<(), SynthesisError>
    {
        for x in state.iter() {
            let elem = x.mul_by_constant(cs.ns(|| format!("partial_product_{}", start_idx_cst)), &P::MDS_CST[start_idx_cst])?;
            start_idx_cst += 1;
            (*res).add_in_place(cs.ns(|| format!("add_partial_product_{}", start_idx_cst)), &elem)?;
        }

        Ok(())
    }

    // Function that does the mix matrix
    fn matrix_mix<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        state: &mut [FpGadget<ConstraintF>],
    ) -> Result<(), SynthesisError>
    {

        // Check that the length of the state vector is t
        assert_eq!(state.len(), P::T);

        // Destination state vector
        let mut new_state = Vec::new();

        // Initialize new destination state vector with zero elements
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(cs.ns(|| format!("hardcode_new_state_elem_{}", i)), &P::ZERO);
            new_state.push(elem);
        }

        // Performs the dot products
        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::dot_prod(cs.ns(|| format!("poseidon_dot_product_{}", i)), &mut new_state[i], state, idx_cst)?;
            idx_cst += P::T;
        }

        // Copy result to the state vector
        for i in 0..P::T {
            state[i] = new_state[i].clone();
        }

        Ok(())
    }
}

impl<ConstraintF, P, SB, SBG> FieldBasedHashGadget<PoseidonHash<ConstraintF, P, SB>, ConstraintF>
    for PoseidonHashGadget<ConstraintF, P, SB, SBG>
        where
            ConstraintF: PrimeField,
            P:           PoseidonParameters<Fr = ConstraintF>,
            SB:          SBox<Field = ConstraintF, Parameters = P>,
            SBG:         SBoxGadget<ConstraintF, SB>,
{
    type DataGadget = FpGadget<ConstraintF>;

    fn enforce_hash_constant_length<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        input: &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>
    // Assumption:
    //     capacity c = 1
    {
        assert_ne!(input.len(), 0, "Input data array does not contain any data.");

        let mut state = Vec::new();
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(
                cs.ns(|| format!("hardcode_state_{}",i)),
                &P::AFTER_ZERO_PERM[i]
            );
            state.push(elem);
        }

        // calculate the number of cycles to process the input dividing in portions of rate elements
        let num_cycles = input.len() / P::R;
        // check if the input is a multiple of the rate by calculating the remainder of the division
        // the remainder of dividing the input length by the rate can be 1 or 0 because we are assuming
        // that the rate is 2
        let rem = input.len() % P::R;

        // index to process the input
        let mut input_idx = 0;
        // iterate of the portions of rate elements
        for i in 0..num_cycles {
            // add the elements to the state vector. Add rate elements
            for j in 0..P::R {
                state[j].add_in_place(cs.ns(|| format!("add_input_{}_{}", i, j)), &input[input_idx])?;
                input_idx += 1;
            }
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| format!("poseidon_perm_{}", i)), &mut state)?;
        }

        // in case the input is not a multiple of the rate, process the remainder part padding zeros
        if rem != 0 {
            for j in 0..rem {
                state[j].add_in_place(cs.ns(|| format!("poseidon_padding_add_{}",j)), &input[input_idx])?;
                input_idx += 1;
            }
            // apply permutation after adding the input vector
            Self::poseidon_perm(cs.ns(|| "poseidon_padding_perm"), &mut state)?;
        }

        // return the first element of the state vector as the hash digest
        Ok(state[0].clone())
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct PoseidonSpongeGadget
<
    ConstraintF: PrimeField,
    P:           PoseidonParameters<Fr = ConstraintF>,
    SB:          SBox<Field = ConstraintF, Parameters = P>,
    SBG:         SBoxGadget<ConstraintF, SB>,
>
{
    pub(crate) mode:    SpongeMode,
    pub(crate) state:   Vec<FpGadget<ConstraintF>>,
    pub(crate) pending: Vec<FpGadget<ConstraintF>>,
    _field:             PhantomData<ConstraintF>,
    _parameters:        PhantomData<P>,
    _sbox:              PhantomData<SB>,
    _sbox_gadget:       PhantomData<SBG>,
}

impl<ConstraintF, P, SB, SBG> PoseidonSpongeGadget<ConstraintF, P, SB, SBG>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>,
        SB:          SBox<Field = ConstraintF, Parameters = P>,
        SBG:         SBoxGadget<ConstraintF, SB>,
{
    fn enforce_permutation<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS
    ) -> Result<(), SynthesisError>
    {
        // add the elements to the state vector. Add rate elements
        for (i, (input, state)) in self.pending.iter().zip(self.state.iter_mut()).enumerate() {
            state.add_in_place(cs.ns(|| format!("add_input_{}_to_state", i)), input)?;
        }

        // apply permutation after adding the input vector
        PoseidonHashGadget::<ConstraintF, P, SB, SBG>::poseidon_perm(
            cs.ns(|| "poseidon_perm"),
            &mut self.state
        )?;

        self.pending.clear();

        Ok(())
    }

    fn enforce_update<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
        input: FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError>
    {
        self.pending.push(input);
        if self.pending.len() == P::R {
            self.enforce_permutation(cs)?;
        }
        Ok(())
    }
}

impl<ConstraintF, P, SB, SBG> AlgebraicSpongeGadget<PoseidonSponge<ConstraintF, P, SB>, ConstraintF>
for PoseidonSpongeGadget<ConstraintF, P, SB, SBG>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>,
        SB:          SBox<Field = ConstraintF, Parameters = P>,
        SBG:         SBoxGadget<ConstraintF, SB>,
{
    type DataGadget = FpGadget<ConstraintF>;

    fn new<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let mut state = Vec::with_capacity(P::T);
        for i in 0..P::T {
            let elem = FpGadget::<ConstraintF>::from_value(
                cs.ns(|| format!("hardcode_state_{}",i)),
                &P::AFTER_ZERO_PERM[i]
            );
            state.push(elem);
        }

        Ok(Self {
            mode: SpongeMode::Absorbing,
            state,
            pending: Vec::with_capacity(P::R),
            _field: PhantomData,
            _parameters: PhantomData,
            _sbox: PhantomData,
            _sbox_gadget: PhantomData
        })
    }

    fn get_state(&self) -> &[FpGadget<ConstraintF>] {
        &self.state
    }

    fn set_state(&mut self, state: Vec<FpGadget<ConstraintF>>) {
        assert_eq!(state.len(), P::T);
        self.state = state;
    }

    fn get_mode(&self) -> &SpongeMode {
        &self.mode
    }

    fn set_mode(&mut self, mode: SpongeMode) {
        self.mode = mode;
    }

    fn enforce_absorb<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        elems: &[Self::DataGadget]
    ) -> Result<(), SynthesisError> {
        if elems.len() > 0 {
            match self.mode {

                SpongeMode::Absorbing => {
                    elems.iter().enumerate().map(|(i, f)| {
                        self.enforce_update(cs.ns(|| format!("update_{}", i)), f.clone())
                    }).collect::<Result<(), SynthesisError>>()?;
                },

                SpongeMode::Squeezing => {
                    self.mode = SpongeMode::Absorbing;
                    self.enforce_absorb(cs, elems)?;
                }
            }
        }
        Ok(())
    }

    fn enforce_squeeze<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        num: usize
    ) -> Result<Vec<Self::DataGadget>, SynthesisError> {
        let mut outputs = Vec::with_capacity(num);

        if num > 0 {
            match self.mode {
                SpongeMode::Absorbing => {

                    if self.pending.len() == 0 {
                        outputs.push(self.state[0].clone());
                    } else {
                        self.enforce_permutation(
                            cs.ns(|| "permutation")
                        )?;

                        outputs.push(self.state[0].clone());
                    }
                    self.mode = SpongeMode::Squeezing;
                    outputs.append(&mut self.enforce_squeeze(
                        cs.ns(|| "squeeze remaining elements"),
                        num - 1
                    )?);
                },

                // If we were squeezing, then squeeze the required number of field elements
                SpongeMode::Squeezing => {
                    for i in 0..num {
                        debug_assert!(self.pending.len() == 0);

                        PoseidonHashGadget::<ConstraintF, P, SB, SBG>::poseidon_perm(
                            cs.ns(|| format!("poseidon_perm_{}", i)),
                            &mut self.state
                        )?;
                        outputs.push(self.state[0].clone());
                    }
                }
            }
        }
        Ok(outputs)
    }
}

impl<ConstraintF, P, SB, SBG> ConstantGadget<PoseidonSponge<ConstraintF, P, SB>, ConstraintF>
for PoseidonSpongeGadget<ConstraintF, P, SB, SBG>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>,
        SB:          SBox<Field = ConstraintF, Parameters = P>,
        SBG:         SBoxGadget<ConstraintF, SB>,
{
    fn from_value<CS: ConstraintSystem<ConstraintF>>(mut cs: CS, value: &PoseidonSponge<ConstraintF, P, SB>) -> Self {
        let state_g = Vec::<FpGadget<ConstraintF>>::from_value(
            cs.ns(|| "hardcode state"),
            &value.get_state().to_vec()
        );

        let pending_g = Vec::<FpGadget<ConstraintF>>::from_value(
            cs.ns(|| "hardcode pending"),
            &value.get_pending().to_vec()
        );

        Self {
            mode: value.get_mode().clone(),
            state: state_g,
            pending: pending_g,
            _field: PhantomData,
            _parameters: PhantomData,
            _sbox: PhantomData,
            _sbox_gadget: PhantomData
        }
    }

    fn get_constant(&self) -> PoseidonSponge<ConstraintF, P, SB> {
        PoseidonSponge::<ConstraintF, P, SB>::new(
            self.mode.clone(),
            self.state.get_constant(),
            self.pending.get_constant(),
        )
    }
}

impl<ConstraintF, P, SB, SBG> From<Vec<FpGadget<ConstraintF>>>
for PoseidonSpongeGadget<ConstraintF, P, SB, SBG>
    where
        ConstraintF: PrimeField,
        P:           PoseidonParameters<Fr = ConstraintF>,
        SB:          SBox<Field = ConstraintF, Parameters = P>,
        SBG:         SBoxGadget<ConstraintF, SB>,
{
    fn from(other: Vec<FpGadget<ConstraintF>>) -> Self {
        assert_eq!(other.len(), P::T);
        Self {
            mode: SpongeMode::Absorbing,
            state: other,
            pending: Vec::with_capacity(P::R),
            _field: PhantomData,
            _parameters: PhantomData,
            _sbox: PhantomData,
            _sbox_gadget: PhantomData
        }
    }
}

#[cfg(test)]
mod test {

    use algebra::PrimeField;
    use crate::crh::test::{constant_length_field_based_hash_gadget_native_test, algebraic_sponge_gadget_native_test};

    fn generate_inputs<F: PrimeField>(num: usize) -> Vec<F>{
        let mut inputs = Vec::with_capacity(num);
        for i in 1..=num {
            let input = F::from(i as u32);
            inputs.push(input);
        }
        inputs
    }

    #[cfg(feature = "mnt4_753")]
    #[test]
    fn poseidon_mnt4_753_gadget_native_test() {
        use crate::mnt4753::*;

        for ins in 1..=5 {
            constant_length_field_based_hash_gadget_native_test::<_, _, MNT4PoseidonHashGadget>(generate_inputs(ins));
            algebraic_sponge_gadget_native_test::<_, _, MNT4PoseidonSpongeGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "mnt6_753")]
    #[test]
    fn poseidon_mnt6_753_gadget_native_test() {
        use crate::mnt6753::*;

        for ins in 1..=5 {
            constant_length_field_based_hash_gadget_native_test::<_, _, MNT6PoseidonHashGadget>(generate_inputs(ins));
            algebraic_sponge_gadget_native_test::<_, _, MNT6PoseidonSpongeGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn crh_bn382_fr_primitive_gadget_test() {
        use crate::bn382::*;

        for ins in 1..=5 {
            constant_length_field_based_hash_gadget_native_test::<_, _, BN382FrPoseidonHashGadget>(generate_inputs(ins));
            algebraic_sponge_gadget_native_test::<_, _, BN382FrPoseidonSpongeGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "bn_382")]
    #[test]
    fn crh_bn382_fq_primitive_gadget_test() {
        use crate::bn382::*;

        for ins in 1..=5 {
            constant_length_field_based_hash_gadget_native_test::<_, _, BN382FqPoseidonHashGadget>(generate_inputs(ins));
            algebraic_sponge_gadget_native_test::<_, _, BN382FqPoseidonSpongeGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn crh_tweedle_fr_primitive_gadget_test() {
        use crate::tweedle::*;

        for ins in 1..=5 {
            constant_length_field_based_hash_gadget_native_test::<_, _, TweedleFrPoseidonHashGadget>(generate_inputs(ins));
            algebraic_sponge_gadget_native_test::<_, _, TweedleFrPoseidonSpongeGadget>(generate_inputs(ins));
        }
    }

    #[cfg(feature = "tweedle")]
    #[test]
    fn crh_tweedle_fq_primitive_gadget_test() {
        use crate::tweedle::*;

        for ins in 1..=5 {
            constant_length_field_based_hash_gadget_native_test::<_, _, TweedleFqPoseidonHashGadget>(generate_inputs(ins));
            algebraic_sponge_gadget_native_test::<_, _, TweedleFqPoseidonSpongeGadget>(generate_inputs(ins));
        }
    }
}