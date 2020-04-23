use algebra::{AffineCurve, ProjectiveCurve, Field, PairingEngine, ToConstraintField};
use proof_systems::groth16::PreparedVerifyingKey;
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;
use r1cs_std::pairing::ConstantPairingGadget;
use crate::{
    NIZKVerifierConstantGadget,
    groth16::{Groth16, Groth16VerifierGadget, PreparedVerifyingKeyGadget, ProofGadget},
};

use std::marker::PhantomData;

#[derive(Derivative)]
#[derivative(Clone(
bound = "CP::G1ConstantGadget: Clone, CP::GTConstantGadget: Clone, CP::G2PreparedConstantGadget: Clone, "
), Eq, PartialEq, Debug)]
pub struct ConstantPreparedVerifyingKeyGadget<
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
    CP: ConstantPairingGadget<PairingE, ConstraintF, P>,
> {
    pub alpha_g1_beta_g2: CP::GTConstantGadget,
    pub gamma_g2_neg_pc:  CP::G2PreparedConstantGadget,
    pub delta_g2_neg_pc:  CP::G2PreparedConstantGadget,
    pub gamma_abc_g1:     Vec<CP::G1ConstantGadget>,
}


pub struct Groth16ConstantVerifierGadget<PairingE, ConstraintF, P, CP>
    where
        PairingE: PairingEngine,
        ConstraintF: Field,
        P: PairingGadget<PairingE, ConstraintF>,
        CP: ConstantPairingGadget<PairingE, ConstraintF, P>,
{
    _pairing_engine:          PhantomData<PairingE>,
    _engine:                  PhantomData<ConstraintF>,
    _pairing_gadget:          PhantomData<P>,
    _constant_pairing_gadget: PhantomData<CP>,
}

impl<PairingE, ConstraintF, P, CP, C, V> NIZKVerifierConstantGadget<
    Groth16<PairingE, C, V>,
    ConstraintF,
    Groth16VerifierGadget<PairingE, ConstraintF, P>,
> for Groth16ConstantVerifierGadget<PairingE, ConstraintF, P, CP>
    where
        PairingE: PairingEngine,
        ConstraintF: Field,
        C: ConstraintSynthesizer<PairingE::Fr>,
        V: ToConstraintField<PairingE::Fr>,
        P: PairingGadget<PairingE, ConstraintF>,
        CP: ConstantPairingGadget<PairingE, ConstraintF, P>,
{
    type PreparedVerificationKeyGadget = PreparedVerifyingKeyGadget<PairingE, ConstraintF, P>;
    type PreparedVerificationKeyConstantGadget = ConstantPreparedVerifyingKeyGadget<PairingE, ConstraintF, P, CP>;

    fn check_verify_with_constant_pvk<'a, CS, I, T>(
        mut cs: CS,
        pvk: &'a PreparedVerifyingKey<PairingE>,
        mut public_inputs: I,
        proof: &ProofGadget<PairingE, ConstraintF, P>,
    ) -> Result<(), SynthesisError>
        where
            CS: ConstraintSystem<ConstraintF>,
            I: Iterator<Item=&'a T>,
            T: 'a + ToBitsGadget<ConstraintF> + ?Sized,
    {
        let pvk =  Self::PreparedVerificationKeyConstantGadget::from_value(cs.ns(|| "hardcode pvk"), pvk);
        let g_ic = {
            let mut cs = cs.ns(|| "Process input");
            let mut g_ic = pvk.gamma_abc_g1[0].clone();
            let mut input_len = 1;
            for (i, (input, b)) in public_inputs
                .by_ref()
                .zip(pvk.gamma_abc_g1.iter().skip(1))
                .enumerate()
                {
                    let input_bits = input.to_bits(cs.ns(|| format!("Input {}", i)))?;
                    g_ic = CP::G1ConstantGadget::mul_bits_fixed_base(
                        &b.get_constant(),
                        cs.ns(|| format!("Mul {}", i)),
                        &g_ic,
                        input_bits.as_slice()
                    )?;
                    input_len += 1;
                }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == pvk.gamma_abc_g1.len() && public_inputs.next().is_none());
            g_ic
        };

        let test_exp = {
            let proof_a_prep = P::prepare_g1(cs.ns(|| "Prepare proof a"), &proof.a)?;
            let proof_b_prep = P::prepare_g2(cs.ns(|| "Prepare proof b"), &proof.b)?;
            let proof_c_prep = P::prepare_g1(cs.ns(|| "Prepare proof c"), &proof.c)?;

            let g_ic_prep = CP::prepare_g1(cs.ns(|| "Prepare g_ic"), &g_ic)?;
            let ml_1 = P::miller_loop(
                cs.ns(|| "Miller loop 1"),
                &[proof_a_prep],
                &[proof_b_prep],
            )?;
            let ml_2 = CP::miller_loop_with_constant_q(
                cs.ns(|| "Miller loop pc"),
                &[g_ic_prep, proof_c_prep],
                &[pvk.gamma_g2_neg_pc, pvk.delta_g2_neg_pc],
            )?;
            ml_1.mul(cs.ns(|| "ML 1 * ML 2"), &ml_2)?
        };

        let test = CP::final_exponentiation(cs.ns(|| "Final Exp"), &test_exp).unwrap();

        test.enforce_equal(cs.ns(|| "Test 1"), &pvk.alpha_g1_beta_g2)?;
        Ok(())
    }
}

impl <PairingE, ConstraintF, P, CP> Into<PreparedVerifyingKeyGadget<PairingE, ConstraintF, P>>
for ConstantPreparedVerifyingKeyGadget<PairingE, ConstraintF, P, CP>
    where
        PairingE: PairingEngine,
        ConstraintF: Field,
        P: PairingGadget<PairingE, ConstraintF>,
        CP: ConstantPairingGadget<PairingE, ConstraintF, P>
{
    fn into(self) -> PreparedVerifyingKeyGadget<PairingE, ConstraintF, P> {
        PreparedVerifyingKeyGadget::<PairingE, ConstraintF, P>{
            alpha_g1_beta_g2: CP::cast_to_gt_gadget(&self.alpha_g1_beta_g2),
            gamma_g2_neg_pc: CP::cast_to_g2_prepared_gadget(&self.gamma_g2_neg_pc),
            delta_g2_neg_pc: CP::cast_to_g2_prepared_gadget(&self.delta_g2_neg_pc),
            gamma_abc_g1: self.gamma_abc_g1.iter().map(
                |c_g1| CP::cast_to_g1_gadget(c_g1)
            ).collect::<Vec<_>>(),
        }
    }
}

impl<PairingE, ConstraintF, P, CP> ConstantGadget<PreparedVerifyingKey<PairingE>, ConstraintF>
for ConstantPreparedVerifyingKeyGadget<PairingE, ConstraintF, P, CP>
    where
        PairingE: PairingEngine,
        ConstraintF: Field,
        P: PairingGadget<PairingE, ConstraintF>,
        CP: ConstantPairingGadget<PairingE, ConstraintF, P>,
{
    fn from_value<CS: ConstraintSystem<ConstraintF>>(mut cs: CS, value: &PreparedVerifyingKey<PairingE>) -> Self
    {

        let alpha_g1_beta_g2 = CP::GTConstantGadget::from_value(cs.ns(|| "hardcode alpha_g1_beta_g2"), &value.alpha_g1_beta_g2);
        let gamma_g2_neg_pc = CP::G2PreparedConstantGadget::from_value(cs.ns(|| "hardcode gamma_g2_neg_pc"), &value.gamma_g2_neg_pc);
        let delta_g2_neg_pc = CP::G2PreparedConstantGadget::from_value(cs.ns(|| "hardcode delta_g2_neg_pc"), &value.delta_g2_neg_pc);

        let ic = value.gamma_abc_g1
            .iter()
            .enumerate()
            .map(|(i, query_i)| {
                CP::G1ConstantGadget::from_value(cs.ns(|| format!("hardcode query_{}", i)),
                                                 &query_i.into_projective()
                )
            })
            .collect::<Vec<_>>();

        Self {
            alpha_g1_beta_g2,
            gamma_g2_neg_pc,
            delta_g2_neg_pc,
            gamma_abc_g1: ic,
        }
    }

    fn get_constant(&self) -> PreparedVerifyingKey<PairingE> {
        PreparedVerifyingKey::<PairingE>{
            alpha_g1_beta_g2: self.alpha_g1_beta_g2.get_constant(),
            gamma_g2_neg_pc: self.gamma_g2_neg_pc.get_constant(),
            delta_g2_neg_pc: self.delta_g2_neg_pc.get_constant(),
            gamma_abc_g1: self.gamma_abc_g1.iter().map(
                |g1| g1.get_constant().into_affine()
            ).collect::<Vec<_>>(),
        }
    }
}

impl<PairingE, ConstraintF, P, CP> CondSelectGadget<ConstraintF>
for ConstantPreparedVerifyingKeyGadget<PairingE, ConstraintF, P, CP>
    where
        PairingE: PairingEngine,
        ConstraintF: Field,
        P: PairingGadget<PairingE, ConstraintF>,
        CP: ConstantPairingGadget<PairingE, ConstraintF, P>,
{
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self
    ) -> Result<Self, SynthesisError> {
        let alpha_g1_beta_g2 = CP::GTConstantGadget::conditionally_select(
            cs.ns(|| "alpha_g1_beta_g2"),
            cond,
            &first.alpha_g1_beta_g2,
            &second.alpha_g1_beta_g2,
        )?;

        let gamma_g2_neg_pc = CP::G2PreparedConstantGadget::conditionally_select(
            cs.ns(|| "gamma_g2_neg_pc"),
            cond,
            &first.gamma_g2_neg_pc,
            &second.gamma_g2_neg_pc,
        )?;

        let delta_g2_neg_pc = CP::G2PreparedConstantGadget::conditionally_select(
            cs.ns(|| "delta_g2_neg_pc"),
            cond,
            &first.delta_g2_neg_pc,
            &second.delta_g2_neg_pc,
        )?;

        let mut gamma_abc_g1 = Vec::new();

        //TODO: Is this needed ? If yes, should we enforce it in the circuit too?
        assert_eq!(first.gamma_abc_g1.len(), second.gamma_abc_g1.len());

        for (i, (first, second)) in
            first.gamma_abc_g1.iter().zip(second.gamma_abc_g1.iter()).enumerate() {
            let val = CP::G1ConstantGadget::conditionally_select(
                cs.ns(|| format!("gamma_abc_g1_{}", i)),
                cond,
                &first,
                &second,
            )?;
            gamma_abc_g1.push(val);
        }

        Ok(Self {
            alpha_g1_beta_g2,
            gamma_g2_neg_pc,
            delta_g2_neg_pc,
            gamma_abc_g1,
        })
    }

    fn cost() -> usize {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use algebra::{
        curves::mnt4753::MNT4 as BaseCurve,
        fields::{
            mnt4753::Fr as BaseField,
            mnt6753::Fr as VerificationField
        },
        BitIterator, PrimeField,
    };

    use proof_systems::groth16::*;

    use r1cs_std::{
        pairing::mnt4753::{
            MNT4753PairingGadget as PairingGadget,
            MNT4753ConstantPairingGadget as ConstantPairingGadget,
        },
        boolean::Boolean, test_constraint_system::TestConstraintSystem
    };

    use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

    use crate::nizk::groth16::{
        hardcoded::*, *
    };

    use rand::{thread_rng, Rng};

    type TestProofSystem = Groth16<BaseCurve, Bench<BaseField>, BaseField>;
    type TestVerifierGadget = Groth16VerifierGadget<BaseCurve, VerificationField, PairingGadget>;
    type TestConstantVerifierGadget = Groth16ConstantVerifierGadget<
        BaseCurve,
        VerificationField,
        PairingGadget,
        ConstantPairingGadget
    >;
    type TestProofGadget = ProofGadget<BaseCurve, VerificationField, PairingGadget>;

    struct Bench<F: Field> {
        inputs: Vec<Option<F>>,
        num_constraints: usize,
    }

    impl<F: Field> ConstraintSynthesizer<F> for Bench<F> {
        fn generate_constraints<CS: ConstraintSystem<F>>(
            self,
            cs: &mut CS,
        ) -> Result<(), SynthesisError> {
            assert!(self.inputs.len() >= 2);
            assert!(self.num_constraints >= self.inputs.len());

            let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
            for (i, input) in self.inputs.into_iter().enumerate() {
                let input_var = cs.alloc_input(
                    || format!("Input {}", i),
                    || input.ok_or(SynthesisError::AssignmentMissing),
                )?;
                variables.push((input, input_var));
            }

            for i in 0..self.num_constraints {
                let new_entry = {
                    let (input_1_val, input_1_var) = variables[i];
                    let (input_2_val, input_2_var) = variables[i + 1];
                    let result_val = input_1_val
                        .and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
                    let result_var = cs.alloc(
                        || format!("Result {}", i),
                        || result_val.ok_or(SynthesisError::AssignmentMissing),
                    )?;
                    cs.enforce(
                        || format!("Enforce constraint {}", i),
                        |lc| lc + input_1_var,
                        |lc| lc + input_2_var,
                        |lc| lc + result_var,
                    );
                    (result_val, result_var)
                };
                variables.push(new_entry);
            }
            Ok(())
        }
    }

    #[test]
    fn groth16_verifier_hardcoded_vk_test() {
        let num_inputs = 5;
        let num_constraints = num_inputs;
        let rng = &mut thread_rng();
        let mut inputs: Vec<Option<BaseField>> = Vec::with_capacity(num_inputs);
        for _ in 0..num_inputs {
            inputs.push(Some(rng.gen()));
        }
        let params = {
            let c = Bench::<BaseField> {
                inputs: vec![None; num_inputs],
                num_constraints,
            };

            generate_random_parameters(c, rng).unwrap()
        };

        let pvk = prepare_verifying_key(&params.vk);

        {
            let proof = {
                // Create an instance of our circuit (with the
                // witness)
                let c = Bench {
                    inputs: inputs.clone(),
                    num_constraints,
                };
                // Create a groth16 proof with our parameters.
                create_random_proof(c, &params, rng).unwrap()
            };

            // assert!(!verify_proof(&pvk, &proof, &[a]).unwrap());
            let mut cs = TestConstraintSystem::<VerificationField>::new();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let mut cs = cs.ns(|| "Allocate Input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let mut input_bits = BitIterator::new(input.into_repr()).collect::<Vec<_>>();
                    // Input must be in little-endian, but BitIterator outputs in big-endian.
                    input_bits.reverse();

                    let input_bits =
                        Vec::<Boolean>::alloc_input(cs.ns(|| format!("Input {}", i)), || {
                            Ok(input_bits)
                        })
                            .unwrap();
                    input_gadgets.push(input_bits);
                }
            }

            let proof_gadget =
                TestProofGadget::alloc(cs.ns(|| "Proof"), || Ok(proof.clone())).unwrap();
            println!("Time to verify!\n\n\n\n");
            <TestConstantVerifierGadget as
            NIZKVerifierConstantGadget<TestProofSystem, VerificationField, TestVerifierGadget>
            >::check_verify_with_constant_pvk(
                cs.ns(|| "Verify"),
                &pvk,
                input_gadgets.iter(),
                &proof_gadget,
            ).unwrap();
            if !cs.is_satisfied() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            // cs.print_named_objects();
            assert!(cs.is_satisfied());

            //Negative test. Change inputs:

            let mut inputs: Vec<Option<BaseField>> = Vec::with_capacity(num_inputs);
            for _ in 0..num_inputs {
                inputs.push(Some(rng.gen()));
            }

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let mut cs = cs.ns(|| "Allocate wrong input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let mut input_bits = BitIterator::new(input.into_repr()).collect::<Vec<_>>();
                    // Input must be in little-endian, but BitIterator outputs in big-endian.
                    input_bits.reverse();

                    let input_bits =
                        Vec::<Boolean>::alloc_input(cs.ns(|| format!("Input {}", i)), || {
                            Ok(input_bits)
                        })
                            .unwrap();
                    input_gadgets.push(input_bits);
                }
            }

            <TestConstantVerifierGadget as
            NIZKVerifierConstantGadget<TestProofSystem, VerificationField, TestVerifierGadget>
            >::check_verify_with_constant_pvk(
                cs.ns(|| "Negative Verify"),
                &pvk,
                input_gadgets.iter(),
                &proof_gadget,
            ).unwrap();

            assert!(!cs.is_satisfied());
        }
    }

    #[test]
    fn groth16_verifier_hardcode_and_conditionally_select_test() {

        let num_inputs = 5;
        let num_constraints = 100;
        let rng = &mut thread_rng();
        let mut inputs: Vec<Option<BaseField>> = Vec::with_capacity(num_inputs);
        for _ in 0..num_inputs {
            inputs.push(Some(rng.gen()));
        }

        let mut create_test_values = || -> (Proof<BaseCurve>, PreparedVerifyingKey<BaseCurve>) {
            let params = {
                let c = Bench::<BaseField> {
                    inputs: vec![None; num_inputs],
                    num_constraints,
                };
                generate_random_parameters(c, rng).unwrap()
            };
            let pvk = prepare_verifying_key(&params.vk);

            let proof = {
                let c = Bench::<BaseField> {
                    inputs: inputs.clone(),
                    num_constraints,
                };
                create_random_proof(c, &params, rng).unwrap()
            };

            (proof, pvk)
        };

        let (proof1, pvk_1) = create_test_values();

        let (proof2, pvk_2) = create_test_values();

        //Test for proof1 and pvk1
        {
            let mut cs = TestConstraintSystem::<VerificationField>::new();

            let inputs: Vec<_> = inputs.iter().map(|&input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let mut cs = cs.ns(|| "Allocate Input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let mut input_bits = BitIterator::new(input.into_repr()).collect::<Vec<_>>();
                    // Input must be in little-endian, but BitIterator outputs in big-endian.
                    input_bits.reverse();

                    let input_bits =
                        Vec::<Boolean>::alloc_input(cs.ns(|| format!("Input {}", i)), || {
                            Ok(input_bits)
                        })
                            .unwrap();
                    input_gadgets.push(input_bits);
                }
            }

            let proof_gadget =
                TestProofGadget::alloc(cs.ns(|| "Proof"), || Ok(proof1.clone())).unwrap();

            <TestConstantVerifierGadget as
            NIZKVerifierConstantGadget<TestProofSystem, VerificationField, TestVerifierGadget>
            >::check_verify_with_constant_pvk_selection(
                cs.ns(|| "select pvk and verify proof 1"),
                &pvk_1,
                &pvk_2,
                &Boolean::constant(true),
                input_gadgets.iter(),
                &proof_gadget
            ).unwrap();

            if !cs.is_satisfied() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            assert!(cs.is_satisfied());

            //Negative test
            <TestConstantVerifierGadget as
            NIZKVerifierConstantGadget<TestProofSystem, VerificationField, TestVerifierGadget>
            >::check_verify_with_constant_pvk_selection(
                cs.ns(|| "negative select pvk and verify proof 1"),
                &pvk_1,
                &pvk_2,
                &Boolean::constant(false),
                input_gadgets.iter(),
                &proof_gadget
            ).unwrap();

            assert!(!cs.is_satisfied());
        }

        //Test for proof2 and pvk2
        {
            let mut cs = TestConstraintSystem::<VerificationField>::new();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let mut cs = cs.ns(|| "Allocate Input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let mut input_bits = BitIterator::new(input.into_repr()).collect::<Vec<_>>();
                    // Input must be in little-endian, but BitIterator outputs in big-endian.
                    input_bits.reverse();

                    let input_bits =
                        Vec::<Boolean>::alloc_input(cs.ns(|| format!("Input {}", i)), || {
                            Ok(input_bits)
                        })
                            .unwrap();
                    input_gadgets.push(input_bits);
                }
            }

            let proof_gadget =
                TestProofGadget::alloc(cs.ns(|| "Proof"), || Ok(proof2.clone())).unwrap();

            <TestConstantVerifierGadget as
            NIZKVerifierConstantGadget<TestProofSystem, VerificationField, TestVerifierGadget>
            >::check_verify_with_constant_pvk_selection(
                cs.ns(|| "select pvk and verify proof 2"),
                &pvk_1,
                &pvk_2,
                &Boolean::constant(false),
                input_gadgets.iter(),
                &proof_gadget
            ).unwrap();

            if !cs.is_satisfied() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            assert!(cs.is_satisfied());

            //Negative test
            <TestConstantVerifierGadget as
            NIZKVerifierConstantGadget<TestProofSystem, VerificationField, TestVerifierGadget>
            >::check_verify_with_constant_pvk_selection(
                cs.ns(|| "negative select pvk and verify proof 2"),
                &pvk_1,
                &pvk_2,
                &Boolean::constant(true),
                input_gadgets.iter(),
                &proof_gadget
            ).unwrap();

            assert!(!cs.is_satisfied());
        }
    }
}