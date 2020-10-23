use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use r1cs_std::Assignment;

/*use r1cs_std::{
    fields::fp::FpGadget,
    alloc::AllocGadget,
};
use r1cs_std::fields::FieldGadget;

type MNT4753FrGadget = FpGadget<MNT4753Fr>;*/


// Let's create a circuit to prove that I know a x s.t. x^2= y
// Witness: x
// Public Input: y
// Statement: x * x == y

// Let's create a struct that holds our witnesses
struct SimpleSquareCircuit {
    x: Option<MNT4753Fr>,
}

// The circuit, made up of witness, public inputs and constraints.
// The abstraction used by Zexe is to provide a "ConstraintSynthesizer" interface with a
// generate_constraints() function that holds the internal logic of the circuit.
// The Generator and the Prover will run the generate_constraints() function, "save" all the inputs,
// witnesses and constraints that they will "meet", and produce proving key, verifying key,
// and the proof respectively (i.e. the generate_constraint() functions takes as input a ConstraintSystem
// Prover and Generator will implement ConstraintSystem trait and provide their own implementation
// for alloc(), alloc_input() and enforce() functions).

// We implement ConstraintSynthesizer for our SimpleSquareCircuit. The struct for which
// ConstraintSynthesizer will be implemented is a sort of description of our circuit in terms
// of the witnesses, public inputs it is made of, plus maybe constants or other variables useful
// to the generation of the constraints
// Note that we need to be sure that the value of the variables and the arithmetics is performed
// on the scalar field of the curve we are using to compute the proof (e.g. the pairing-friendly
// curve).
impl ConstraintSynthesizer<MNT4753Fr> for SimpleSquareCircuit {
    fn generate_constraints<CS: ConstraintSystem<MNT4753Fr>>(self, cs: &mut CS) -> Result<(), SynthesisError> {

        // Create a witness variable, and assign to it the value of x
        let x_val = self.x;

        // "alloc" is used to define a new witness
        let x_var = cs.alloc(
            // a function whose output must be a (unique) name to assign to the variable
            // this is just for debugging purposes: internally each variable has assigned
            // a unique index, used to identify it.
            || "x",
            // a function whose output is the value to assign to the variable. Note that
            // at (pk, vk) generation time this function is totally ignored: the Generator
            // algorithm doesn't need to know the value of the variables, but just the
            // architecture of the circuit in terms of what are the inputs, witnesses and the
            // way these are connected together through the constraints. The prover MUST
            // supply a value for each (witness and input) variable instead, because it has
            // to prove that he knows a variables assignment that satisfies the circuit.
            || x_val.ok_or(SynthesisError::AssignmentMissing)
        )?;

        // Create a public input variable, and assign to it the value of x^2. The prover needs to
        // compute the value of the input variables anyway (in our case x_val * x_val).
        // At verifier time, of course, y_var will be assigned with the value supplied by
        // the verifier, and not with the value supplied by the function below like in the previous
        // case.
        let y_var = cs.alloc_input(
            || "y",
            || Ok(x_val.get()? * &x_val.get()?)
        )?;

        // Finally, let's enforce that x * x == y, through a R1CS constraint.
        // Notation here is a little bit confusing, just ignore for the moment
        // the "lc" that the functions take as input: that is related to Generator
        // and Prover internal implementation details.
        cs.enforce(
            || "x * x == y",
            |lc| lc + x_var, // A
            |lc| lc + x_var, // B
            |lc| lc + y_var // C
        );

        ////////// WITH GADGETS
        /*let x_g = MNT4753FrGadget::alloc(
            cs.ns(|| "alloc x"),
            || x_val.ok_or(SynthesisError::AssignmentMissing)
        )?;


        let y_g = MNT4753FrGadget::alloc_input(
            cs.ns(|| "alloc input y"),
            || Ok(x_val.get()? * &x_val.get()?)
        )?;

        x_g.mul_equals(cs.ns(|| "x * x = y"), &x_g, &y_g)?;*/

        Ok(())
    }
}

use algebra::UniformRand;
use proof_systems::groth16::{generate_random_parameters, create_random_proof, prepare_verifying_key, verify_proof};
use algebra::curves::mnt4753::MNT4;
use rand::rngs::OsRng;

fn main() {
    let mut rng = OsRng::default();

    // Let's create proving key and verifying key for our circuit. "params" is a struct
    // holding them.
    let params = {

        // We don't need to supply a value for x: the Generator doesn't need it and
        // doesn't use it anyway.
        let circuit = SimpleSquareCircuit{ x: None };

        // Generate proving key and verifying key using a OsRng to generate secret values.
        // Normally, the (pk, vk) generation requires a MPC protocol to generate secret values
        // being sure no-one cheated. Secret values must be tossed away (toxic waste) immediately
        // after, and this seems to be a moment they particularly enjoy in ceremonies.
        generate_random_parameters::<MNT4, _, _>(circuit, &mut rng).unwrap()
    };

    // Let's create our proof
    let x = MNT4753Fr::rand(&mut rng); // Generate a random x

    let proof = {

        // This time we need to supply a value for x
        let circuit = SimpleSquareCircuit{ x: Some(x) };

        // Let's create the proof and generate the prover secrets values (if ZK is required)
        // through a OsRng
        create_random_proof(circuit, &params, &mut rng).unwrap()
    };

    // Let's verify our proof

    // We can do some precomputations on the vk points, that will be useful to make the
    // pairing operations more efficient.
    let pvk = prepare_verifying_key(&params.vk);

    let y = x * &x; // Compute our public input

    // Verify our proof using pvk and y
    assert!(verify_proof(&pvk, &proof, &[y]).unwrap());

    // Let's be sure that we put the constraints correctly: compute y = x^3
    // and check that the proof verification fails
    assert!(!verify_proof(&pvk, &proof, &[y * &x]).unwrap());
}
