use algebra::PrimeField;
use crate::FieldBasedHashParameters;

//TODO: scalar_mul and matrix_mix are fine in this trait or should be private impl functions ?
pub trait SBox {
    type Field: PrimeField;
    type Parameters: FieldBasedHashParameters<Fr = Self::Field>;

    // Function that does the scalar multiplication
    fn scalar_mul(res: &mut Self::Field, state: &mut [Self::Field], start_idx_cst: usize);

    // Function that does the mix matrix
    fn matrix_mix(state: &mut Vec<Self::Field>);

    // Apply this SBox to the state, if performing a full round
    fn apply_full(state: &mut Vec<Self::Field>, last: bool);

    // Apply this SBox to the state, if performing a partial round
    fn apply_partial(state: &mut Vec<Self::Field>);
}

pub trait BatchSBox: SBox {

    fn apply_full_batch(vec_state: &mut [Vec<Self::Field>], last: bool) {
        vec_state.iter_mut().for_each(|s| Self::apply_full(s, last));
    }

    fn apply_partial_batch(vec_state: &mut [Vec<Self::Field>]) {
        vec_state.iter_mut().for_each(|s| Self::apply_partial(s));
    }
}