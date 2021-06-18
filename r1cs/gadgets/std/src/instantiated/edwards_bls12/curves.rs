use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::{
    fields::edwards_bls12::fq::Fq,
    curves::edwards_bls12::EdwardsParameters,
};

use crate::edwards_bls12::FqGadget;

pub type EdwardsBlsGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsParameters, EdwardsBlsGadget>();
}
