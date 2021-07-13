use algebra::fields::tweedle::{Fq, Fr};
use primitives::crh::parameters::{
    TweedleFqPoseidonParameters, TweedleFqQuinticSbox,
    TweedleFrPoseidonParameters, TweedleFrQuinticSbox,
};
use crate::crh::{
    sbox::QuinticSBoxGadget,
    poseidon::{PoseidonHashGadget, PoseidonSpongeGadget}
};

type TweedleFqQuinticSboxGadget = QuinticSBoxGadget<Fq, TweedleFqQuinticSbox>;
pub type TweedleFqPoseidonHashGadget = PoseidonHashGadget<
    Fq,
    TweedleFqPoseidonParameters,
    TweedleFqQuinticSbox,
    TweedleFqQuinticSboxGadget
>;
pub type TweedleFqPoseidonSpongeGadget = PoseidonSpongeGadget<
    Fq,
    TweedleFqPoseidonParameters,
    TweedleFqQuinticSbox,
    TweedleFqQuinticSboxGadget
>;

type TweedleFrQuinticSboxGadget = QuinticSBoxGadget<Fr, TweedleFrQuinticSbox>;
pub type TweedleFrPoseidonHashGadget = PoseidonHashGadget<
    Fr,
    TweedleFrPoseidonParameters,
    TweedleFrQuinticSbox,
    TweedleFrQuinticSboxGadget
>;
pub type TweedleFrPoseidonSpongeGadget = PoseidonSpongeGadget<
    Fr,
    TweedleFrPoseidonParameters,
    TweedleFrQuinticSbox,
    TweedleFrQuinticSboxGadget
>;
