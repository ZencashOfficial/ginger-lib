use crate::pairing::mnt4::{
    MNT4PairingGadget, MNT4ConstantPairingGadget
};
use algebra::curves::mnt4753::MNT4_753Parameters;

pub type MNT4753PairingGadget = MNT4PairingGadget<MNT4_753Parameters>;
pub type MNT4753ConstantPairingGadget = MNT4ConstantPairingGadget<MNT4_753Parameters>;