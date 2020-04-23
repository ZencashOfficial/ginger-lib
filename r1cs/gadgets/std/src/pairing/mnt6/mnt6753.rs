use crate::pairing::mnt6::{
    MNT6PairingGadget, MNT6ConstantPairingGadget,
};
use algebra::curves::mnt6753::MNT6_753Parameters;

pub type MNT6753PairingGadget = MNT6PairingGadget<MNT6_753Parameters>;
pub type MNT6753ConstantPairingGadget = MNT6ConstantPairingGadget<MNT6_753Parameters>;
