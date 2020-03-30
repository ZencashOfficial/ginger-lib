use algebra::{Field, ToBytes, FromBytes};
use rand::Rng;
use std::{hash::Hash, fmt::Debug};
use crate::Error;

pub mod ecvrf;

pub trait FieldBasedVrf {
    type Data: Field;
    type PublicKey: ToBytes + Hash + Eq + Clone + Default + Send + Sync;
    type SecretKey: ToBytes + Clone + Default;
    type Proof: Copy + Clone + Default + Send + Sync + Debug + Eq + PartialEq + ToBytes + FromBytes;
    type GHParams: Clone + Default;

    fn keygen<R: Rng>(
        rng: &mut R,
    ) -> (Self::PublicKey, Self::SecretKey);

    fn prove<R: Rng>
    (
        rng:     &mut R,
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        sk:      &Self::SecretKey,
        message: &[Self::Data],
    ) -> Result<Self::Proof, Error>;

    fn verify
    (
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        message: &[Self::Data],
        proof:   &Self::Proof,
    ) -> Result<Self::Data, Error>;

    fn keyverify(pk: &Self::PublicKey) -> bool;
}