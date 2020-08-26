use crate::bits::{boolean::Boolean, uint8::UInt8};
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};

pub mod boolean;
pub mod uint32;
pub mod uint8;

/* Provides the interfaces for the conversion circuits ("unpacking") for
field element gadgets into a vector of Booleans.
*/
pub trait ToBitsGadget<ConstraintF: Field> {
    /* Interface for insecure but inexpensive unpacking, does not enforce the resulting Boolean
    vector to be the integer representation of the field element, unless some extra conditions are
    met.
    */
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError>;

    /* Interface for the secure unpacking of field elements.
    */
    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError>;
}

/* Provides the interfaces for the conversion circuits of field elements into a vectors of UInt8.
*/
pub trait ToBytesGadget<ConstraintF: Field> {
    /* Interface for insecure but inexpensive unpacking, does not enforce the resulting UInt8
    vector to be the integer representation of the field element, unless some extra conditions are
    met.
    */
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError>;

    /* Interface for the secure unpacking of field elements.
    */
    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError>;
}

/* Provides the interface for the secure conversion circuit ("packing") of a slice of Boolean
into a field element gadget.
*/
pub trait FromBitsGadget<ConstraintF: Field>
    where
        Self: Sized
{
    /// Given a bit representation `bits` of bit len not bigger than CAPACITY
    /// (i.e. MODULUS - 1) of `Self` in *big endian* form, reconstructs a `Self`.
    fn from_bits<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        bits: &[Boolean],
    ) -> Result<Self, SynthesisError>;
}

/* Provides the interface used for point compression
*/
pub trait ToCompressedBitsGadget<ConstraintF: Field> {

    /// Enforce compression of an element through serialization of the x coordinate and storing
    /// a sign bit for the y coordinate. For GT elements we assume x <-> c1 and y <-> c0 to avoid
    /// confusion. When enforcing byte serialization of a field element, "x_in_field" and "y_in_field"
    /// flags could be set in order to enforce too that their bit representation is under the
    /// field modulus (default behaviour is both set to false).
    fn to_compressed<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError>;
}


/* Implementations of the above traits (except the ToCompressedBitsGadget) for Booleans,
slices and vectors of Booleans.
*/

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for Boolean {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for [Boolean] {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }
}
impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for Vec<Boolean> {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }
}

/* Implementations of the above traits (except the ToCompressedBitsGadget) for UInt8,
slices and vectors of UInt8 gadgets.
*/
impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for [UInt8] {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut result = Vec::with_capacity(&self.len() * 8);
        for byte in self {
            result.extend_from_slice(&byte.into_bits_le());
        }
        Ok(result)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        self.to_bits(cs)
    }
}

impl<ConstraintF: Field> ToBytesGadget<ConstraintF> for [UInt8] {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

impl<'a, ConstraintF: Field> ToBytesGadget<ConstraintF> for &'a [UInt8] {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

impl<'a, ConstraintF: Field, T: 'a + ToBytesGadget<ConstraintF>> ToBytesGadget<ConstraintF>
for &'a T
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        (*self).to_bytes(cs)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}
