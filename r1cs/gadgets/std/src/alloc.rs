use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::borrow::Borrow;

pub trait AllocGadget<V, ConstraintF: Field>
where
    Self: Sized,
    V: ?Sized,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_without_check<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
        where
            F: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<V>, { Self::alloc(cs, f) }

    fn alloc_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>,
    {
        Self::alloc(cs, f)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_input_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>,
    {
        Self::alloc_input(cs, f)
    }
}

impl<I, ConstraintF: Field, A: AllocGadget<I, ConstraintF>> AllocGadget<[I], ConstraintF>
    for Vec<A>
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc(&mut cs.ns(|| format!("value_{}", i)), || {
                Ok(value)
            })?);
        }
        Ok(vec)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc_input(
                &mut cs.ns(|| format!("value_{}", i)),
                || Ok(value),
            )?);
        }
        Ok(vec)
    }

    fn alloc_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc_checked(
                &mut cs.ns(|| format!("value_{}", i)),
                || Ok(value),
            )?);
        }
        Ok(vec)
    }

    fn alloc_input_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc_input_checked(
                &mut cs.ns(|| format!("value_{}", i)),
                || Ok(value),
            )?);
        }
        Ok(vec)
    }
}

/// Get a Gadget from the corresponding constant. At low level, the constant
/// will be the coefficient of the CS::one() variable.
pub trait ConstantGadget<V, ConstraintF: Field>
    where
        Self: Sized,
        V: Sized ,
{
    fn from_value<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        value: &V
    ) -> Self;

    fn get_constant(&self) -> V;
}

pub trait FromGadget<T, ConstraintF: Field>
    where
        Self: Sized,
        T: Sized,
{
    fn from<CS: ConstraintSystem<ConstraintF>>(_: T, cs: CS) -> Result<Self, SynthesisError>;
}