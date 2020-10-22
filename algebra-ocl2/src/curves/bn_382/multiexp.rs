use bit_vec::{self, BitVec};
use futures::Future;
use log::{info, warn};
use std::io;
use std::iter;
use std::sync::Arc;

use algebra::{BigInteger, BigInteger384, Field, PrimeField, AffineCurve, ProjectiveCurve};
use algebra::curves::bn_382::{G1Affine, G1Projective};
use algebra::fields::bn_382::Fr;

use super::multicore::Worker;
use super::gpu::{LockedMultiexpKernel, MultiexpKernel};
use crate::SynthesisError;

/// An object that builds a source of bases.
pub trait SourceBuilder: Send + Sync + 'static + Clone {
    type Source: Source;

    fn new(self) -> Self::Source;
    fn get(self) -> (Arc<Vec<G1Affine>>, usize);
}

/// A source of bases, like an iterator.
pub trait Source {
    /// Parses the element from the source. Fails if the point is at infinity.
    fn add_assign_mixed(
        &mut self,
        to: &mut G1Projective,
    ) -> Result<(), SynthesisError>;

    /// Skips `amt` elements from the source, avoiding deserialization.
    fn skip(&mut self, amt: usize) -> Result<(), SynthesisError>;
}

impl SourceBuilder for (Arc<Vec<G1Affine>>, usize) {
    type Source = (Arc<Vec<G1Affine>>, usize);

    fn new(self) -> (Arc<Vec<G1Affine>>, usize) {
        (self.0.clone(), self.1)
    }

    fn get(self) -> (Arc<Vec<G1Affine>>, usize) {
        (self.0.clone(), self.1)
    }
}

impl Source for (Arc<Vec<G1Affine>>, usize) {
    fn add_assign_mixed(
        &mut self,
        to: &mut G1Projective,
    ) -> Result<(), SynthesisError> {

        if self.0.len() <= self.1 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "expected more bases from source",
            )
            .into());
        }

        if self.0[self.1].is_zero() {
            return Err(SynthesisError::UnexpectedIdentity);
        }

        to.add_assign_mixed(&self.0[self.1]);

        self.1 += 1;

        Ok(())
    }

    fn skip(&mut self, amt: usize) -> Result<(), SynthesisError> {
        if self.0.len() <= self.1 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "expected more bases from source",
            )
            .into());
        }

        self.1 += amt;

        Ok(())
    }
}

pub trait QueryDensity {
    /// Returns whether the base exists.
    type Iter: Iterator<Item = bool>;

    fn iter(self) -> Self::Iter;
    fn get_query_size(self) -> Option<usize>;
}

#[derive(Clone)]
pub struct FullDensity;

impl AsRef<FullDensity> for FullDensity {
    fn as_ref(&self) -> &FullDensity {
        self
    }
}

impl<'a> QueryDensity for &'a FullDensity {
    type Iter = iter::Repeat<bool>;

    fn iter(self) -> Self::Iter {
        iter::repeat(true)
    }

    fn get_query_size(self) -> Option<usize> {
        None
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct DensityTracker {
    pub bv: BitVec,
    pub total_density: usize,
}

impl<'a> QueryDensity for &'a DensityTracker {
    type Iter = bit_vec::Iter<'a>;

    fn iter(self) -> Self::Iter {
        self.bv.iter()
    }

    fn get_query_size(self) -> Option<usize> {
        Some(self.bv.len())
    }
}

impl DensityTracker {
    pub fn new() -> DensityTracker {
        DensityTracker {
            bv: BitVec::new(),
            total_density: 0,
        }
    }

    pub fn add_element(&mut self) {
        self.bv.push(false);
    }

    pub fn inc(&mut self, idx: usize) {
        if !self.bv.get(idx).unwrap() {
            self.bv.set(idx, true);
            self.total_density += 1;
        }
    }

    pub fn get_total_density(&self) -> usize {
        self.total_density
    }

    /// Extend by concatenating `other`. If `is_input_density` is true, then we are tracking an input density,
    /// and other may contain a redundant input for the `One` element. Coalesce those as needed and track the result.
    pub fn extend(&mut self, other: Self, is_input_density: bool) {
        if other.bv.is_empty() {
            // Nothing to do if other is empty.
            return;
        }

        if self.bv.is_empty() {
            // If self is empty, assume other's density.
            self.total_density = other.total_density;
            self.bv = other.bv;
            return;
        }

        if is_input_density {
            // Input densities need special handling to coalesce their first inputs.

            if other.bv[0] {
                // If other's first bit is set,
                if self.bv[0] {
                    // And own first bit is set, then decrement total density so the final sum doesn't overcount.
                    self.total_density -= 1;
                } else {
                    // Otherwise, set own first bit.
                    self.bv.set(0, true);
                }
            }
            // Now discard other's first bit, having accounted for it above, and extend self by remaining bits.
            self.bv.extend(other.bv.iter().skip(1));
        } else {
            // Not an input density, just extend straightforwardly.
            self.bv.extend(other.bv);
        }

        // Since any needed adjustments to total densities have been made, just sum the totals and keep the sum.
        self.total_density += other.total_density;
    }
}

fn multiexp_inner<Q, D, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<BigInteger384>>,
    mut skip: u32,
    c: u32,
    handle_trivial: bool,
) -> Box<dyn Future<Item = G1Projective, Error = SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    S: SourceBuilder,
{
    // Perform this region of the multiexp
    let this = {
        let bases = bases.clone();
        let exponents = exponents.clone();
        let density_map = density_map.clone();

        pool.compute(move || {
            // Accumulate the result
            let mut acc = G1Projective::zero();

            // Build a source for the bases
            let mut bases = bases.new();

            // Create space for the buckets
            let mut buckets = vec![G1Projective::zero(); (1 << c) - 1];

            let zero = Fr::zero().into_repr();
            let one = Fr::one().into_repr();

            // Sort the bases into buckets
            for (&exp, density) in exponents.iter().zip(density_map.as_ref().iter()) {
                if density {
                    if exp == zero {
                        bases.skip(1)?;
                    } else if exp == one {
                        if handle_trivial {
                            bases.add_assign_mixed(&mut acc)?;
                        } else {
                            bases.skip(1)?;
                        }
                    } else {
                        let mut exp = exp;
                        exp.divn(1 << skip);
                        let exp = exp.as_ref()[0] % (1 << c);

                        if exp != 0 {
                            bases.add_assign_mixed(&mut buckets[(exp - 1) as usize])?;
                        } else {
                            bases.skip(1)?;
                        }
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G1Projective::zero();
            for exp in buckets.into_iter().rev() {
                running_sum.add_assign_mixed(&exp.into_affine());
                acc.add_assign_mixed(&running_sum.into_affine());
            }

            Ok(acc)
        })
    };

    skip += c;

    if skip >= 382 {
        // There isn't another region.
        Box::new(this)
    } else {
        // There's another region more significant. Calculate and join it with
        // this region recursively.
        Box::new(
            this.join(multiexp_inner(
                pool,
                bases,
                density_map,
                exponents,
                skip,
                c,
                false,
            ))
            .map(move |(this, mut higher)| {
                for _ in 0..c {
                    higher.double_in_place();
                }

                higher.add_assign_mixed(&this.into_affine());

                higher
            }),
        )
    }
}

/// Perform multi-exponentiation. The caller is responsible for ensuring the
/// query size is the same as the number of exponents.
pub fn multiexp_c<Q, D, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<BigInteger384>>,
    kern: &mut Option<LockedMultiexpKernel>,
    c: usize
) -> Box<dyn Future<Item = G1Projective, Error = SynthesisError>>
    where
            for<'a> &'a Q: QueryDensity,
            D: Send + Sync + 'static + Clone + AsRef<Q>,
            S: SourceBuilder,
{
    if let Some(ref mut kern) = kern {
        if let Ok(p) = kern.with(|k: &mut MultiexpKernel| {
            let mut exps = vec![exponents[0]; exponents.len()];
            let mut n = 0;
            for (&e, d) in exponents.iter().zip(density_map.as_ref().iter()) {
                if d {
                    exps[n] = e;
                    n += 1;
                }
            }

            let (bss, skip) = bases.clone().get();
            k.multiexp_c(pool, bss, Arc::new(exps.clone()), skip, n, c)
        }) {
            return Box::new(pool.compute(move || Ok(p)));
        }
    }

    let c = if exponents.len() < 32 {
        3u32
    } else {
        (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    if let Some(query_size) = density_map.as_ref().get_query_size() {
        // If the density map has a known query size, it should not be
        // inconsistent with the number of exponents.
        assert!(query_size == exponents.len());
    }

    let future = multiexp_inner(pool, bases, density_map, exponents, 0, c, true);
    {
        // Do not give the control back to the caller till the
        // multiexp is done. We may want to reacquire the GPU again
        // between the multiexps.
        let result = future.wait();
        Box::new(pool.compute(move || result))
    }
}

/// Perform multi-exponentiation. The caller is responsible for ensuring the
/// query size is the same as the number of exponents.
pub fn multiexp<Q, D, S>(
    pool: &Worker,
    bases: S,
    density_map: D,
    exponents: Arc<Vec<BigInteger384>>,
    kern: &mut Option<LockedMultiexpKernel>,
) -> Box<dyn Future<Item = G1Projective, Error = SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    S: SourceBuilder,
{
    if let Some(ref mut kern) = kern {
        if let Ok(p) = kern.with(|k: &mut MultiexpKernel| {
            let mut exps = vec![exponents[0]; exponents.len()];
            let mut n = 0;
            for (&e, d) in exponents.iter().zip(density_map.as_ref().iter()) {
                if d {
                    exps[n] = e;
                    n += 1;
                }
            }

            let (bss, skip) = bases.clone().get();
            k.multiexp(pool, bss, Arc::new(exps.clone()), skip, n)
        }) {
            return Box::new(pool.compute(move || Ok(p)));
        }
    }

    let c = if exponents.len() < 32 {
        3u32
    } else {
        (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    if let Some(query_size) = density_map.as_ref().get_query_size() {
        // If the density map has a known query size, it should not be
        // inconsistent with the number of exponents.
        assert!(query_size == exponents.len());
    }

    let future = multiexp_inner(pool, bases, density_map, exponents, 0, c, true);
    {
        // Do not give the control back to the caller till the
        // multiexp is done. We may want to reacquire the GPU again
        // between the multiexps.
        let result = future.wait();
        Box::new(pool.compute(move || result))
    }
}

pub fn create_multiexp_kernel(_log_d: usize, priority: bool) -> Option<MultiexpKernel>
{
    match MultiexpKernel::create(priority) {
        Ok(k) => {
            info!("GPU Multiexp kernel instantiated!");
            Some(k)
        }
        Err(e) => {
            warn!("Cannot instantiate GPU Multiexp kernel! Error: {}", e);
            None
        }
    }
}
