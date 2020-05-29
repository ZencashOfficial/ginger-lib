//! Base field Fq, exponent field Fr, and embedding field Fq12 for the BLS12-381.
//!
//! The construction of the degree 12 extension of Fq (with q=1 mod 6) is based on a
//! non-square and non-cube from a quadratic extension.

pub mod fr;
pub use self::fr::*;

pub mod fq;
pub use self::fq::*;

pub mod fq2;
pub use self::fq2::*;

pub mod fq6;
pub use self::fq6::*;

pub mod fq12;
pub use self::fq12::*;

#[cfg(test)]
mod tests;
