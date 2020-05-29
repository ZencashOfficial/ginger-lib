//! Base field Fq, exponent field Fr and embedding field Fq4 for the MNT4-298.
//!
//! The construction of the degree 4 extension of Fq (with q=1 mod 4) is based on a
//! non-square 13 by successively adding roots of X^4-13.

pub mod fr;
pub use self::fr::*;

pub mod fq;
pub use self::fq::*;

pub mod fq2;
pub use self::fq2::*;

pub mod fq4;
pub use self::fq4::*;

#[cfg(test)]
mod tests;
