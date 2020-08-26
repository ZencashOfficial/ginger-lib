//! Base field Fq, exponent field Fr and embedding field Fq6 for the MNT6-753.
//!
//! The construction of the degree 6 extension of Fq (with q=1 mod 12) is based on a
//! non-square and non-cube 11 by successively adding roots of X^6 - 11.

pub mod fr;
pub use self::fr::*;

pub mod fq;
pub use self::fq::*;

pub mod fq3;
pub use self::fq3::*;

pub mod fq6;
pub use self::fq6::*;

#[cfg(test)]
mod tests;
