//! Montgomery arithmetics for large integers and models of extension fields.
//!
//! - Fp256, Fp320, Fp384 capturing moduli of length 4,5 and 6 words of 64 bits,
//! - Fp768 and Fp832 capturing moduli of length 12 and 13 words of 64 bits,
//! - Quadratic and cubic extensions of prime fields, degree 6 and 12 extension by towering.

pub mod fp_256;
pub use self::fp_256::*;

pub mod fp_320;
pub use self::fp_320::*;

pub mod fp_384;
pub use self::fp_384::*;

pub mod fp_768;
pub use self::fp_768::*;

pub mod fp_832;
pub use self::fp_832::*;

pub mod fp2;
pub use self::fp2::*;

pub mod fp3;
pub use self::fp3::*;

pub mod fp4;
pub use self::fp4::*;

pub mod fp6_2over3;
pub use self::fp6_2over3::*;

pub mod fp6_3over2;

pub mod fp12_2over3over2;
