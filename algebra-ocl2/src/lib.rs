use std::io;

#[cfg(feature = "gpu")]
pub mod ffgen;

#[cfg(feature = "gpu")]
pub mod gpu;

#[cfg(feature = "gpu")]
pub mod curves;

// pub mod msm;
// pub use self::msm::*;

#[cfg(all(test, feature = "gpu"))]
mod tests;

/// This is an error that could occur during circuit synthesis contexts,
/// such as CRS generation, proving or verification.
#[derive(thiserror::Error, Debug)]
pub enum SynthesisError {
    /// During synthesis, we lacked knowledge of a variable assignment.
    #[error("an assignment for a variable could not be computed")]
    AssignmentMissing,
    /// During synthesis, we divided by zero.
    #[error("division by zero")]
    DivisionByZero,
    /// During synthesis, we constructed an unsatisfiable constraint system.
    #[error("unsatisfiable constraint system")]
    Unsatisfiable,
    /// During synthesis, our polynomials ended up being too high of degree
    #[error("polynomial degree is too large")]
    PolynomialDegreeTooLarge,
    /// During proof generation, we encountered an identity in the CRS
    #[error("encountered an identity element in the CRS")]
    UnexpectedIdentity,
    /// During proof generation, we encountered an I/O error with the CRS
    #[error("encountered an I/O error: {0}")]
    IoError(#[from] io::Error),
    /// During verification, our verifying key was malformed.
    #[error("malformed verifying key")]
    MalformedVerifyingKey,
    /// During CRS generation, we observed an unconstrained auxiliary variable
    #[error("auxiliary variable was unconstrained")]
    UnconstrainedVariable,
    /// During GPU multiexp/fft, some GPU related error happened
    #[error("encountered a GPU error: {0}")]
    GPUError(#[from] gpu::error::GPUError),
}