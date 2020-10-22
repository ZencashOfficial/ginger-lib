/// Code adapted from https://github.com/filecoin-project/bellman

#[cfg(feature = "gpu")]
pub mod gpu;

#[cfg(feature = "gpu")]
pub mod msm;

#[cfg(all(test, features = "gpu"))]
mod tests;