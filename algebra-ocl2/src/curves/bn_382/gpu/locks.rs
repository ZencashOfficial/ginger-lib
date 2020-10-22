use fs2::FileExt;
use log::{debug, info, warn};
use std::fs::File;
use std::path::PathBuf;

const GPU_LOCK_NAME: &str = "bellman.gpu.lock";
const PRIORITY_LOCK_NAME: &str = "bellman.priority.lock";
fn tmp_path(filename: &str) -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(filename);
    p
}

/// `GPULock` prevents two kernel objects to be instantiated simultaneously.
#[derive(Debug)]
pub struct GPULock(File);
impl GPULock {
    pub fn lock() -> GPULock {
        debug!("Acquiring GPU lock...");
        let f = File::create(tmp_path(GPU_LOCK_NAME)).unwrap();
        f.lock_exclusive().unwrap();
        debug!("GPU lock acquired!");
        GPULock(f)
    }
}
impl Drop for GPULock {
    fn drop(&mut self) {
        debug!("GPU lock released!");
    }
}

/// `PrioriyLock` is like a flag. When acquired, it means a high-priority process
/// needs to acquire the GPU really soon. Acquiring the `PriorityLock` is like
/// signaling all other processes to release their `GPULock`s.
/// Only one process can have the `PriorityLock` at a time.
#[derive(Debug)]
pub struct PriorityLock(File);
impl PriorityLock {
    pub fn lock() -> PriorityLock {
        debug!("Acquiring priority lock...");
        let f = File::create(tmp_path(PRIORITY_LOCK_NAME)).unwrap();
        f.lock_exclusive().unwrap();
        debug!("Priority lock acquired!");
        PriorityLock(f)
    }
    pub fn wait(priority: bool) {
        if !priority {
            File::create(tmp_path(PRIORITY_LOCK_NAME))
                .unwrap()
                .lock_exclusive()
                .unwrap();
        }
    }
    pub fn should_break(priority: bool) -> bool {
        !priority
            && File::create(tmp_path(PRIORITY_LOCK_NAME))
                .unwrap()
                .try_lock_exclusive()
                .is_err()
    }
}
impl Drop for PriorityLock {
    fn drop(&mut self) {
        debug!("Priority lock released!");
    }
}

use crate::gpu::error::{GPUError, GPUResult};
use crate::curves::bn_382::multiexp::create_multiexp_kernel;
use super::multiexp::MultiexpKernel;

macro_rules! locked_kernel {
    ($class:ident, $kern:ident, $func:ident, $name:expr) => {
        pub struct $class
        {
            log_d: usize,
            priority: bool,
            kernel: Option<$kern>,
        }

        impl $class
        {
            pub fn new(log_d: usize, priority: bool) -> $class {
                $class {
                    log_d,
                    priority,
                    kernel: None,
                }
            }

            fn init(&mut self) {
                if self.kernel.is_none() {
                    PriorityLock::wait(self.priority);
                    info!("GPU is available for {}!", $name);
                    self.kernel = $func(self.log_d, self.priority);
                }
            }

            fn free(&mut self) {
                if let Some(_kernel) = self.kernel.take() {
                    warn!(
                        "GPU acquired by a high priority process! Freeing up {} kernels...",
                        $name
                    );
                }
            }

            pub fn with<F, R>(&mut self, mut f: F) -> GPUResult<R>
            where
                F: FnMut(&mut $kern) -> GPUResult<R>,
            {
                if std::env::var("BELLMAN_NO_GPU").is_ok() {
                    return Err(GPUError::GPUDisabled);
                }

                self.init();

                loop {
                    if let Some(ref mut k) = self.kernel {
                        match f(k) {
                            Err(GPUError::GPUTaken) => {
                                self.free();
                                self.init();
                            }
                            Err(e) => {
                                warn!("GPU {} failed! Falling back to CPU... Error: {}", $name, e);
                                return Err(e);
                            }
                            Ok(v) => return Ok(v),
                        }
                    } else {
                        return Err(GPUError::KernelUninitialized);
                    }
                }
            }
        }
    };
}

locked_kernel!(
    LockedMultiexpKernel,
    MultiexpKernel,
    create_multiexp_kernel,
    "Multiexp"
);
