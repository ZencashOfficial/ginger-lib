use crate::gpu::{GPUError, GPUResult, utils, kernel_multiexp};
use algebra::{
    AffineCurve, ProjectiveCurve, PairingEngine,
    PrimeField,
};
use log::{error, info};
use rust_gpu_tools::*;
use std::any::TypeId;

use super::get_cpu_utilization;

const MAX_WINDOW_SIZE: usize = 10;
const LOCAL_WORK_SIZE: usize = 256;
const MEMORY_PADDING: f64 = 0.2f64; // Let 20% of GPU memory be free
    
fn calc_num_groups(core_count: usize, num_windows: usize) -> usize {
    // Observations show that we get the best performance when num_groups * num_windows ~= 2 * CUDA_CORES
    2 * core_count / num_windows
}

fn calc_window_size(n: usize, scalar_bits: usize, core_count: usize) -> usize {
    // window_size = ln(n / num_groups)
    // num_windows = scalar_bits / window_size
    // num_groups = 2 * core_count / num_windows = 2 * core_count * window_size / scalar_bits
    // window_size = ln(n / num_groups) = ln(n * scalar_bits / (2 * core_count * window_size))
    // window_size = ln(scalar_bits * n / (2 * core_count)) - ln(window_size)
    //
    // Thus we need to solve the following equation:
    // window_size + ln(window_size) = ln(scalar_bits * n / (2 * core_count))
    let lower_bound = (((scalar_bits * n) as f64) / ((2 * core_count) as f64)).ln();
    for w in 0..MAX_WINDOW_SIZE {
        if (w as f64) + (w as f64).ln() > lower_bound {
            return w;
        }
    }

    MAX_WINDOW_SIZE
}

fn calc_best_chunk_size(max_window_size: usize, core_count: usize, scalar_bits: usize) -> usize {
    // Best chunk-size (N) can also be calculated using the same logic as calc_window_size:
    // n = e^window_size * window_size * 2 * core_count / scalar_bits
    (((max_window_size as f64).exp() as f64)
        * (max_window_size as f64)
        * 2f64
        * (core_count as f64)
        / (scalar_bits as f64))
        .ceil() as usize
}

fn calc_chunk_size<E>(mem: u64, core_count: usize) -> usize
where
    E: PairingEngine,
{
    let aff_size = std::mem::size_of::<E::G1Affine>() + std::mem::size_of::<E::G2Affine>();
    let scalar_size = std::mem::size_of::<<E::Fr as PrimeField>::BigInt>();
    let proj_size = std::mem::size_of::<E::G1Projective>() + std::mem::size_of::<E::G2Projective>();
    ((((mem as f64) * (1f64 - MEMORY_PADDING)) as usize)
        - (2 * core_count * ((1 << MAX_WINDOW_SIZE) + 1) * proj_size))
        / (aff_size + scalar_size)
}

// Multiscalar kernel for a single GPU
pub struct SingleMSMKernel<E>
where
    E: PairingEngine,
{
    pub program: opencl::Program,

    pub core_count: usize,
    pub n: usize,

    _phantom: std::marker::PhantomData<<E::Fr as PrimeField>::BigInt>,
}

impl<E> SingleMSMKernel<E>
where
    E: PairingEngine,
{
    pub fn create(d: opencl::Device) -> GPUResult<SingleMSMKernel<E>> {

        let src = kernel_multiexp::<E>(true);

        let scalar_bits = std::mem::size_of::<<E::Fr as PrimeField>::BigInt>() * 8;
        let core_count = utils::get_core_count(&d);
        let mem = d.memory();
        let max_n = calc_chunk_size::<E>(mem, core_count);
        let best_n = calc_best_chunk_size(MAX_WINDOW_SIZE, core_count, scalar_bits);
        let n = std::cmp::min(max_n, best_n);

        Ok(SingleMSMKernel {
            program: opencl::Program::from_opencl(d, &src)?,
            core_count,
            n,
            _phantom: std::marker::PhantomData,
        })
    }

    pub fn msm<G>(
        &self,
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        n: usize,
    ) -> GPUResult<G::Projective>
    where
        G: AffineCurve
    {
        let scalar_bits = std::mem::size_of::<<G::ScalarField as PrimeField>::BigInt>() * 8;
        let window_size = calc_window_size(n as usize, scalar_bits, self.core_count);
        let num_windows = ((scalar_bits as f64) / (window_size as f64)).ceil() as usize;
        let num_groups = calc_num_groups(self.core_count, num_windows);
        let bucket_len = 1 << window_size;

        // Each group will have `num_windows` threads and as there are `num_groups` groups, there will
        // be `num_groups` * `num_windows` threads in total.
        // Each thread will use `num_groups` * `num_windows` * `bucket_len` buckets.

        let mut base_buffer = self.program.create_buffer::<G>(n)?;
        base_buffer.write_from(0, bases)?;
        let mut scalar_buffer = self
            .program
            .create_buffer::<<G::ScalarField as PrimeField>::BigInt>(n)?;
        scalar_buffer.write_from(0, scalars)?;

        let bucket_buffer = self
            .program
            .create_buffer::<G::Projective>(2 * self.core_count * bucket_len)?;
        let result_buffer = self
            .program
            .create_buffer::<G::Projective>(2 * self.core_count)?;

        // Make global work size divisible by `LOCAL_WORK_SIZE`
        let mut global_work_size = num_windows * num_groups;
        global_work_size +=
            (LOCAL_WORK_SIZE - (global_work_size % LOCAL_WORK_SIZE)) % LOCAL_WORK_SIZE;

        let kernel = self.program.create_kernel(
            if TypeId::of::<G>() == TypeId::of::<E::G1Affine>() {
                "G1_bellman_multiexp"
            } else if TypeId::of::<G>() == TypeId::of::<E::G2Affine>() {
                "G2_bellman_multiexp"
            } else {
                return Err(GPUError::Simple("Only E::G1 and E::G2 are supported!"));
            },
            global_work_size,
            None,
        );

        call_kernel!(
            kernel,
            &base_buffer,
            &bucket_buffer,
            &result_buffer,
            &scalar_buffer,
            n as u32,
            num_groups as u32,
            num_windows as u32,
            window_size as u32
        )?;

        let mut results = vec![G::Projective::zero(); num_groups * num_windows];
        result_buffer.read_into(0, &mut results)?;

        // Using the algorithm below, we can calculate the final result by accumulating the results
        // of those `NUM_GROUPS` * `NUM_WINDOWS` threads.
        let mut acc = G::Projective::zero();
        let mut bits = 0;
        for i in 0..num_windows {
            let w = std::cmp::min(window_size, scalar_bits - bits);
            for _ in 0..w {
                acc.double_in_place();
            }
            for g in 0..num_groups {
                acc.add_assign_mixed(&results[g * num_windows + i].into_affine());
            }
            bits += w; // Process the next window
        }

        Ok(acc)
    }
}

pub struct MSMWorker<E>
where
    E: PairingEngine,
{
    kernels: Vec<SingleMSMKernel<E>>
}

impl<E> MSMWorker<E>
where
    E: PairingEngine,
{
    pub fn create() -> GPUResult<MSMWorker<E>> {

        let devices = opencl::Device::all()?;

        let kernels: Vec<_> = devices
            .into_iter()
            .map(|d| (d.clone(), SingleMSMKernel::<E>::create(d)))
            .filter_map(|(device, res)| {
                if let Err(ref e) = res {
                    error!(
                        "Cannot initialize kernel for device '{}'! Error: {}",
                        device.name(),
                        e
                    );
                }
                res.ok()
            })
            .collect();

        if kernels.is_empty() {
            return Err(GPUError::Simple("No working GPUs found!"));
        }
        info!(
            "Multiexp: {} working device(s) selected. (CPU utilization: {})",
            kernels.len(),
            get_cpu_utilization()
        );
        for (i, k) in kernels.iter().enumerate() {
            info!(
                "Multiexp: Device {}: {} (Chunk-size: {})",
                i,
                k.program.device().name(),
                k.n
            );
        }
        Ok(MSMWorker {
            kernels
        })
    }

    pub fn get_kernels(&self) -> &[SingleMSMKernel<E>]
    {
        self.kernels.as_slice()
    }
}
