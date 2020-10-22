use crate::gpu::error::{GPUError, GPUResult};
use crate::gpu::utils;
use crate::curves::bn_382::multicore::Worker;
use crate::curves::bn_382::multiexp::{multiexp as cpu_multiexp, FullDensity};
use super::locks;
use crossbeam::thread;
use futures::Future;
use log::{error, info};
use rust_gpu_tools::*;
use std::sync::Arc;

use algebra::{ProjectiveCurve, BigInteger384};
use algebra::curves::bn_382::{G1Affine, G1Projective, G2Affine, G2Projective};

const MAX_WINDOW_SIZE: usize = 10;
const LOCAL_WORK_SIZE: usize = 256;
const MEMORY_PADDING: f64 = 0.2f64; // Let 20% of GPU memory be free

pub fn get_cpu_utilization() -> f64 {
    use std::env;
    env::var("BELLMAN_CPU_UTILIZATION")
        .and_then(|v| match v.parse() {
            Ok(val) => Ok(val),
            Err(_) => {
                error!("Invalid BELLMAN_CPU_UTILIZATION! Defaulting to 0...");
                Ok(0f64)
            }
        })
        .unwrap_or(0f64)
        .max(0f64)
        .min(1f64)
}

// Multiexp kernel for a single GPU
pub struct SingleMultiexpKernel
{
    program: opencl::Program,

    core_count: usize,
    n: usize,

    priority: bool,
    _phantom: std::marker::PhantomData<BigInteger384>,
}

fn calc_num_groups(core_count: usize, num_windows: usize) -> usize {
    // Observations show that we get the best performance when num_groups * num_windows ~= 2 * CUDA_CORES
    2 * core_count / num_windows
}

fn calc_window_size(n: usize, exp_bits: usize, core_count: usize) -> usize {
    // window_size = ln(n / num_groups)
    // num_windows = exp_bits / window_size
    // num_groups = 2 * core_count / num_windows = 2 * core_count * window_size / exp_bits
    // window_size = ln(n / num_groups) = ln(n * exp_bits / (2 * core_count * window_size))
    // window_size = ln(exp_bits * n / (2 * core_count)) - ln(window_size)
    //
    // Thus we need to solve the following equation:
    // window_size + ln(window_size) = ln(exp_bits * n / (2 * core_count))
    let lower_bound = (((exp_bits * n) as f64) / ((2 * core_count) as f64)).ln();
    for w in 0..MAX_WINDOW_SIZE {
        if (w as f64) + (w as f64).ln() > lower_bound {
            return w;
        }
    }

    MAX_WINDOW_SIZE
}

fn calc_best_chunk_size(max_window_size: usize, core_count: usize, exp_bits: usize) -> usize {
    // Best chunk-size (N) can also be calculated using the same logic as calc_window_size:
    // n = e^window_size * window_size * 2 * core_count / exp_bits
    (((max_window_size as f64).exp() as f64)
        * (max_window_size as f64)
        * 2f64
        * (core_count as f64)
        / (exp_bits as f64))
        .ceil() as usize
}

fn calc_chunk_size(mem: u64, core_count: usize) -> usize
{
    let aff_size = std::mem::size_of::<G1Affine>() + std::mem::size_of::<G2Affine>();
    let exp_size = std::mem::size_of::<BigInteger384>();
    let proj_size = std::mem::size_of::<G1Projective>() + std::mem::size_of::<G2Projective>();
    ((((mem as f64) * (1f64 - MEMORY_PADDING)) as usize)
        - (2 * core_count * ((1 << MAX_WINDOW_SIZE) + 1) * proj_size))
        / (aff_size + exp_size)
}

impl SingleMultiexpKernel
{
    pub fn create(d: opencl::Device, priority: bool) -> GPUResult<SingleMultiexpKernel> {
        let src = super::sources::kernel(d.brand() == opencl::Brand::Nvidia);

        let exp_bits = std::mem::size_of::<BigInteger384>() * 8;
        let core_count = utils::get_core_count(&d);
        let mem = d.memory();
        let max_n = calc_chunk_size(mem, core_count);
        let best_n = calc_best_chunk_size(MAX_WINDOW_SIZE, core_count, exp_bits);
        let n = std::cmp::min(max_n, best_n);

        Ok(SingleMultiexpKernel {
            program: opencl::Program::from_opencl(d, &src)?,
            core_count,
            n,
            priority,
            _phantom: std::marker::PhantomData,
        })
    }

    pub fn multiexp_c(
        &mut self,
        bases: &[G1Affine],
        exps: &[BigInteger384],
        n: usize,
        c: usize
    ) -> GPUResult<G1Projective>
    {
        if locks::PriorityLock::should_break(self.priority) {
            return Err(GPUError::GPUTaken);
        }

        let exp_bits = std::mem::size_of::<BigInteger384>() * 8;
        //let window_size = calc_window_size(n as usize, exp_bits, self.core_count);
        let window_size = c;
        let num_windows = ((exp_bits as f64) / (window_size as f64)).ceil() as usize;
        let num_groups = calc_num_groups(self.core_count, num_windows);
        let bucket_len = 1 << window_size;

        // Each group will have `num_windows` threads and as there are `num_groups` groups, there will
        // be `num_groups` * `num_windows` threads in total.
        // Each thread will use `num_groups` * `num_windows` * `bucket_len` buckets.

        let mut base_buffer = self.program.create_buffer::<G1Affine>(n)?;
        base_buffer.write_from(bases)?;
        let mut exp_buffer = self
            .program
            .create_buffer::<BigInteger384>(n)?;
        exp_buffer.write_from(exps)?;

        let bucket_buffer = self
            .program
            .create_buffer::<G1Projective>(2 * self.core_count * bucket_len)?;
        let result_buffer = self
            .program
            .create_buffer::<G1Projective>(2 * self.core_count)?;

        // Make global work size divisible by `LOCAL_WORK_SIZE`
        let mut global_work_size = num_windows * num_groups;
        global_work_size +=
            (LOCAL_WORK_SIZE - (global_work_size % LOCAL_WORK_SIZE)) % LOCAL_WORK_SIZE;

        let kernel = self.program.create_kernel(
            "G1_bellman_multiexp",
            global_work_size,
            None,
        );

        call_kernel!(
            kernel,
            &base_buffer,
            &bucket_buffer,
            &result_buffer,
            &exp_buffer,
            n as u32,
            num_groups as u32,
            num_windows as u32,
            window_size as u32
        )?;

        let mut results = vec![G1Projective::zero(); num_groups * num_windows];
        result_buffer.read_into(&mut results)?;

        // Using the algorithm below, we can calculate the final result by accumulating the results
        // of those `NUM_GROUPS` * `NUM_WINDOWS` threads.
        let mut acc = G1Projective::zero();
        let mut bits = 0;
        for i in 0..num_windows {
            let w = std::cmp::min(window_size, exp_bits - bits);
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

    pub fn multiexp(
        &mut self,
        bases: &[G1Affine],
        exps: &[BigInteger384],
        n: usize,
    ) -> GPUResult<G1Projective>
    {
        if locks::PriorityLock::should_break(self.priority) {
            return Err(GPUError::GPUTaken);
        }

        let exp_bits = std::mem::size_of::<BigInteger384>() * 8;
        let window_size = calc_window_size(n as usize, exp_bits, self.core_count);
        let num_windows = ((exp_bits as f64) / (window_size as f64)).ceil() as usize;
        let num_groups = calc_num_groups(self.core_count, num_windows);
        let bucket_len = 1 << window_size;

        // Each group will have `num_windows` threads and as there are `num_groups` groups, there will
        // be `num_groups` * `num_windows` threads in total.
        // Each thread will use `num_groups` * `num_windows` * `bucket_len` buckets.

        let mut base_buffer = self.program.create_buffer::<G1Affine>(n)?;
        base_buffer.write_from(bases)?;
        let mut exp_buffer = self
            .program
            .create_buffer::<BigInteger384>(n)?;
        exp_buffer.write_from(exps)?;

        let bucket_buffer = self
            .program
            .create_buffer::<G1Projective>(2 * self.core_count * bucket_len)?;
        let result_buffer = self
            .program
            .create_buffer::<G1Projective>(2 * self.core_count)?;

        // Make global work size divisible by `LOCAL_WORK_SIZE`
        let mut global_work_size = num_windows * num_groups;
        global_work_size +=
            (LOCAL_WORK_SIZE - (global_work_size % LOCAL_WORK_SIZE)) % LOCAL_WORK_SIZE;

        let kernel = self.program.create_kernel(
            "G1_bellman_multiexp",
            global_work_size,
            None,
        );

        call_kernel!(
            kernel,
            &base_buffer,
            &bucket_buffer,
            &result_buffer,
            &exp_buffer,
            n as u32,
            num_groups as u32,
            num_windows as u32,
            window_size as u32
        )?;

        let mut results = vec![G1Projective::zero(); num_groups * num_windows];
        result_buffer.read_into(&mut results)?;

        // Using the algorithm below, we can calculate the final result by accumulating the results
        // of those `NUM_GROUPS` * `NUM_WINDOWS` threads.
        let mut acc = G1Projective::zero();
        let mut bits = 0;
        for i in 0..num_windows {
            let w = std::cmp::min(window_size, exp_bits - bits);
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

// A struct that containts several multiexp kernels for different devices
pub struct MultiexpKernel
{
    kernels: Vec<SingleMultiexpKernel>,
    _lock: locks::GPULock, // RFC 1857: struct fields are dropped in the same order as they are declared.
}

impl MultiexpKernel
{
    pub fn create(priority: bool) -> GPUResult<MultiexpKernel> {
        let lock = locks::GPULock::lock();

        let devices = opencl::Device::all()?;

        let kernels: Vec<_> = devices
            .into_iter()
            .map(|d| (d.clone(), SingleMultiexpKernel::create(d, priority)))
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
        Ok(MultiexpKernel {
            kernels,
            _lock: lock,
        })
    }

    pub fn multiexp_c(
        &mut self,
        pool: &Worker,
        bases: Arc<Vec<G1Affine>>,
        exps: Arc<Vec<BigInteger384>>,
        skip: usize,
        n: usize,
        c: usize
    ) -> GPUResult<G1Projective>
    {
        let num_devices = self.kernels.len();
        // Bases are skipped by `self.1` elements, when converted from (Arc<Vec<G>>, usize) to Source
        // https://github.com/zkcrypto/bellman/blob/10c5010fd9c2ca69442dc9775ea271e286e776d8/src/multiexp.rs#L38
        let bases = &bases[skip..(skip + n)];
        let exps = &exps[..n];

        let cpu_n = ((n as f64) * get_cpu_utilization()) as usize;
        let n = n - cpu_n;
        let (cpu_bases, bases) = bases.split_at(cpu_n);
        let (cpu_exps, exps) = exps.split_at(cpu_n);

        let chunk_size = ((n as f64) / (num_devices as f64)).ceil() as usize;

        match thread::scope(|s| -> Result<G1Projective, GPUError> {
            let mut acc = G1Projective::zero();
            let mut threads = Vec::new();
            if n > 0 {
                for ((bases, exps), kern) in bases
                    .chunks(chunk_size)
                    .zip(exps.chunks(chunk_size))
                    .zip(self.kernels.iter_mut())
                {
                    threads.push(s.spawn(
                        move |_| -> Result<G1Projective, GPUError> {
                            let mut acc = G1Projective::zero();
                            for (bases, exps) in bases.chunks(kern.n).zip(exps.chunks(kern.n)) {
                                let result = kern.multiexp_c(bases, exps, bases.len(), c)?;
                                acc.add_assign_mixed(&result.into_affine());
                            }
                            Ok(acc)
                        },
                    ));
                }
            }

            let cpu_acc = cpu_multiexp(
                &pool,
                (Arc::new(cpu_bases.to_vec()), 0),
                FullDensity,
                Arc::new(cpu_exps.to_vec()),
                &mut None,
            );

            let mut results = vec![];
            for t in threads {
                results.push(t.join());
            }
            for r in results {
                acc.add_assign_mixed(&r??.into_affine());
            }

            acc.add_assign_mixed(&cpu_acc.wait().unwrap().into_affine());

            Ok(acc)
        }) {
            Ok(res) => res,
            Err(e) => Err(GPUError::from(e)),
        }
    }


    pub fn multiexp(
        &mut self,
        pool: &Worker,
        bases: Arc<Vec<G1Affine>>,
        exps: Arc<Vec<BigInteger384>>,
        skip: usize,
        n: usize,
    ) -> GPUResult<G1Projective>
    {
        let num_devices = self.kernels.len();
        // Bases are skipped by `self.1` elements, when converted from (Arc<Vec<G>>, usize) to Source
        // https://github.com/zkcrypto/bellman/blob/10c5010fd9c2ca69442dc9775ea271e286e776d8/src/multiexp.rs#L38
        let bases = &bases[skip..(skip + n)];
        let exps = &exps[..n];

        let cpu_n = ((n as f64) * get_cpu_utilization()) as usize;
        let n = n - cpu_n;
        let (cpu_bases, bases) = bases.split_at(cpu_n);
        let (cpu_exps, exps) = exps.split_at(cpu_n);

        let chunk_size = ((n as f64) / (num_devices as f64)).ceil() as usize;

        match thread::scope(|s| -> Result<G1Projective, GPUError> {
            let mut acc = G1Projective::zero();
            let mut threads = Vec::new();
            if n > 0 {
                for ((bases, exps), kern) in bases
                    .chunks(chunk_size)
                    .zip(exps.chunks(chunk_size))
                    .zip(self.kernels.iter_mut())
                {
                    threads.push(s.spawn(
                        move |_| -> Result<G1Projective, GPUError> {
                            let mut acc = G1Projective::zero();
                            for (bases, exps) in bases.chunks(kern.n).zip(exps.chunks(kern.n)) {
                                let result = kern.multiexp(bases, exps, bases.len())?;
                                acc.add_assign_mixed(&result.into_affine());
                            }
                            Ok(acc)
                        },
                    ));
                }
            }

            let cpu_acc = cpu_multiexp(
                &pool,
                (Arc::new(cpu_bases.to_vec()), 0),
                FullDensity,
                Arc::new(cpu_exps.to_vec()),
                &mut None,
            );

            let mut results = vec![];
            for t in threads {
                results.push(t.join());
            }
            for r in results {
                acc.add_assign_mixed(&r??.into_affine());
            }

            acc.add_assign_mixed(&cpu_acc.wait().unwrap().into_affine());

            Ok(acc)
        }) {
            Ok(res) => res,
            Err(e) => Err(GPUError::from(e)),
        }
    }
}
