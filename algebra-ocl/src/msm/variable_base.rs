use crate::gpu::GPUError;
use algebra::{
    AffineCurve, ProjectiveCurve, PairingEngine,
    PrimeField,
};
use crossbeam::thread;
use super::get_cpu_utilization;
use super::MSMWorker;

pub struct VariableBaseMSM;

impl VariableBaseMSM
{
    fn msm_inner<G, E>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        msm_worker: &MSMWorker<E>
    ) -> G::Projective
    where
        G: AffineCurve,
        E: PairingEngine,
        G::Projective: ProjectiveCurve<Affine = G>
    {
        let zero = G::Projective::zero();

        let n = bases.len();
        let num_devices = msm_worker.get_kernels().len();
        
        let cpu_n = ((n as f64) * get_cpu_utilization()) as usize;
        let n = n - cpu_n;
        let (cpu_bases, bases) = bases.split_at(cpu_n);
        let (cpu_scalars, scalars) = scalars.split_at(cpu_n);

        let chunk_size = ((n as f64) / (num_devices as f64)).ceil() as usize;

        match thread::scope(|s| -> Result<G::Projective, GPUError> {
            let mut acc = G::Projective::zero();
            let mut threads = Vec::new();

            if n > 0 {
                for ((bases, scalars), kern) in bases
                    .chunks(chunk_size)
                    .zip(scalars.chunks(chunk_size))
                    .zip(msm_worker.get_kernels().iter())
                {
                    threads.push(s.spawn(
                        move |_| -> Result<G::Projective, GPUError> {
                            let mut acc = G::Projective::zero();
                            for (bases, scalars) in bases.chunks(kern.n).zip(scalars.chunks(kern.n)) {
                                let result = kern.msm(bases, scalars, bases.len())?;
                                acc.add_assign_mixed(&result.into_affine());
                            }
                            Ok(acc)
                        },
                    ));
                }
            }

            if cpu_n > 0 {
                threads.push(s.spawn(
                    move |_| -> Result<G::Projective, GPUError> {
                        let acc = algebra::msm::VariableBaseMSM::multi_scalar_mul(cpu_bases, cpu_scalars);
                        Ok(acc)
                    }
                ))
            }

            let mut results = vec![];
            for t in threads {
                results.push(t.join());
            }
            for r in results {
                acc.add_assign_mixed(&r??.into_affine());
            }

            Ok(acc)
        }) {
            Ok(res) => res.unwrap(),
            Err(_) => zero
        }
    }

    pub fn multi_scalar_mul<G, E>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        msm_worker: &MSMWorker<E>
    ) -> G::Projective 
    where
        G: AffineCurve,
        E: PairingEngine,
        G::Projective: ProjectiveCurve<Affine = G>,
    {
        Self::msm_inner(bases, scalars, msm_worker)
    }
}
