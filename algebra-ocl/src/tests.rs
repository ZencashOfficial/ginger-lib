#[test]
fn test_algebra_ocl_gpu_kernel() {
    use crate::gpu::kernel_multiexp;
    use algebra::curves::bn_382::Bn382;

    let _ = env_logger::try_init();
    
    println!("{}", kernel_multiexp::<Bn382>(true));
}

#[test]
fn gpu_bn_382_msm_test() {
    use std::time::Instant;
    use algebra::{
        fields::bn_382::Fr,
        curves::bn_382::{BN382, G1Projective},
        UniformRand, ProjectiveCurve, PrimeField
    };
    use crate::msm::MSMWorker;

    let _ = env_logger::try_init();

    const START_LOG_D: usize = 10;
    const MAX_LOG_D: usize = 23;

    let mut rng = &mut rand::thread_rng();
    let msm_worker = MSMWorker::<Bn382>::create().unwrap();

    println!("Start bases generation...");

    let mut bases = (0..(1 << 10))
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();
    for _ in 10..START_LOG_D {
        bases = [bases.clone(), bases.clone()].concat();
    }

    println!("Bases are generated!");

    for log_d in START_LOG_D..(MAX_LOG_D + 1) {

        let samples = 1 << log_d;
        println!("Testing Multiexp for {} elements...", samples);

        let scalars = (0..samples)
            .map(|_| Fp::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();

        let now = Instant::now();
        let gpu = crate::msm::VariableBaseMSM::multi_scalar_mul(bases.as_slice(), scalars.as_slice(), &msm_worker);
        let gpu_dur = now.elapsed().as_secs() * 1000 as u64 + now.elapsed().subsec_millis() as u64;
        println!("GPU took {}ms.", gpu_dur);

        let now = Instant::now();
        let cpu = algebra_core::msm::VariableBaseMSM::multi_scalar_mul(bases.as_slice(), scalars.as_slice());
        let cpu_dur = now.elapsed().as_secs() * 1000 as u64 + now.elapsed().subsec_millis() as u64;
        println!("CPU took {}ms.", cpu_dur);

        println!("Speedup: x{}", cpu_dur as f32 / gpu_dur as f32);

        assert_eq!(cpu, gpu);

        println!("============================");

        bases = [bases.clone(), bases.clone()].concat();
    }
}
