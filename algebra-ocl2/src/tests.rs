use std::sync::Arc;
use std::time::Instant;
//use algebra_core::{UniformRand, ProjectiveCurve, PrimeField, BigInteger384, FromBytes};
use algebra::{BigInteger384, FromBytes};
use algebra::curves::bn_382::G1Affine;
use std::fs::File;

#[test]
fn gpu_bn_382_source_generate() {
    use crate::curves::bn_382::gpu::kernel;

    let src = kernel(true);

    println!("{}", src);
}

#[test]
fn gpu_bn_382_multiexp_test() {
    //use algebra::bn_382::{G1Projective, Fp};
    use futures::Future;
    use crate::curves::bn_382::{
        multiexp::{
            multiexp,
            FullDensity
        },
        multicore::Worker,
        gpu::LockedMultiexpKernel
    };

    const START_LOG_D: usize = 23;
    const MAX_LOG_D: usize = 23;

    //let mut rng = &mut rand::thread_rng();
    let mut kern = Some(LockedMultiexpKernel::new(MAX_LOG_D, false));
    let pool = Worker::new();

    // println!("Start bases generation...");
    //
    // let mut bases = (0..(1 << 10))
    //     .map(|_| G1Projective::rand(&mut rng).into_affine())
    //     .collect::<Vec<_>>();
    // for _ in 10..START_LOG_D {
    //     bases = [bases.clone(), bases.clone()].concat();
    // }
    //
    // println!("Bases are generated!");

    for log_d in START_LOG_D..(MAX_LOG_D + 1) {

        let samples = 1 << log_d;
        let (v1, g1) = load_data(samples);

        for c in 3..14 {

            println!("Testing Multiexp for {} elements... c={}", samples, c);

            //let (v1, g1) = load_data(samples);

            // let scalars = (0..samples)
            //     .map(|_| Fp::rand(&mut rng).into_repr())
            //     .collect::<Vec<_>>();

            let g = Arc::new(g1.clone());
            let v = Arc::new(v1.clone());

            let now = Instant::now();
            //let gpu = multiexp_c(&pool, (g.clone(), 0), FullDensity, v.clone(), &mut kern, c)
            let gpu = multiexp(&pool, (g.clone(), 0), FullDensity, v.clone(), &mut kern)
                .wait()
                .unwrap();
            let gpu_dur = now.elapsed().as_secs() * 1000 as u64 + now.elapsed().subsec_millis() as u64;
            println!("GPU took {}ms.", gpu_dur);

            // let now = Instant::now();
            // let cpu = algebra_core::msm::VariableBaseMSM::multi_scalar_mul_affine_sd(g1.as_slice(), v1.as_slice());
            // let cpu_dur = now.elapsed().as_secs() * 1000 as u64 + now.elapsed().subsec_millis() as u64;
            // println!("CPU took {}ms.", cpu_dur);
            //
            // println!("Speedup: x{}", cpu_dur as f32 / gpu_dur as f32);
            //
            // assert_eq!(cpu, gpu);

            println!("============================");
        }
       // bases = [bases.clone(), bases.clone()].concat();
    }
}

//const SAMPLES: usize = 1<<23;

fn load_data(samples:i64) -> (Vec<BigInteger384>,Vec<G1Affine>) {

    let mut fs = File::open("/home/mkaihara/Documents/dev/HL_Zexe_Coda/algebra-ocl2/src/scalars_bases").unwrap();
    let mut v = Vec::with_capacity(samples as usize);
    let mut g = Vec::with_capacity(samples as usize);

    for _i in 0..samples {
        let elem = BigInteger384::read(&mut fs).unwrap();
        v.push(elem);
        let elem = G1Affine::read(&mut fs).unwrap();
        g.push(elem);
    }
    (v,g)
}
