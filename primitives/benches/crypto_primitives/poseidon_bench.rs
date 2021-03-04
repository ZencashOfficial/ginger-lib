use criterion::{criterion_main, criterion_group, Criterion};

mod variable_base_msm_affine {

    use std::fs::File;
    use std::path::Path;

    use algebra::{
        curves::tweedle::dee::{Projective as G1Projective, Affine as G1Affine},
        BigInteger256, UniformRand, ProjectiveCurve, FromBytes, ToBytes
    };
    use algebra_utils::msm::VariableBaseMSM;
    use criterion::{BenchmarkId, BatchSize, Criterion};
    use rand::{SeedableRng};

    const DATA_PATH: &str = "/tmp/variable_msm_points.dat";
    const MAX_NUM_POINTS: usize = 2097152; //2^21
    const BASE_POINTS_SEED: u64 = 4001728647;
    const COEFFICIENTS_SEED: u64 = 1299915800;

    fn generate_base_points(max_num_elements: usize) {
        let mut random_generator = rand_chacha::ChaChaRng::seed_from_u64(BASE_POINTS_SEED);
        let mut file = File::create(DATA_PATH).unwrap();

        let start = std::time::Instant::now();

        for _ in 0..max_num_elements {
            let element: G1Affine = G1Projective::rand(&mut random_generator).into_affine();

            match element.write(&mut file) {
                Ok(_) => {},
                Err(msg) => { panic!("Cannot save base points to file: {}", msg)}
            }
        }

        println!("Vector generation time: {:?}", start.elapsed());
    }

    fn generate_coefficients(num_elements: usize) -> Vec<BigInteger256> {
        let mut random_generator = rand_chacha::ChaChaRng::seed_from_u64(COEFFICIENTS_SEED);
        let mut coefficients: Vec<BigInteger256> = Vec::with_capacity(num_elements);

        for _ in 0..num_elements {
            coefficients.push(BigInteger256::rand(&mut random_generator));
        }

        coefficients
    }

    fn load_base_points(num_elements: usize) -> Vec<G1Affine> {
        if !Path::new(DATA_PATH).exists() {
            generate_base_points(MAX_NUM_POINTS);
        }

        let mut fs = File::open(DATA_PATH).unwrap();
        let mut base_points: Vec<G1Affine> = Vec::with_capacity(num_elements);

        for _ in 0..num_elements {
            base_points.push(G1Affine::read(&mut fs).unwrap());
        }

        base_points
    }

    pub fn benchmark(c: &mut Criterion) {
        let mut group = c.benchmark_group("Variable base MSM affine, size 2^");
        let base_points = load_base_points(MAX_NUM_POINTS);
        let num_scalars_pow = (12..=21).collect::<Vec<_>>();

        for pow in num_scalars_pow {
            group.bench_with_input(BenchmarkId::from_parameter(pow), &pow, |b, pow| {
                b.iter_batched(|| {
                    let coefficients = generate_coefficients(2usize.pow(*pow));
    
                    (&base_points, coefficients)
                },
                |(base_points, coefficients)| {
                    VariableBaseMSM::multi_scalar_mul(&base_points[0..2usize.pow(*pow)], coefficients.as_slice());
                },
                BatchSize::PerIteration);
            });
        }
    }
}

criterion_group!(
name = poseidon_benchmark;
config = Criterion::default().sample_size(50);
targets = variable_base_msm_affine::benchmark
);

criterion_main!(poseidon_benchmark);