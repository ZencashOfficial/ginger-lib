mod variable_base_msm_affine {

    use std::fs::File;
    use std::path::Path;

    use algebra::{
        curves::tweedle::dee::{Projective as G1Projective, Affine as G1Affine},
        BigInteger256, UniformRand, ProjectiveCurve, FromBytes, ToBytes
    };
    use rand::{SeedableRng};

    pub const MAX_NUM_POINTS: usize = 1 << 21;

    const DATA_PATH: &str = "/tmp/variable_msm_points.dat";
    const BASE_POINTS_SEED: u64 = 4001728647;
    const COEFFICIENTS_SEED: u64 = 1299915800;

    pub fn generate_base_points(max_num_elements: usize) {
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

    pub fn generate_coefficients(num_elements: usize) -> Vec<BigInteger256> {
        let mut random_generator = rand_chacha::ChaChaRng::seed_from_u64(COEFFICIENTS_SEED);
        let mut coefficients: Vec<BigInteger256> = Vec::with_capacity(num_elements);

        for _ in 0..num_elements {
            coefficients.push(BigInteger256::rand(&mut random_generator));
        }

        coefficients
    }

    pub fn load_base_points(num_elements: usize) -> Vec<G1Affine> {
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
}

mod bowe_hopwood {

    use primitives::crh::pedersen::PedersenWindow;

    pub const BASE_POINTS_SEED: u64 = 2411988112;
    pub const COEFFICIENTS_SEED: u64 = 2959228773;

    #[derive(Clone)]
    pub struct BenchmarkWindow {}

    impl PedersenWindow for BenchmarkWindow {
        const WINDOW_SIZE: usize = 64;
        const NUM_WINDOWS: usize = 21846; // to serve a bitvector of max size 2^22 bits.
    }
}

mod poseidon {

    use rand::Rng;

    use algebra::{fields::tweedle::Fr as TweedleFr, biginteger::BigInteger256, field_new, FromBits};
    use primitives::{
        crh::poseidon::parameters::tweedle::{TweedleFrPoseidonHash, TweedleFrBatchPoseidonHash},
        merkle_tree::field_based_mht::{
            FieldBasedMerkleTreeParameters, FieldBasedMerkleTreePrecomputedEmptyConstants,
            BatchFieldBasedMerkleTreeParameters, FieldBasedOptimizedMHT
        }
    };

    pub type TweedlePoseidonMHT = FieldBasedOptimizedMHT<TweedleFieldBasedMerkleTreeParams>;

    pub const TWEEDLE_MHT_POSEIDON_PARAMETERS: FieldBasedMerkleTreePrecomputedEmptyConstants<'static, TweedleFrPoseidonHash> =
    FieldBasedMerkleTreePrecomputedEmptyConstants {
        nodes: &[
            field_new!(TweedleFr, BigInteger256([0, 0, 0, 0])),
            field_new!(TweedleFr, BigInteger256([1685121506851839463, 16981548520505137274, 8913412180966583038, 3862086594483176853])),
            field_new!(TweedleFr, BigInteger256([11170234764586973796, 10725839658946531096, 4089526305768698754, 1460429405770838413])),
            field_new!(TweedleFr, BigInteger256([17916607642221592882, 8899562254223372899, 5718113766989047716, 1076406490778760958])),
            field_new!(TweedleFr, BigInteger256([16814393994791617149, 4003549585878224960, 91189797674671164, 2489268348396502646])),
            field_new!(TweedleFr, BigInteger256([3460562161953312097, 9272175926981975155, 3088781530878338181, 408762626573460268])),
            field_new!(TweedleFr, BigInteger256([4239231910949981686, 8869286273077200517, 15652728562557699831, 3576102463490196815])),
            field_new!(TweedleFr, BigInteger256([13380559706935466086, 9731495717823444557, 13485044794068855700, 920123060019082593])),
            field_new!(TweedleFr, BigInteger256([11296480562263487099, 13017529970429918980, 15492229384387414211, 1785749344954309826])),
            field_new!(TweedleFr, BigInteger256([9738206838161957894, 11701328267258081479, 11712994190976813369, 484156460220120380])),
            field_new!(TweedleFr, BigInteger256([3442436595719325694, 2431894943368234580, 17363188394577578270, 1614283916602080385])),
            field_new!(TweedleFr, BigInteger256([3972068898698943381, 11117259520378692208, 4253450360934634403, 239101287040501995])),
            field_new!(TweedleFr, BigInteger256([16580610996142154741, 7643908918402835608, 4538907510080444171, 4263606647859070795])),
            field_new!(TweedleFr, BigInteger256([8592329819226182315, 13382013060597509488, 8836783333517823500, 3786991317438607835])),
            field_new!(TweedleFr, BigInteger256([1773237737270872796, 2410057326158883883, 2737246428764964492, 143982332310609966])),
            field_new!(TweedleFr, BigInteger256([2707616787051394019, 11722910454912726756, 3333932516671235313, 2519035104169477350])),
            field_new!(TweedleFr, BigInteger256([5813486669979523104, 622781981755528836, 3795306211510301966, 2280896473682328864])),
            field_new!(TweedleFr, BigInteger256([9884651142463613651, 11192489360191969729, 3533957099232207802, 4365455770905630450])),
            field_new!(TweedleFr, BigInteger256([11940960536276443993, 12237559331804864182, 855310162320223754, 3250995830672065852])),
            field_new!(TweedleFr, BigInteger256([9844345352101824289, 11007606007155221061, 5337768567583173935, 1757712574894454062])),
            field_new!(TweedleFr, BigInteger256([17397097032371215469, 2560985447691462466, 5873153672621324816, 3004853082856123488])),
            field_new!(TweedleFr, BigInteger256([2906639260900593014, 12565622307387416057, 18208249875710924485, 3076942007405729816])),
            field_new!(TweedleFr, BigInteger256([12332293864653714593, 6336761349198557987, 14481492470913974358, 4413622281835599417])),
            field_new!(TweedleFr, BigInteger256([79692058390698902, 14504666139340959355, 16768104803024975328, 369341246831639552])),
            field_new!(TweedleFr, BigInteger256([2621547940131718999, 17743506990112205174, 6971752532124278629, 3241396862924145590])),
            field_new!(TweedleFr, BigInteger256([11163086857151067005, 12984510300292820262, 4717994806295806048, 4495298772085010977])),
            field_new!(TweedleFr, BigInteger256([15531065475439967178, 10139480523449757505, 6054301480319923629, 683911869248719980])),
            field_new!(TweedleFr, BigInteger256([3902755205937944407, 13093206721129864441, 16246765226559434352, 2735622140097895633])),
            field_new!(TweedleFr, BigInteger256([8690158721593369984, 730875283828647497, 3841607634391808388, 2371594890750177334])),
            field_new!(TweedleFr, BigInteger256([3026815083929086599, 14419444751560137298, 3651949377921605765, 3671437283458210474])),
            field_new!(TweedleFr, BigInteger256([14126118267082632123, 6827034892773261391, 5556113879686593587, 4539405425966083639])),
            field_new!(TweedleFr, BigInteger256([11865934732736629922, 2696435346237195560, 18058422909608662479, 2071268813115074451])),
            field_new!(TweedleFr, BigInteger256([13347092608138479115, 4066036543997219427, 7959630727080163169, 3601068210152443865])),
            field_new!(TweedleFr, BigInteger256([15138377926095811795, 15380039026895908568, 271121507124550600, 1024255194234800819])),
            field_new!(TweedleFr, BigInteger256([10134471144513638072, 16831793127198349974, 6339097251062257493, 541654095494825270])),
            field_new!(TweedleFr, BigInteger256([15800878851654609259, 9774427126525231605, 14726151358032257215, 2809506730656175605])),
            field_new!(TweedleFr, BigInteger256([12250046100357603655, 3019385177324988069, 844813400383977349, 2799285112129997345])),
            field_new!(TweedleFr, BigInteger256([5411213645516041177, 16651961682340139834, 2931464647892840002, 4405490370974890431])),
            field_new!(TweedleFr, BigInteger256([18074358143639506426, 14767712014724128887, 14587381883459330445, 4481815717901469638])),
            field_new!(TweedleFr, BigInteger256([3551361605813546439, 13978661106348185458, 11925279340031940548, 2940381839804134085]))
        ],
        merkle_arity: 2,
        max_height: 40,
    };

    pub const SEED: u64 = 406518596;

    const MAX_USABLE_BITS: usize = 253;

    #[derive(Clone, Debug)]
    pub struct TweedleFieldBasedMerkleTreeParams;

    impl FieldBasedMerkleTreeParameters for TweedleFieldBasedMerkleTreeParams {
        type Data = TweedleFr;
        type H = TweedleFrPoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> = Some(TWEEDLE_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for TweedleFieldBasedMerkleTreeParams {
        type BH = TweedleFrBatchPoseidonHash;
    }

    fn generate_random_field_element(size: usize, random_generator: & mut rand_chacha::ChaChaRng) -> Vec<bool> {

        let mut bit_vector: Vec<bool> = vec!(false; size);
        
        let distr = rand::distributions::Uniform::new_inclusive(0, 1);

        for i in 0..bit_vector.len() as usize {
            if random_generator.sample(distr) == 1 {
                bit_vector[i] = true;
            }
        }

        bit_vector
    }
    
    pub fn generate_random_leaves(size: usize, random_generator: & mut rand_chacha::ChaChaRng) -> Vec<TweedleFr> {
        let mut leaves: Vec<TweedleFr> = Vec::with_capacity(size);

        for _ in 0..size {
            let bit_vector = generate_random_field_element(MAX_USABLE_BITS, random_generator);
            leaves.push(TweedleFr::read_bits(bit_vector).unwrap());
        }

        leaves
    }
}

use criterion::{criterion_main, criterion_group, BenchmarkId, BatchSize, Criterion};
use rand::{Rng, SeedableRng};

use algebra::curves::tweedle::dee::Projective as TweedleDeeProjective;
use algebra_utils::msm::VariableBaseMSM;
use primitives::{
    crh::{bowe_hopwood::{BoweHopwoodPedersenCRH, BoweHopwoodPedersenParameters}},
    FixedLengthCRH, merkle_tree::field_based_mht::FieldBasedMerkleTree
};

pub fn benchmark(c: &mut Criterion) {

    // MSM AFFINE VARIABLES
    let mut group = c.benchmark_group("Variable base MSM affine, size 2^");
    let base_points = variable_base_msm_affine::load_base_points(variable_base_msm_affine::MAX_NUM_POINTS);
    let num_scalars_pow = (12..=21).collect::<Vec<_>>();

    // BOWE HOPWOOD VARIABLES
    let mut bowe_hopwood_random_base_points_generator = rand_chacha::ChaChaRng::seed_from_u64(bowe_hopwood::BASE_POINTS_SEED);
    let mut bowe_hopwood_random_coefficients_generator = rand_chacha::ChaChaRng::seed_from_u64(bowe_hopwood::COEFFICIENTS_SEED);
    let bowe_hopwood_params: BoweHopwoodPedersenParameters<TweedleDeeProjective> = <BoweHopwoodPedersenCRH<TweedleDeeProjective, bowe_hopwood::BenchmarkWindow> as FixedLengthCRH>::setup(&mut bowe_hopwood_random_base_points_generator).unwrap();

    // POSEIDON VARIABLES
    let mut poseidon_random_generator = rand_chacha::ChaChaRng::seed_from_u64(poseidon::SEED);

    for pow in num_scalars_pow {

        // MSM AFFINE BENCHMARK
        group.bench_with_input(BenchmarkId::new("MSM AFFINE - 2^", pow), &pow, |b, pow| {
            b.iter_batched(|| {
                let coefficients = variable_base_msm_affine::generate_coefficients(1 << *pow);

                (&base_points, coefficients)
            },
            |(base_points, coefficients)| {
                VariableBaseMSM::multi_scalar_mul(&base_points[0..1 << *pow], coefficients.as_slice());
            },
            BatchSize::PerIteration);
        });

        // BOWE HOPWOOD BENCHMARK
        group.bench_with_input(BenchmarkId::new("BOWE HOPWOOD - Tweedle FE 2^", pow), &pow, |b, pow| {
            b.iter_batched(|| {

                let size = (1 << *pow) / 8;
                let mut input: Vec<u8> = Vec::with_capacity(size);

                for _ in 0..size {
                    let bowe_hopwood_distribution = rand::distributions::Uniform::new(0, 1 << 8);
                    input.push(bowe_hopwood_random_coefficients_generator.sample(bowe_hopwood_distribution) as u8);
                }

                input
            },
            |input| {
                <BoweHopwoodPedersenCRH<TweedleDeeProjective, bowe_hopwood::BenchmarkWindow> as FixedLengthCRH>::evaluate(&bowe_hopwood_params, input.as_slice()).unwrap();
            },
            BatchSize::PerIteration);
        });

        // POSEIDON BENCHMARK
        group.bench_with_input(BenchmarkId::new("POSEIDON - 2^", pow), &pow, |b, pow| {
            b.iter_batched(|| {

                let num_leaves = 1 << *pow;
                let mt = poseidon::TweedlePoseidonMHT::init(
                    *pow as usize,
                    num_leaves,
                );

                let leaves = poseidon::generate_random_leaves(num_leaves, &mut poseidon_random_generator);

                (mt, leaves)
            },
            |(mut mt, leaves)| {
                leaves.iter().for_each(|&leaf| { mt.append(leaf); });
                mt.finalize_in_place().root().unwrap();
            },
            BatchSize::PerIteration);
        });
    }
}

criterion_group!(
name = poseidon_benchmark;
config = Criterion::default().sample_size(50);
targets = benchmark
);

criterion_main!(poseidon_benchmark);