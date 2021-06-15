use algebra::{AffineCurve, ToConstraintField, serialize::*};
use poly_commit::{
    PolynomialCommitment,
    ipa_pc::InnerProductArgPC
};
use proof_systems::darlin::{
    FinalDarlin,
    tests::{
        get_keys,
        final_darlin::generate_test_data as generate_final_darlin_test_data
    },
};
use digest::Digest;
use criterion::*;
use rand::{thread_rng, SeedableRng};
use blake2::Blake2s;
use rand_xorshift::XorShiftRng;

fn bench_single_verifier<G1: AffineCurve, G2: AffineCurve, D: Digest>(
    c: &mut Criterion,
    bench_name: &str,
    num_constraints: usize,
    segment_size_vals: Vec<usize>,
)
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let rng = &mut XorShiftRng::seed_from_u64(1234567890u64);
    let mut group = c.benchmark_group(bench_name);

    for segment_size in segment_size_vals.into_iter() {

        //Generate DLOG keys
        let params_g1 = InnerProductArgPC::<G1, D>::setup(segment_size - 1).unwrap();
        let params_g2 = InnerProductArgPC::<G2, D>::setup(segment_size - 1).unwrap();

        let (
            _, verifier_key_g1,
            _, verifier_key_g2
        ) = get_keys::<_, _, D>(&params_g1, &params_g2);

        let (final_darlin_pcd, index_vk) = generate_final_darlin_test_data::<G1, G2, D, _>(
            num_constraints,
            segment_size,
            &params_g1,
            &params_g2,
            1,
            rng
        );

        let (proof, vk, usr_ins) = (&final_darlin_pcd[0].final_darlin_proof, &index_vk[0], &final_darlin_pcd[0].usr_ins);

        {
            let proof_serialized_size = proof.serialized_size();
            let mut proof_serialized = Vec::with_capacity(proof_serialized_size);
            CanonicalSerialize::serialize(proof, &mut proof_serialized).unwrap();
            assert_eq!(proof_serialized.len(), proof_serialized_size);
            println!("Proof size: {:?}", proof_serialized_size);

            let vk_serialized_size = vk.serialized_size();
            let mut vk_serialized = Vec::with_capacity(vk_serialized_size);
            CanonicalSerialize::serialize(&vk, &mut vk_serialized).unwrap();
            assert_eq!(vk_serialized.len(), vk_serialized_size);
            println!("Vk size: {:?}", vk_serialized_size);
        }

        group.bench_with_input(BenchmarkId::from_parameter(segment_size), &segment_size, |bn, _segment_size| {
            bn.iter(|| {
                assert!(FinalDarlin::<G1, G2, D>::verify(
                    vk,
                    &verifier_key_g1,
                    &verifier_key_g2,
                    usr_ins,
                    proof,
                    &mut thread_rng()
                ).unwrap());
            });
        });
    }

    group.finish();
}

fn bench_single_verifier_ss_tweedle(c: &mut Criterion) {

    use algebra::curves::tweedle::{
        dee::Affine as TweedleDee,
        dum::Affine as TweedleDum,
    };

    /*bench_single_verifier::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = 1 << 19, |K| = 1 << 20, segment_size: ",
        1 << 19,
        vec![1 << 14, 1 << 15, 1 << 16, 1 << 17, 1 << 18],
    );*/

    bench_single_verifier::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = 1 << 20, |K| = 1 << 21, segment_size",
        1 << 20,
        vec![1 << 15, 1 << 16, 1 << 17, 1 << 18],
    );
}

criterion_group!(
name = single_verification_ss;
config = Criterion::default().sample_size(10);
targets = bench_single_verifier_ss_tweedle
);

criterion_main!(single_verification_ss);