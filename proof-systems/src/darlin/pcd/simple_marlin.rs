//! Simple Marlin "proof carrying data". This corresponds to non-recursive applications.
use algebra::{AffineCurve, SemanticallyValid, serialize::*};
use digest::Digest;
use marlin::{VerifierKey as MarlinVerifierKey, Proof, Marlin, AHPForR1CS};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, VerifierKey as DLogVerifierKey
    },
    rng::FiatShamirRng,
};
use crate::darlin::{
    pcd::{PCD, error::PCDError},
    accumulators::{
        dlog::{DLogItem, DLogItemAccumulator}, ItemAccumulator
    },
};
use poly_commit::ipa_pc::Commitment;
use std::ops::Deref;
use std::marker::PhantomData;

#[derive(Derivative)]
#[derivative(Clone(bound = ""), Debug(bound = ""), Eq(bound = ""), PartialEq(bound = ""))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct MarlinProof<G: AffineCurve, D: Digest>(pub Proof<G::ScalarField, InnerProductArgPC<G, D>>);

/// Check that vk.index_info and vk.index_comms are consistent with specified segment_size,
/// max_domain_h_size and max_domain_k size
pub fn is_vk_consistent_with<G: AffineCurve, D: Digest>(
    vk: &MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>,
    segment_size: usize,
    max_domain_h_size: usize,
    max_domain_k_size: usize,
) -> bool
{
    let num_constraints = vk.index_info.num_constraints.next_power_of_two();
    let num_variables = vk.index_info.num_variables.next_power_of_two();
    let domain_h_size = std::cmp::max(num_constraints, num_variables);
    let domain_k_size = vk.index_info.num_non_zero.next_power_of_two();
    let num_segs = (domain_k_size as f64/segment_size.next_power_of_two() as f64).ceil() as usize;

    domain_h_size <= max_domain_h_size.next_power_of_two() &&
        domain_k_size <= max_domain_k_size.next_power_of_two() &&
        vk.index_comms.iter().all(|comm| comm.comm.len() <= num_segs)
}

/// Check that proof is consistent with specified vk, min_segment_size, max_segment_size,
/// max_domain_h_size and max_domain_k_size
pub fn is_proof_consistent_with<G: AffineCurve, D: Digest>(
    proof: &MarlinProof<G, D>,
    vk: &MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>,
    min_segment_size: usize,
    max_segment_size: usize,
    max_domain_h_size: usize,
    max_domain_k_size: usize,
) -> bool
{
    let segment_size: usize = 1 << proof.pc_proof.proof.l_vec.len();
    let zk_bound = if proof.pc_proof.proof.hiding_comm.is_some() { 1 } else { 0 };
    let w_z_a_b_segs = ((max_domain_h_size + 2 * zk_bound) as f64/segment_size as f64).ceil() as usize;
    let t_segs = ((max_domain_h_size as f64/segment_size as f64)).ceil() as usize;
    let z_1_segs = ((max_domain_h_size + 3 * zk_bound) as f64/segment_size as f64).ceil() as usize;
    let h_1_segs = ((2 * max_domain_h_size + 4 * zk_bound - 2) as f64/segment_size as f64).ceil() as usize;
    let z_2_segs = (max_domain_k_size as f64/segment_size as f64).ceil() as usize;
    let h_2_segs =  ((3 * max_domain_k_size - 3) as f64/segment_size as f64).ceil() as usize;

    is_vk_consistent_with(vk, segment_size, max_domain_h_size, max_domain_k_size) &&
        segment_size >= min_segment_size &&
        segment_size <= max_segment_size &&
        proof.commitments[0][0].comm.len() <= w_z_a_b_segs && // w poly
        proof.commitments[0][1].comm.len() <= w_z_a_b_segs && // z_a poly
        proof.commitments[0][2].comm.len() <= w_z_a_b_segs && // z_b poly
        proof.commitments[1][0].comm.len() <= t_segs && // t poly
        proof.commitments[1][1].comm.len() <= z_1_segs && // z_1 poly
        proof.commitments[1][2].comm.len() <= h_1_segs && // h_1 poly
        proof.commitments[2][0].comm.len() <= z_2_segs && // z_2 poly
        proof.commitments[2][1].comm.len() <= h_2_segs // h_2 poly
}

/// Given segment_size, density, zk, proof_type, return the maximum value of |H| s.t. proof size will be <= max_proof_size
/// and the corresponding proof size and vk size
pub fn compute_max_domain_h_size(
    segment_size: usize,
    density: usize,
    zk: bool,
    max_proof_size: usize,
    proof_type: &str
) -> (usize, usize, usize)
{
    let zk_bound: usize = if zk { 1 } else { 0 };
    let segment_size = segment_size.next_power_of_two();
    let mut max_supported_h_size = 0;
    let mut max_supported_proof_size = 0;
    let mut max_supported_vk_size = 0;

    loop {
        let h = 1 << max_supported_h_size;
        let k = (h * density).next_power_of_two(); // |K|/|H| = density
        let w_z_a_b_segs = ((h + 2 * zk_bound) as f64/segment_size as f64).ceil() as usize;
        let t_segs = ((h as f64/segment_size as f64)).ceil() as usize;
        let z_1_segs = ((h + 3 * zk_bound) as f64/segment_size as f64).ceil() as usize;
        let h_1_segs = ((2 * h + 4 * zk_bound - 2) as f64/segment_size as f64).ceil() as usize;
        let z_2_segs = (k as f64/segment_size as f64).ceil() as usize;
        let h_2_segs =  ((3 * k - 3) as f64/segment_size as f64).ceil() as usize;

        if w_z_a_b_segs > 255 || t_segs > 255 || z_1_segs > 255 || h_1_segs > 255 || z_2_segs > 255 || h_2_segs > 255 {
            max_supported_h_size += 1;
            continue;
        }

        let num_segments = 3 * w_z_a_b_segs + t_segs + z_1_segs + h_1_segs + h_2_segs + z_2_segs;

        let num_evaluations = 22; // indexer polys (12) + prover polys (8) + 2 (z_1 and z_2 are queried at 2 different points) 

        let pc_proof_size = 1 // l_vec_len
            + 2 * algebra::log2_floor(segment_size) * 33 // l_vec and r_vec elems
            + 33 // G_final
            + 32 // c_final
            + 1 // Hiding comm is Some or None
            + if zk { 33 } else { 0 } // If zk we will have the hiding comm
            + 1 // Rand is Some or None
            + if zk { 32 } else { 0 }; // If zk we will have the rand

        let batch_poly_segs = ((3 * k - 4) as f64/segment_size as f64).ceil() as usize;
        let pc_batch_proof_size = (num_evaluations - 2) * 32 // 32 bytes to serialize 1 field element
            + 1 // 1 byte to encode length of evaluations vec
            + 33 * batch_poly_segs // num segs of the highest degree polynomial as the batch poly will have this degree too
            + 1 // 1 byte to encode length of segments vec
            + pc_proof_size as usize;

        let proof_size = num_segments * 33 // 33 bytes used for point compressed representation
            + 8 // 1 byte for each poly to encode shifted comm being Some or None
            + 8 // 1 byte for each poly to encode length of segments vector
            + num_evaluations * 32
            + pc_batch_proof_size
            + if proof_type == "darlin" {
                2 * // 2 deferred accumulators
                (
                    33 // G_final
                    + 1 // xi_s len
                    + algebra::log2_floor(segment_size) * 16 // xi_s (only 128 bits long)
                ) 
            } else { 0 } as usize
        ;

        if proof_size > max_proof_size {
            break (max_supported_h_size - 1, max_supported_proof_size, max_supported_vk_size)
        }

        max_supported_proof_size = proof_size;

        let indexer_polys_num_segs = (k as f64/segment_size as f64).ceil() as usize;
        max_supported_vk_size = 24 // index_info
            + 1 // indexer comms vec len
            + indexer_polys_num_segs * 33 * 12 // segment commitments for each indexer poly
            + 12 // comms vec len for each indexer poly
            + 12 // shifted comm some or none for each indexer poly
        ;

        max_supported_h_size += 1;
    }
}

impl<G: AffineCurve, D: Digest> Deref for MarlinProof<G, D> {
    type Target = Proof<G::ScalarField, InnerProductArgPC<G, D>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<G: AffineCurve, D: Digest> SemanticallyValid for MarlinProof<G, D> {
    fn is_valid(&self) -> bool {
        // Check commitments number and validity
        let num_rounds = 3;
        let comms_per_round = vec![3, 3, 2];

        // Check commitments are grouped into correct num_rounds
        if self.commitments.len() != num_rounds { return false };

        // Check that each round has the expected number of commitments
        for i in 0..comms_per_round.len() {
            if self.commitments[i].len() != comms_per_round[i] { return false };
        }

        // Check evaluations num
        let num_polys = AHPForR1CS::<G::ScalarField>::PROVER_POLYNOMIALS.len() +
            AHPForR1CS::<G::ScalarField>::INDEXER_POLYNOMIALS.len();
        let evaluations_num = num_polys + 2;

        self.commitments.is_valid() &&  // Check that each commitment is valid
            self.evaluations.len() == evaluations_num && // Check correct number of evaluations
            self.evaluations.is_valid() && // Check validity of each evaluation
            self.prover_messages.len() == num_rounds &&// Check correct number of prover messages
            self.prover_messages.is_valid() && // Check prover messages are valid
            // Check opening proof
            self.pc_proof.is_valid() &&
            self.pc_proof.batch_values.len() == num_polys
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct SimpleMarlinPCD<'a, G: AffineCurve, D: Digest> {
    pub proof:                     MarlinProof<G, D>,
    pub usr_ins:                   Vec<G::ScalarField>,
    _lifetime:                     PhantomData<&'a ()>,
}

/// As every PCD, the `SimpleMarlinPCD` comes as a proof plus "statement".
impl<'a, G, D> SimpleMarlinPCD<'a, G, D>
    where
        G: AffineCurve,
        D: Digest + 'a,
{
    pub fn new(
        // A normal (coboundary) Marlin proof
        proof:   MarlinProof<G, D>,
        // The "statement" of the proof. Typically the full public inputs
        usr_ins: Vec<G::ScalarField>
    ) -> Self
    {
        Self { proof, usr_ins, _lifetime: PhantomData }
    }
}

/// To verify the PCD of a simple Marlin we only need the `MarlinVerifierKey` (or, the 
/// IOP verifier key) of the circuit, and the two dlog committer keys for G1 and G2.
pub struct SimpleMarlinPCDVerifierKey<'a, G: AffineCurve, D: Digest>(
    pub &'a MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>,
    pub &'a DLogVerifierKey<G>
);

impl<'a, G: AffineCurve, D: Digest> AsRef<DLogVerifierKey<G>> for SimpleMarlinPCDVerifierKey<'a, G, D> {
    fn as_ref(&self) -> &DLogVerifierKey<G> {
        &self.1
    }
}

impl<'a, G, D> PCD for SimpleMarlinPCD<'a, G, D>
    where
        G: AffineCurve,
        D: Digest + 'a,
{
    type PCDAccumulator = DLogItemAccumulator<G, D>;
    type PCDVerifierKey = SimpleMarlinPCDVerifierKey<'a, G, D>;

    fn succinct_verify(
        &self,
        vk: &Self::PCDVerifierKey,
    ) -> Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError>
    {
        let succinct_time = start_timer!(|| "Marlin succinct verifier");

        // Verify the IOP/AHP 
        let (query_set, evaluations, labeled_comms, mut fs_rng) = Marlin::<G::ScalarField, InnerProductArgPC<G, D>, D>::verify_ahp(
            &vk.1,
            &vk.0,
            self.usr_ins.as_slice(),
            &self.proof,
        ).map_err(|e| {
            end_timer!(succinct_time);
            PCDError::FailedSuccinctVerification(format!("{:?}", e))
        })?;

        // Absorb evaluations and sample new challenge
        fs_rng.absorb(&self.proof.evaluations);

        // Succinct verify DLOG proof
        let (xi_s, g_final) = InnerProductArgPC::<G, D>::succinct_batch_check_individual_opening_challenges(
            &vk.1,
            &labeled_comms,
            &query_set,
            &evaluations,
            &self.proof.pc_proof,
            &mut fs_rng,
        ).map_err(|e| {
            end_timer!(succinct_time);
            PCDError::FailedSuccinctVerification(e.to_string())
        })?;

        // Successfull verification: return current accumulator
        let acc = DLogItem::<G> {
            g_final: Commitment::<G> {  comm: vec![g_final], shifted_comm: None  },
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(acc)
    }
}