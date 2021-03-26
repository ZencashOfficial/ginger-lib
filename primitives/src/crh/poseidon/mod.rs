extern crate rand;
extern crate rayon;

use algebra::PrimeField;

use std::marker::PhantomData;

use crate::crh::{
    FieldBasedHash,
    FieldBasedHashParameters,
};

pub mod batched_crh;

pub mod parameters;
pub use self::parameters::*;

pub mod sbox;
pub use self::sbox::*;

pub trait PoseidonParameters: 'static + FieldBasedHashParameters + Clone {
    const T: usize;  // Number of S-Boxes
    const R_F:i32;   // Number of half full rounds (R_f in the Poseidon paper)
    const R_P:i32;   // Number of partial rounds
    const ZERO:Self::Fr;   // The zero element in the field
    const C2:Self::Fr;     // The constant to add in the position corresponding to the capacity
    const AFTER_ZERO_PERM: &'static[Self::Fr]; // State vector after a zero permutation
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix
}

#[derive(Debug)]
pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>, SB: PoseidonSBox<P>>{
    state: Vec<F>,
    pending: Vec<F>,
    _parameters: PhantomData<P>,
    _sbox: PhantomData<SB>,
}

impl<F, P, SB> PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    #[inline]
    fn apply_permutation(&mut self) {
        for (input, state) in self.pending.iter().zip(self.state.iter_mut()) {
            *state += input;
        }
        // we do not use domain separation for now
        // domain separator for arity 2 Merkle tree hashing with all 2 leafs present.
        // self.state[P::R] += &P::C2;  
        Self::poseidon_perm(&mut self.state);
    }

    #[inline]
    fn _finalize(&self) -> F {
        let mut state = self.state.clone();
        for (input, s) in self.pending.iter().zip(state.iter_mut()) {
            *s += input;
        }
        // we do not use domain separation for now
        // state[P::R] += &P::C2;
        Self::poseidon_perm(&mut state);
        state[0]
    }

    pub(crate) fn poseidon_perm (state: &mut Vec<F>) {

        // index that goes over the round constants
        let mut round_cst_idx = 0;
        //println!("Full rounds:");
        // First full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                //println!("{:?}", state);
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            SB::apply_full(state, false)
        }

        // Partial rounds
        for _i in 0..P::R_P {

            //println!("Partial rounds:");
            // Add the round constants to the state vector
            for d in state.iter_mut() {
                //println!("{:?}", state);
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply S-BOX only to the first element of the state vector
            SB::apply_partial(state);
        }

        // Second full rounds
        for _i in 0..P::R_F {

            //println!("Full rounds:");
            // Add the round constants
            for d in state.iter_mut() {
                //println!("{:?}", state);
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            SB::apply_full(state, false);
        }

    }
}

impl<F, P, SB> FieldBasedHash for PoseidonHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonSBox<P>,
{
    type Data = F;
    type Parameters = P;

    fn init(personalization: Option<&[Self::Data]>) -> Self {
        let mut state = Vec::with_capacity(P::T);
        for i in 0..P::T {
            state.push(P::AFTER_ZERO_PERM[i]);
        }
        let mut instance = Self {
            state,
            pending: Vec::with_capacity(P::R),
            _parameters: PhantomData,
            _sbox: PhantomData,
        };

        // If personalization Vec is not multiple of the rate, we pad it with zero field elements.
        // This will allow eventually to precompute the constants of the initial state. This
        // is exactly as doing H(personalization, padding, ...). NOTE: this way of personalizing
        // the hash is not mentioned in https://eprint.iacr.org/2019/458.pdf
        if personalization.is_some(){
            let personalization = personalization.unwrap();

            for &p in personalization.into_iter(){
                instance.update(p);
            }

            let padding = if personalization.len() % P::R != 0 {
                P::R - ( personalization.len() % P::R )
            } else {
                0
            };

            for _ in 0..padding {
                instance.update(F::zero());
            }
            assert_eq!(instance.pending.len(), 0);
        }
        instance
    }

    // Note: `Field` implements the `Copy` trait, therefore invoking this function won't
    // cause a moving of ownership for `input`, but just a copy. Another copy is
    // performed below in `self.pending.push(input);`
    // We can reduce this to one copy by passing a reference to `input`, but from an
    // interface point of view this is not logically correct: someone calling this
    // functions will likely not use the `input` anymore in most of the cases
    // (in the other cases he can just clone it).
    fn update(&mut self, input: Self::Data) -> &mut Self {
        self.pending.push(input);
        if self.pending.len() == P::R {
            self.apply_permutation();
            self.pending.clear();
        }
        self
    }

    fn finalize(&self) -> Self::Data {
        if !self.pending.is_empty() {
            self._finalize()
        } else {
            self.state[0]
        }
    }

    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self {
        let new_instance = Self::init(personalization);
        *self = new_instance;
        self
    }
}

#[cfg(test)]
mod test {
    use algebra::{
        fields::{
            mnt4753::Fr as MNT4753Fr,
            mnt6753::Fr as MNT6753Fr,
        },
    };
    use std::str::FromStr;
    use algebra::biginteger::BigInteger768;
    use crate::crh::{
        poseidon::parameters::{
            mnt4753::MNT4PoseidonHash,
            mnt6753::MNT6PoseidonHash,
        },
        FieldBasedHash,
    };

    #[test]
    fn test_poseidon_hash_mnt4() {

        // Regression test
        let expected_output = MNT4753Fr::new(BigInteger768([120759599714708995, 15132412086599307425, 1270378153255747692, 3280164418217209635, 5680179791594071572, 2475152338055275001, 9455820118751334058, 6363436228419696186, 3538976751580678769, 14987158621073838958, 10703097083485496843, 48481977539350]));
        let mut poseidon_digest = MNT4PoseidonHash::init(None);
        let input = [MNT4753Fr::from_str("1").unwrap(), MNT4753Fr::from_str("2").unwrap()];

        let output = poseidon_digest
            .update(input[0])
            .update(input[1])
            .finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT4753.");

        // Test finalize() holding the state and allowing updates in between different calls to it
        poseidon_digest
            .reset(None)
            .update(input[0].clone());
        poseidon_digest.finalize();
        poseidon_digest.update(input[1].clone());
        assert_eq!(output, poseidon_digest.finalize());

        //Test finalize() being idempotent
        assert_eq!(output, poseidon_digest.finalize());
    }

    #[test]
    fn test_poseidon_hash_mnt4_single_element() {
        let expected_output = MNT4753Fr::new(BigInteger768([10133114337753187244, 13011129467758174047, 14520750556687040981, 911508844858788085, 1859877757310385382, 9602832310351473622, 8300303689130833769, 981323167857397563, 5760566649679562093, 8644351468476031499, 10679665778836668809, 404482168782668]));
        let mut poseidon_digest = MNT4PoseidonHash::init(None);
        poseidon_digest.update(MNT4753Fr::from_str("1").unwrap());
        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    fn test_poseidon_hash_mnt4_three_element() {
        let expected_output = MNT4753Fr::new(BigInteger768([5991160601160569512, 9804741598782512164, 8257389273544061943, 15170134696519047397, 9908596892162673198, 7815454566429677811, 9000639780203615183, 8443915450757188195, 1987926952117715938, 17724141978374492147, 13890449093436164383, 191068391234529]));
        let mut poseidon_digest = MNT4PoseidonHash::init(None);

        for i in 1..=3{
            poseidon_digest.update(MNT4753Fr::from(i as u64));
        }

        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT4753.");
    }

    #[test]
    fn test_poseidon_hash_mnt6() {
        let expected_output = MNT6753Fr::new(BigInteger768([8195238283171732026, 13694263410588344527, 1885103367289967816, 17142467091011072910, 13844754763865913168, 14332001103319040991, 8911700442280604823, 6452872831806760781, 17467681867740706391, 5384727593134901588, 2343350281633109128, 244405261698305]));
        let mut poseidon_digest = MNT6PoseidonHash::init(None);
        let input = [MNT6753Fr::from_str("1").unwrap(), MNT6753Fr::from_str("2").unwrap()];

        let output = poseidon_digest
            .update(input[0])
            .update(input[1])
            .finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT6753.");

        // Test finalize() holding the state and allowing updates in between different calls to it
        poseidon_digest
            .reset(None)
            .update(input[0].clone());
        poseidon_digest.finalize();
        poseidon_digest.update(input[1].clone());
        assert_eq!(output, poseidon_digest.finalize());

        //Test finalize() being idempotent
        assert_eq!(output, poseidon_digest.finalize());
    }

    #[test]
    fn test_poseidon_hash_mnt6_single_element() {
        let expected_output = MNT6753Fr::new(BigInteger768([9820480440897423048, 13953114361017832007, 6124683910518350026, 12198883805142820977, 16542063359667049427, 16554395404701520536, 6092728884107650560, 1511127385771028618, 14755502041894115317, 9806346309586473535, 5880260960930089738, 191119811429922]));
        let mut poseidon_digest = MNT6PoseidonHash::init(None);
        let input = MNT6753Fr::from_str("1").unwrap();
        poseidon_digest.update(input);
        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT6753.");
    }

    #[test]
    fn test_poseidon_hash_mnt6_three_element() {
        let expected_output = MNT6753Fr::new(BigInteger768([13800884891843937189, 3814452749758584714, 14612220153016028606, 15886322817426727111, 12444362646204085653, 5214641378156871899, 4248022398370599899, 5982332416470364372, 3842784910369906888, 11445718704595887413, 5723531295320926061, 101830932453997]));
        let mut poseidon_digest = MNT6PoseidonHash::init(None);

        for i in 1..=3{
            let input = MNT6753Fr::from(i as u64);
            poseidon_digest.update(input);
        }

        let output = poseidon_digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match for MNT6753.");
    }

    use algebra::{
        fields::bn_382::{
            Fr as BN382Fr, Fq as BN382Fq,
        },
        biginteger::BigInteger384,
    };
    use crate::crh::parameters::bn382::*;

    #[test]
    fn test_poseidon_hash_bn382_fr() {
        let expected_output = BN382Fr::new(BigInteger384([5374955110091081208, 9708994766202121080, 14988884941712225891, 5210165913215347951, 13114182334648522197, 392522167697949297]));

        let mut digest = BN382FrPoseidonHash::init(None);
        digest.update(BN382Fr::from_str("0").unwrap());
        digest.update(BN382Fr::from_str("0").unwrap());
        let output = digest.finalize();

        assert_eq!(output, expected_output, "Outputs do not match for BN382Fr");
    }

    #[test]
    fn test_poseidon_hash_bn382_fq() {
        let expected_output = BN382Fq::new(BigInteger384([10704305393280846886, 13510271104066299406, 8759721062701909552, 14597420682011858322, 7770486455870140465, 1389855295932765543]));

        let mut digest = BN382FqPoseidonHash::init(None);
        digest.update(BN382Fq::from_str("1").unwrap());
        digest.update(BN382Fq::from_str("2").unwrap());
        let output = digest.finalize();

        assert_eq!(output, expected_output, "Outputs do not match for BN382Fq");
    }

    use algebra::{
        fields::tweedle::{
            fr::Fr as tweedleFr, fq::Fq as tweedleFq,
        },
        biginteger::BigInteger256 as BigInteger,
        Field, 
    };
    use crate::crh::parameters::*;

    #[test]
    fn test_poseidon_permutation_tweedle_fr() {
        // test vectors are computed using the sage script of  
        let test_ins = vec![
            vec![tweedleFr::zero(); 3],
            vec![ 
                tweedleFr::new(BigInteger([0x4b6843cb99b405a7,0x9144f963e6be5b2f,0xd6cd171f511f42ae,0x3236f7bf7169c287,])),
                tweedleFr::new(BigInteger([0x5b3a70c1157f1212,0x12dd719b00c4df79,0xbd19cef4c96d3866,0x30ee9801b0360968,])),
                tweedleFr::new(BigInteger([0x8bc11448cb71832,0x1d9f276875cc842a,0x42fd3aabdab2b6ea,0xc0424292f53f9bd,])),
            ],
            vec![
                tweedleFr::new(BigInteger([0x6041a866e4558c17,0x3a786f845d979e7f,0x5d6b1332c8f927b0,0x16fbd749aad1f957,])),
                tweedleFr::new(BigInteger([0x8f35db97f2816fa3,0x80d16d81b03fdada,0x1379a31c4728cc76,0x38be09a3a66ff4c2,])),
                tweedleFr::new(BigInteger([0xe566ac71952ead83,0x37b7521f0ecab296,0x7ec918e8f984ab0b,0x12c52d236863265b,])),
            ],
            vec![
                tweedleFr::new(BigInteger([0x4c81a064348dd3e7,0xf0a3dca1cdacff2b,0xd364cf53239ffde,0x26d11588ffb19887,])),
                tweedleFr::new(BigInteger([0xfdfe267f262ca737,0x91eaa0df803a9dc1,0x254ed4b68d1fd2cd,0x20799470abb1adf7,])),
                tweedleFr::new(BigInteger([0xabf88cf4a0dfef16,0x659d666776f371e5,0xf52159ce419ae46e,0x231a73b9cd584d51,])),
            ],
            vec![
                tweedleFr::new(BigInteger([0x375748ce5dce48f3,0xd056583705f71847,0x9ee4a117c59be931,0x11480580568757d7,])),
                tweedleFr::new(BigInteger([0x74e3326574f2c3cf,0xaaee18cf2c9d6835,0x433459495ad3cd5,0x16187918f8142b98,])),
                tweedleFr::new(BigInteger([0xa30844cb51aa9776,0xef0e1428779ee89a,0x36745453414361b6,0xb5034593923cb50,])),
            ],
        ];
        
        let test_outs = vec![
            vec![
                tweedleFr::new(BigInteger([0x24474df3c64b6467,0xc0c46698dadb69fb,0x89479fdaa2d4a3d3,0x1038e2503c010653])),
                tweedleFr::new(BigInteger([0xdbd0f46c618cbe3f,0xd4c1353c7911ca30,0xc7d9c14b932bea56,0x3e8a462618bedbfc])),
                tweedleFr::new(BigInteger([0x24022ec1dbfc566b,0xed45883272d605ba,0x147a173707407b8b,0xebc6f1edf47ed55])),
            ],
            vec![
                tweedleFr::new(BigInteger([0x20fa7060162fe3fc,0x1eac682e3c5faa3c,0x9c5bceb01d30f9ba,0x2b3eabe68941720,])),
                tweedleFr::new(BigInteger([0x7c1779b6127cf06c,0x9f32d6ad1f28799e,0x4731c51b69013a3a,0x3d1ed3edc77b97ae,])),
                tweedleFr::new(BigInteger([0xe952f46225b93cf7,0x4ec21947f91d9ca7,0xcd24b17ef0c76ae2,0x23e7221b28b210f3,])),
            ],
            vec![
                tweedleFr::new(BigInteger([0xcf79d1777224c062,0x52d39a745e1c3c6a,0xd941fc13d0c346be,0xb699418e5ef9834,])),
                tweedleFr::new(BigInteger([0xe2396e61bfda0277,0xeaed6a4e968ea195,0x3f08907457a3236f,0x3c63cead78221117,])),
                tweedleFr::new(BigInteger([0xe90a5898600235d9,0xaa897713652ae7b0,0xc5a256f19615f569,0x3bc13c2bfd63df7e,])),
            ],
            vec![
                tweedleFr::new(BigInteger([0x66d425888381162f,0xc63132cb5028045b,0x53871f80990aedbc,0x363207508e5b631b,])),
                tweedleFr::new(BigInteger([0xfe5d1e98b1530279,0x27c024f1831bc8e5,0x3d1cd5426dab4aee,0x82529e25a79f601,])),
                tweedleFr::new(BigInteger([0xb0326dacb304b76e,0x28014ea114458515,0xee0838406e061d3b,0x16413d5dc81bd4e1,]))
            ],
            vec![
                tweedleFr::new(BigInteger([0xd51f314c6d11765c,0xae73a2e659ab5842,0xe87e6d45db318fdf,0x2fc5b82cbe80b30e,])),
                tweedleFr::new(BigInteger([0xb3f2983fa2e73b8c,0xf674add1d5cfba09,0x4e378c14f155766e,0x3ef85cc9fe99b0dd,])),
                tweedleFr::new(BigInteger([0xf61caec9a1007cf8,0x81f1b58a6c336283,0x24b14f97c35bc880,0x398b61120ab4fcc,])),
            ],
        ];

        for i in 0..test_ins.len(){
            let mut state = test_ins[i].clone();
            FrPoseidonHash::poseidon_perm(&mut state);
            assert_eq!(state, test_outs[i],"Output does not match tweedleFr, test vector {:?}", i);
        }
        
    }

    #[test]
    fn test_poseidon_permutation_tweedle_fq() {
        let test_ins = vec![
            vec![tweedleFq::zero(); 3],
            vec![
                tweedleFq::new(BigInteger([0xf442d87fbf3ecd8c,0x75c6ac589dfecc42,0x83effdb327a6dfef,0xbad6a789fe330dd,])),
                tweedleFq::new(BigInteger([0xc3abc1cdf4935f44,0xda31dc8bb5c04e2e,0x3417d59cf717c1eb,0x361003499c842ef6,])),
                tweedleFq::new(BigInteger([0x41792959d4240eec,0xeca34dfbd3be7ceb,0x96ab342287a71ec5,0x2cd32749a0efc9f3,])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x8a42a591ccdfea4f,0x7fac93729f624247,0x6e3610f78e23ee57,0x20e05d3a71be99b2,])),
                tweedleFq::new(BigInteger([0x193bc9eb2db311b3,0xa11070475db1a67f,0xc44fc3968a8c8ddc,0x24d8ea97bb23fd7f,])),
                tweedleFq::new(BigInteger([0x10b4ddca66c0b564,0xa8728f932bb40ee6,0xba11d94c27eabac6,0xd3a13af304669ac,])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x6c255f28b95bd1e6,0x69504c6efad8bdb7,0x8f1ad7644227b830,0x19939485311e6a40,])),
                tweedleFq::new(BigInteger([0x55b0edfcd04e8978,0x53ebab97a87e54fb,0x677f2fc34fa46db3,0x26d2700a1a4bb573,])),
                tweedleFq::new(BigInteger([0xdfc52a04cf954d1e,0x72113c3d0ce20967,0xb72c15b74381d6a0,0xb0819d8686c00f9,])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x723df3e990dca9d1,0xc92ae0972587d963,0xd90188b7354a2258,0xc0aa1bfc0aa6c97,])),
                tweedleFq::new(BigInteger([0x446d4bab99ee5bff,0xdba6ba04335e124e,0x347790bc7c8eb68a,0x598e0584f1b2f97,])),
                tweedleFq::new(BigInteger([0x1559a27f16b7abfc,0x164beb53557e0c40,0x3ff0f79d7fdb5c20,0xb362cb07a9be10f,])),
            ],
        ];
        
        let test_outs = vec![
            vec![
                tweedleFq::new(BigInteger([0xb46ad2d321d1cefd, 0xbb1203b20c8e17ab, 0x7413c99cee248276, 0x10f3026539fe9e4])),
                tweedleFq::new(BigInteger([0x3489f06408fb2c18, 0x668c8ed351f01494, 0x944306c7ff00efae, 0x1e58dc79b4fcae6c])),
                tweedleFq::new(BigInteger([0xad5e1e851822cc05, 0x84d9ce927f225195, 0xdd4ab88f002f38dc, 0x19f0f4588a618830])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x123848facefaabd9,0xea61b7b6e2aba395,0xbac0891688b620bd,0x32cb9ee8f31699a7,])),
                tweedleFq::new(BigInteger([0xcf43eb1a4b301683,0x3fd3d32fd2e7e72a,0x8ab0423c10f4691d,0xf1a94dcc8177c3f,])),
                tweedleFq::new(BigInteger([0xd1a73d73eef198f,0x488a01394d86adf0,0xe08f1081d924d423,0xb9cff6f01f3ee6d,])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x27551948e97f0eb,0x84f7c9012be83f48,0x15e7b16a83d00ce8,0x36fccebf901a1145,])),
                tweedleFq::new(BigInteger([0xc2e6b8fecab6bb12,0x1414ee62b4279c08,0x2fc8cf5ff3f57c39,0x3631f456f95d8003,])),
                tweedleFq::new(BigInteger([0xee418c2224438111,0xa5fcfec23d9921fd,0x803975a67437d3ee,0x32c246f59e9ab177,])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x47f1223f3c2a4b1b,0xc8da66fcb1167479,0x2084bba52ace1e43,0xf17072bd192088e,])),
                tweedleFq::new(BigInteger([0xcf104b943e3dbc6c,0x95082a7b3681ac7e,0x54bbe8f15ecca259,0x1d5381c193218ac0,])),
                tweedleFq::new(BigInteger([0xa4614e0357afcbbb,0xb6bb7a99d65c17f1,0x39a3ed06dce8d892,0x3fde42201e273c66,])),
            ],
            vec![
                tweedleFq::new(BigInteger([0x3ba1131ea507bc64,0xdb5361627ca99dde,0x1f2d47f63d54ae70,0x3a6e4e5385a085c2,])),
                tweedleFq::new(BigInteger([0xb4b6cc077add652,0xc47e64b8b1316ff,0x975ab087fcecdd58,0x3641f4c76b53c84e,])),
                tweedleFq::new(BigInteger([0x5f6190e71f11bf05,0x8be71bf7461f1ac4,0x513e0394d4bde9d1,0x3ecc4ef09355de45,])),
            ],
        ];

        for i in 0..test_ins.len(){
            let mut state = test_ins[i].clone();
            FqPoseidonHash::poseidon_perm(&mut state);
            assert_eq!(state, test_outs[i],"Output does not match tweedleFq, test vector {:?}", i);
        }

    }

}
