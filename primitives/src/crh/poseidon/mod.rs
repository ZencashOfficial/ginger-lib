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
        println!("Full rounds:");
        // First full rounds
        for _i in 0..P::R_F {
            println!("{:?}", state.clone());
            // Add the round constants to the state vector
            for d in state.iter_mut() {
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
    fn test_poseidon_permutation_mnt4(){
        // test vectors are computed via the script in ./parameters/evidence
        let test_ins = vec![
            vec![MNT4753Fr::zero(); 3],
            vec![
                MNT4753Fr::new(BigInteger768([0xf770348fbe4e29b6,0xfefd6b30dfb52494,0xec61827e5cf9425,0xc6288db72079112c,0xd70e11f75c351bac,0x2e4657caf8648c8e,0x7f9f3a94358aa2f7,0xee7f886bb42e8eab,0xe5ae5d4ec1b0796f,0xd056464cb38777c6,0xf3d7cd676c74ae38,0x120d49a741c34,])),
                MNT4753Fr::new(BigInteger768([0x96de60f9741f78b7,0xa98cc9495bb4615e,0xc4b3aeadfd321c2c,0x40e4b75eb8fe1116,0x1396ee290297e819,0x9762744e4cfded19,0xbedcef99b43ee15a,0x8b84865c31d378a0,0xf5468754aa4a4c4e,0xfd715c8245c2e124,0x31cb5bb04a339986,0xdaf306180aed,])),
                MNT4753Fr::new(BigInteger768([0x7e874134d509e406,0x729d013268020212,0x8b362dd530097799,0xae5054da3ad04250,0xe2e7413bd0fcbe5f,0xad08673f2f925bee,0xfb93f0ee8900d97e,0x2c1d037343b00151,0xd3dac3f2b1139f55,0x154e788ae1aca4cc,0x663269814fb52d57,0x676d9c4d8329,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xa26b0bc72724d615,0x729202dca25403d4,0x1b2ff6dc78c46b5e,0xed529329c88557ec,0xa7264c3cd1f1ca2d,0xa9f0e2b1e57c800f,0x2322b96082d360ec,0x138d00037c082f1c,0x6c25792c21edce0a,0x75723fc00d8d1bc3,0xf60868fea31de240,0x14e224d41e354,])),
                MNT4753Fr::new(BigInteger768([0x21c229d68cde6f3f,0xf96b852ba3677e55,0x815b51e9b5e329c2,0xedec4ec2b77a9d36,0x44e0217411a0dea9,0x724a35de8cbd3141,0x8008cb5f0106c484,0x921855777c1c9cd3,0xd87d5b5babb7c9ab,0x603fc082a06ed8c4,0xe589b5a1adea946e,0x129d1f84a0c66,])),
                MNT4753Fr::new(BigInteger768([0x80794339ccdf973f,0x8f537759fc1b1aca,0x7997a170b362d649,0x7b1cddf6db6ca199,0x6b25316a81753330,0xa143d6d50bd07ebf,0x4d65e4fd6f8587d6,0x572c858cf606bd90,0x245465ba33e044b1,0x86f9aaa423b9390,0x8ee2bbed6bda13a6,0x7fa83fcd7a59,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0x275345cd3949fba9,0xaa492ccf37b80d9,0xdd9c6b17371c879a,0x846303d5b851d739,0x8d2b1b900c8c2227,0x780824b721514171,0xe08b4ffffb8a4f71,0xc69a0eb1b3f3ad,0x409578a5de88b1df,0xef2b552006465afb,0x2539560ecdf8147,0x134fe3e183dcd,])),
                MNT4753Fr::new(BigInteger768([0xf7f3c59f70e5b72a,0xec1ae7ed077f2d99,0xbbf075b432e1a2d8,0xf32012c620b8cd09,0x81e964a2687b8654,0x43082373cc23c4f6,0x494428fd5d2b9d5,0xed89d49a5f32ca1a,0x8d2c7f6937d4bc08,0x8aa8316d21567c0c,0x5e2c9cde56f4c802,0x6422f65bc889,])),
                MNT4753Fr::new(BigInteger768([0x44238a7e541cdf0,0xc09a1bda2e310a6d,0xef2001005bbaf873,0x1fd97ee19fea97eb,0xce43458dee7839cd,0x735d8cff80565348,0xca740dd90f883e06,0x8825f23c63c39a44,0xe80c50eb3548e408,0xddc815aae7e6a432,0x519048208b84f07f,0x50d352305dca,])),
            ],
        ];
        
        let test_outs = vec![
            vec![
                // MNT4753Fr::new(BigInteger768([0x52626a743486529e,0x9ad956d530d8401e,0xef5eb8658628e928,0xeb5f8d670ce9a52e,0xf9e3512d1adde5be,0xca64a145072b2892,0xfdd8382f678a4884,0xfc64e6b3be14b979,0xaeae6274012a2b15,0xf6c694b5b96d56ac,0x8dd50ddfc9ec4417,0x1b605a49b041a,])),
                //MNT4753Fr::new(BigInteger768([0xb23e0e9f9cbeeeff,0x92b368da903a94cf,0x7dce7bd7bbbd42ee,0xd2130eac5a761a7,0x49373fbcf5b77fc2,0x11ae17d13a7f6f0b,0x2866bfc968ce04ce,0x5bb1a38f29603f30,0x17bab9b817fa37ae,0xd80efe15dd474cac,0x69045036fd2183a1,0xc78e7c434288,])),
                //MNT4753Fr::new(BigInteger768([0x72bbe685d4019804,0x4dc4bedc884b1f8f,0x25a7be14ef623727,0xac801a19d210202c,0x305be6132e688893,0xfa875b479565d67c,0x9e083411ac146578,0xca9ed4ed4acec81,0xf5ce61769933636b,0xde7aac3fe677a913,0x20ab7e29781eb546,0x19f269e3ece48,]))
                MNT4753Fr::new(BigInteger768([0x4f54c026da6ed8f0,0x12700bf5ad94f6c9,0x23a3fa62e9c042c1,0x2394c785581c75e7,0x839626f16bd60d08,0xb29828eef68c9bd4,0xd1479004b0f71d2,0x9d1a0dffdd1e7b00,0x9f1df2af9215e68c,0xc562186972253d2e,0xf6b8c66a6f3999b0,0xa040e4e0ff92,])),
                MNT4753Fr::new(BigInteger768([0xb0258a782c08064,0x6a04841f8be4990a,0xda027778a67d713b,0xb88be63be3cac9b4,0xde929c2510a321e5,0xc0d9dd704886213e,0xfbe0efc728d44f11,0x77c8d6422b5eb368,0x2827d5af4fe0fbad,0xb90c8793bc2a9d21,0xf9ce1fdde5140214,0x15a64a6345311,])),
                MNT4753Fr::new(BigInteger768([0xde9731dd4ad29db3,0x86caaccf88b402a1,0xe5e77eee08fca8a2,0x1dd9e752e50aad07,0x2d0f73cfb9508a83,0xb2b6ab08f14d96eb,0x224833c17d87490d,0x4e7d2e13141aaa55,0x1796b61e1cc3563,0xdbeb6f5ed60179f,0xb0633f07c680eda2,0x601b999b7143,])),                   
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xe749d7517ebe099b,0xc6abeacc602cf0bf,0x958f4b91f3c3b22d,0x9a295b36c4a6ea9e,0xd3925085d5ae2179,0xf23a8b4284968652,0x8018232a8a8fd30b,0x34533842150d4c6a,0xf0531c8f2f4a3dd4,0xeaab2b7956c6e7cb,0x9fc2b52eb516b457,0x7e2c759defce,])),
                MNT4753Fr::new(BigInteger768([0xfc5dab1dedb49656,0x78deb85913893c98,0x6088942fdbff357e,0xb3c15f514de46072,0x5dc205c3ccd4df39,0x591d9320bec689a6,0x99a7765caae47a86,0x2fcfe60a560fa3ed,0x43e2f302b5852456,0x5b4087eaa01f39c6,0xcc7db3f671985b7d,0x1272366ae322b,])),
                MNT4753Fr::new(BigInteger768([0xc23a10d72a73058e,0x7125f89599d62e8e,0x944ffd3948d3b453,0xc1513ee7ef29c1d2,0xdf1ddf8a25a2233,0x193c0cac56b49055,0xcb23ffde25ea2bd6,0x6d4a4ad2f3e415af,0x7da1b50b3731057,0x30f2f41a6746bd09,0x2a3cfda1f9885424,0xe6f1af34a223,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0xbfcb18d74e65c563,0x722359395bfeb077,0xb8e0b7abddb9a694,0xc830a386c2854b6b,0x53d7c0704e145ce,0xbe91d2a17d6f8874,0x2b49e38e1b99292a,0xc7e2cb48c2be1151,0xa5e54b3a714aad54,0xf634e385fe3d9b90,0x66f9a11a59535867,0x1425351d064a2,])),
                MNT4753Fr::new(BigInteger768([0x4a28ff3c4fecbb8d,0x60a639f0a2a002d9,0x5149d27ed99128c1,0x6dacfe4ce235b503,0xf21ef2fe6f344e69,0xbac70a5d64a033de,0x54f1cb89e291c8e6,0x2548230a2b8eeb67,0x763440a89ffdc8de,0x3ac6435a7c2b7922,0xacb97881f998663d,0x8ae31b1e760f,])),
                MNT4753Fr::new(BigInteger768([0x9dfe82b5a7baefa5,0x14bff3144e3c4f00,0xcbb47c1db66e74c4,0x8c3d330245b24464,0x3be7110fcc0f2674,0xb4a9281c6d349356,0xa4894a010cef488c,0x2abe0a21b8a83ca7,0xf9e9d807e418b54,0x439e4046be879838,0x3204e13287f737d5,0x3098a5738444,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0x470bac44ae262597,0x37c75eb3f00758fb,0xae77bbd563b5fac6,0xa22469cb36563eb5,0x4db9a5ea229af500,0xf6848cf2a64ad4a5,0x3a4611a0ed9e6243,0xf63fb5b6489325dd,0x1a9c90dd1544863f,0xdab1cb220fdf73d4,0xb9ec40309591932b,0x141777a73c602,])),
                MNT4753Fr::new(BigInteger768([0xedab7a7bd3a0061b,0x32d0ba278e569bec,0x83a9e0f317060812,0x29acd35e4d33cdb6,0x3f13496b623a9cde,0xa565606e05e4a5d,0xba87579c189af741,0x45bcb5fbad648a4e,0x32e8658135401638,0xbc853abb54e732b5,0xc37855ec443e12d3,0x1ad1ff8f54ad6,])),
                MNT4753Fr::new(BigInteger768([0xaba94817dccf0311,0x601cdff2f1e54d9e,0x6a0d8ab8a097a5b6,0x51d8c83d12239512,0x92f9ef537fc921e8,0x688b9fe86605c3ae,0x250ebdd755ad043c,0x29d412ee38a1e765,0xb31f5447678264b4,0x7d053f0ea44d854b,0x5d83d881795db690,0x397b9db5b588,])),
            ],
            vec![
                MNT4753Fr::new(BigInteger768([0x911b5559a3eeb52d,0x482afb0b1b566e49,0x3983c4efc4fb37da,0x3288b81e77372d01,0xc69bd18751793c34,0x103f732ca150f840,0xbe72b866f7fd8512,0x19f4e9f908c9d1bf,0xb7976427cfc0fe4e,0xc9f43b7c2ad54601,0x3f2eb373787a291,0x9d3dd62a7475,])),
                MNT4753Fr::new(BigInteger768([0x799693496d2180d4,0x9c8364f338a500b7,0x37a57ca5674e1252,0x2c19b0502325bead,0x32b30a126f41f5ac,0x8bcd51ff52cedf29,0x9e04cb66d8d16160,0x59e8aaadbc99fab6,0xbd046f342e99d386,0x4488dd3ce29590aa,0xdcc2bb0149b02eaa,0x1543d162aa244,])),
                MNT4753Fr::new(BigInteger768([0xbb41e5acd82643f9,0x4042aec0d83f7624,0x2c14ed2f563bb21e,0x9cee7ec494eb57e9,0x41eec6c2b0056ac2,0xd1ea7cfa30f223ef,0xf148c377c2fba415,0xb3b56ee96972c9cb,0x82c3e44086911217,0x9ef750feb5842cc6,0x9f33c28feb810dc0,0x727b9f80e6df,])),
            ],
        ];

        for i in 0..test_ins.len(){
            let mut state = test_ins[i].clone();
            MNT4PoseidonHash::poseidon_perm(&mut state);
            assert_eq!(state, test_outs[i],"Output does not match mnt4-753Fr, test vector {:?}", i);
        }
        
    }

    /*
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
    */

    #[test]
    fn test_poseidon_permutation_mnt6(){
        // test vectors are computed via the script in ./parameters/evidence
        let test_ins = vec![
            vec![MNT6753Fr::zero(); 3],
            vec![
                MNT6753Fr::new(BigInteger768([0x2045f548c283a386,0x9f90d7b623ef9965,0x634e1e0bcd6ce5f1,0xed09fb1cd92e2f48,0xa4b92193ab3a4c,0xc38d823f5e556d81,0x93e8a09384f1d5f0,0xa463757a137a2127,0xc948555766dabe44,0x3246e78f29a70bfe,0x21ebc006f85e213,0x18c2c2170055e,])),
                MNT6753Fr::new(BigInteger768([0x5abb4a33f5026781,0xa3510b40fb1bd1e7,0xce8ae77f3e0e9a1d,0xd1375569096b196a,0x107156721a5241bd,0x82b75d6eb65ccdc,0x9f6a6933bbe7d8ad,0x9335a61a85fe8998,0x5179ec766656404c,0x8052414d46077e77,0xb77841abce4c69c,0x10e71d39ef7ee,])),
                MNT6753Fr::new(BigInteger768([0xf76a1c08fa236882,0x25e1b757eb33ed43,0x1f63d4997a13c8b1,0xe23eae7ea2605b4b,0xe8c20feb190f9dd,0xa63856368a5c24f9,0x114eaf0c94cc670b,0xe858d17f6da22272,0x9b5443cadda8156a,0xfe92bd2a3eefc8b3,0x2c8a4defc4a4ff9,0x19cc15d056674,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x1b940903c57e8e7f,0xbd38cde2e8e16008,0xe18d1abcfe05990a,0x8e86b1ca3a0ee1f5,0x33a31929417f05f9,0x170be227265f62bd,0x29e22c2b9864352a,0x901db3c41b27245e,0xc3bc6e6cfce69e3c,0x498f01eea65c0215,0xbf86a87e3005b3db,0x90f488bd8e09,])),
                MNT6753Fr::new(BigInteger768([0xb2d9ad48cbb812ba,0xc53cb754a7a02d89,0x89f52c6630ad8f86,0xe623c68f3610652f,0x198f83682c814e5d,0xfb78854e850e95fb,0x46e398cb56c27f78,0x81e60dab3991f035,0x3babbc1fe35f4f30,0x8056c683be44ffab,0x167af8aceb070f00,0x1a2572baaf46d,])),
                MNT6753Fr::new(BigInteger768([0x6242acf3bfbe2c6e,0x7afcb4878b2fcab1,0xccdee01e7839e6ff,0x8ebef555a3fcaeb9,0xa627b970cb4d56d2,0xb672bd365dab0d61,0x71f74eef13dab0fd,0x5a138a0bd718f4c3,0x7d08a2cf2ef0747c,0x8a0cdeefcdfded66,0xfe18f6573bbabadb,0x12c02029e0030,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0xf2ca60b5bb454d9f,0xb4ae3ba59e4a711,0x62154368b888061c,0x6214f711b35b4f9,0x5dd4d44dc9d4f0ad,0x4304e1c271f64602,0x80d4e3b0e1025ae3,0x5316732d6accc44d,0x24fc5d7d7bba464e,0x12d10c9485d208a1,0xca6df371c62a8872,0x86ce9f608bae,])),
                MNT6753Fr::new(BigInteger768([0xcdf0f7492613b504,0x455faa0e541fa1e6,0xb77242df6b8a68be,0x3b5435160d723cb6,0x77b8914a813586bf,0xc17dabd68e491d26,0xa85720ce2231af9d,0xd19e81cea5d64c41,0x56c90bfdb43ce182,0x9ff4ff3aba6a9a01,0x8875906afee26342,0x16a993a8df862,])),
                MNT6753Fr::new(BigInteger768([0xad98e2452d8be558,0xed19ce15ee0069d3,0xf889b49a8ad1016e,0x42760a3cbfb291b7,0x3d94e422b333dc5d,0xc27cbbac2884c097,0x851fd495c84543e9,0xf9b100c34675f211,0x11eae122f8ff1706,0xf3eecc4f60743020,0x38fc6ca1e5d1b4a7,0xffa8124e7034,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x376743561f755f41,0xf0a8830457e9879b,0xa134b300b8f2d67b,0x1806324852aa9feb,0xdb98705dbf283859,0x565bca638d85ee57,0x1c6b6fe1fe752b0,0xd9eb23176d5d6110,0x5c5e86b5774422e2,0xd6fdf4c92ea236a1,0xeb2a915f44b72fa3,0x195c5f80dbf29,])),
                MNT6753Fr::new(BigInteger768([0x4c674244dfb68ecc,0x24a104856852ac3f,0x83e10d6c10dd3f4f,0xe99fe1f0d8553d3c,0x2d371b923253d5c0,0x14594932de89a19e,0xfd4589d2f8e53f17,0xe2ba2c7b929a53b3,0x3891f35b974a36ec,0xf17f8749ca140c09,0x6be74c21301f7c9e,0x13de4e1311a04,])),
                MNT6753Fr::new(BigInteger768([0xc366ce203caca4b0,0xe1d195b5bf3af54e,0x24b93c34bd0043ee,0x91559c070b29c53a,0xe866e46830168ff8,0xaeeda2129518cab7,0x37f8bb28ae15d7f3,0x5811fb22acd02c55,0xce7d805057f58acc,0x3a80df0b2af5f4fd,0x4dc7c29c8f6bed72,0xe511723afdb9,])),                
            ]
        ];
        
        let test_outs = vec![
            vec![
                MNT6753Fr::new(BigInteger768([0xef99f18ca1164fb0,0x1bf161755d689806,0x83ee017c500c6964,0x8abab822f92200c0,0x4b64884b9cc7eef9,0x53d4a2f13e17017c,0x551b8da2668dad8a,0x9939a48a0191c96c,0x2e1d80ef403671a0,0xb037bb60fbeb0212,0x6a22eba60581eb12,0x6ec196c9026d,])),
                MNT6753Fr::new(BigInteger768([0x18c4207483ba0f2f,0x6c50abc8aca74de3,0x7c1acfd6686351c,0xf367937c1356e91f,0xcdbf0447592ec1,0xe13763baac982387,0x2e1f904290e7045f,0xb6ffbcccd73c1092,0xfae22550de44cf2c,0x14c26231e52c7eae,0x471836049049f3b7,0xdc46826797ae,])),
                MNT6753Fr::new(BigInteger768([0x2ee4a96e4cda5f6f,0x7442a7b7f51fdbfc,0x23d03839ab7d811,0x1f873a8c0ddfd7a4,0x872f14e24612551a,0xd43181c852d5f78b,0xb2ff35a74130d2cd,0xd64aaa80f389157,0xb954953b8d35d74,0x37aba7a7212e96c,0xcce2fff62e11a3d4,0xfb3f9157120d,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x626e4d0e6e3e1936,0x7c99da459f8385d0,0xbd84a2fb934889a6,0xff40b1979118e180,0x76cb8b37a32cce54,0x6c389f3f88157389,0xb9f0135ec3d92cc2,0xfd6a928e603a79be,0x5472af35b978d0a6,0x109995c9831f98c2,0x976c556bfe34da5a,0xf838693b701,])),
                MNT6753Fr::new(BigInteger768([0x58fb485fd781fcc6,0xd92a60427ce67147,0x2cca412943d39ade,0xc55d3362bac1743,0xcb8dcfa4ae0fcda1,0x25bde06b8f99facd,0x2d30b30add5faa3e,0xbe0ebdda1ba7458d,0x296f6010c1db1c7b,0x506364ec0031a00e,0x24c13847d3fe6ab7,0xea0c23423f1a,])),
                MNT6753Fr::new(BigInteger768([0xc36816e6dafa2f57,0x554255a6e34a34d4,0x29f17ff72b3c5695,0xae97815a3cc05077,0x64a0824e4b9b1aae,0x267cf597a9a556ef,0x8d8c67fc33757cbc,0xad2db4d1a3c73012,0xf3fcee4d169de439,0xfc4632cd5cb31baf,0xe1420a2c4e68de6,0x1bd34ad51cd02,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x160dacef01b59531,0x313dd55066586bd8,0xdcc16189ec00c953,0xcc44967095828982,0x1066ee6f582ba5ea,0x3d879be40c078337,0xb9cb0ef83e1b4a51,0xc9b91de1e758c41,0xe578ceb8440e2bb8,0x3d6f2d210d4278df,0x2bab83b243a3335a,0x1afd20a9dbdc7,])),
                MNT6753Fr::new(BigInteger768([0x3a7ee60628dc201d,0xae1dcd081da757a,0xde1625ce6e93bc19,0xfb1a64dd14c0ae77,0x1bb5eba30eb2f202,0xdf064e762ce2f903,0x9abc764fb4c55d03,0x6db04d43d811c05d,0x87d85ec650763745,0x1bdcd095b0e1ada2,0x8681985565baa005,0x154d78a914323,])),
                MNT6753Fr::new(BigInteger768([0x101437542e4c39d4,0xcbdcf8d57d75fdd2,0x40996ed826c3b401,0xe492943442e0833b,0xf088ed10c7619f8c,0xb8e27256e0a69172,0x7112494180a5924,0x58d0e045a50972e9,0x4285049c582ed300,0xba0daceb8ab6d3c0,0x5ebb479b97c4c24d,0x820fdfe15d33,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x645f79445d3423f1,0x699de15f996c470c,0x3740c3b7e7818751,0xac5c029dba988fd2,0x7342c873ecef9aee,0x4ff8cedd8fa15877,0xa9f8d05cc0c37cdb,0x6342d403e9995fcc,0xcd1206bec9b26855,0x9c7d8a00045eb24d,0x9c63e4f9f6757a65,0x1b358d82afeb4,])),
                MNT6753Fr::new(BigInteger768([0x5c47dc04494f4bd2,0x9c673cd9289d41af,0x162259acba9d8d18,0x62cad4f296328097,0x8aaf9e1700b7c75d,0x55e78bf0544350b2,0x4f68ebcc4892c902,0xdab2889f96fa7b5b,0x2a03de10d75b9f18,0x1ea1e16fc08e4df6,0x6acecbff7d2f538,0x9435d0a83b56,])),
                MNT6753Fr::new(BigInteger768([0x57c48852f8169d69,0x770318c8f24e3ac0,0xa0305f4306f0fbf4,0xf24a6cdad69062c1,0x193310c1c542ab5e,0x34b6461663f4fe2a,0xe7a085a783023999,0xb5ce7b9c96faf8e0,0x7552f4cfa41a306a,0x2f174937af08a752,0x1a0cef0caa379120,0xaf994027adab,])),
            ],
            vec![
                MNT6753Fr::new(BigInteger768([0x9720001bc9352497,0x26db06d9f4127454,0x9cce839d50eab099,0xba25501620cf63a9,0x795125f6eb018f87,0x694e8cec73b544f8,0xdb77a066d8a2cdd5,0x7aabd5789a9eafe3,0x178cc6b3542ceaa6,0xa6ac0cd365b9c275,0x122759efe8da9356,0x8e1dde78adb9,])),
                MNT6753Fr::new(BigInteger768([0xa9c2b63431ec99e7,0xb05d41809af7e5dc,0x2cbd97c762aecb7,0x4d41c4687b6d4477,0x8381b288c0dbf80,0x50d30f6e9cd8073e,0xbd5d9a24ab8be9f5,0x53f6ff54d29bfaf6,0xdfcf47396745930f,0xf9624d429b121957,0x2eff2dd22352fa1c,0x8062baa0e970,])),
                MNT6753Fr::new(BigInteger768([0x686af5fafbfbf6ea,0x1e1c039393b53fbf,0x395bda15104e42d7,0x86bd133dc0ecd7de,0xe6edda60379dd98,0xa4b50608cd0cbda3,0x71914eaa21572,0x716fc727079df56d,0x92d198f1997ebcb0,0x2bc460bbd690afcc,0xed78f65c0b4e499e,0x2bfad26243bd,])),
            ]
        ];

        for i in 0..test_ins.len(){
            let mut state = test_ins[i].clone();
            MNT6PoseidonHash::poseidon_perm(&mut state);
            assert_eq!(state, test_outs[i],"Output does not match mnt6-753Fr, test vector {:?}", i);
        }
        
    }

    /*
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
*/

/*
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
*/
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
        // test vectors are computed on random input using the sage script of 
        // [IAIK](https://extgit.iaik.tugraz.at/krypto/hadeshash), commit 7ecf9a7d4f37e777ea27e4c4d379443151270563. 
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
        // test vectors are computed on random input using the sage script of 
        // [IAIK](https://extgit.iaik.tugraz.at/krypto/hadeshash), commit 7ecf9a7d4f37e777ea27e4c4d379443151270563. 
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
