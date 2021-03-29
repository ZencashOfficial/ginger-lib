use crate::crh::{
    PoseidonParameters,
    FieldBasedHashParameters, PoseidonHash, batched_crh::PoseidonBatchHash,
    PoseidonQuinticSBox,
};
use algebra::fields::tweedle::*;

use algebra::biginteger::BigInteger256 as BigInteger;
use algebra::field_new;

#[derive(Clone)]
/// x^5-POSEIDON-128 parameters for scalar field of the Tweedle Dum (=the Tweedle Dee Fq ).
/// 
/// The number of rounds are computed by ./evidence/calc_round_numbers.py, round constants and matrix 
/// are generated using the script ./evidence/generate_parameters_grain.
pub struct FqPoseidonParameters;

impl FieldBasedHashParameters for FqPoseidonParameters {
    type Fr = Fq;
    const R:usize = 2;  // The rate of the hash function
}

impl PoseidonParameters for FqPoseidonParameters {

    const T:usize = 3;  // Size of the internal state (in field elements)
    const R_F:i32 = 4;  // Half number of full rounds (the R_f in the paper)
    const R_P:i32 = 57; // Number of partial rounds

    // The zero element of the field
    const ZERO:Fq = field_new!(Fq, BigInteger([0x0, 0x0, 0x0, 0x0]));
    // The constant 3 "to add to the position corresponding to the capacity" (Montgomery representation)
    // (we don't use domain separation for now)
    const C2:Fq = field_new!(Fq, 
        BigInteger([
            0x123baa49fffffff5,
            0xd90b134e59472506,
            0xffffffffffffffff,
            0x3fffffffffffffff
        ]));

    // State vector after permutation of zero state vector (Montgomery representation)
    const AFTER_ZERO_PERM: &'static [Fq] = &[
        field_new!(Fq,BigInteger([0xb46ad2d321d1cefd, 0xbb1203b20c8e17ab, 0x7413c99cee248276, 0x10f3026539fe9e4])),
        field_new!(Fq,BigInteger([0x3489f06408fb2c18, 0x668c8ed351f01494, 0x944306c7ff00efae, 0x1e58dc79b4fcae6c])),
        field_new!(Fq,BigInteger([0xad5e1e851822cc05, 0x84d9ce927f225195, 0xdd4ab88f002f38dc, 0x19f0f4588a618830])),
    ];

    // Array of round constants
    const ROUND_CST: &'static[Fq]  = &[
        // Constants in Montgomery representation. 
        // For rounds 4 + 57 + 4 = 65

        field_new!(Fq,BigInteger([0x37c47f15ee35be3c, 0x279b126d65063b38, 0x6246979a2333ffaa, 0x2c64f53b3f9ca9e2])),
        field_new!(Fq,BigInteger([0xd3be6dd8e168b95f, 0x53ca28c6a541f42c, 0xf37a81809f8d4327, 0x2d995cf3ee78ffd8])),
        field_new!(Fq,BigInteger([0x78d761183f904ecd, 0xe29ff6871cf0f602, 0x6b8a6121e35e2d09, 0x22a7cd9b62e652ec])),
        field_new!(Fq,BigInteger([0xe483a4d8a2eb834, 0x2a1ac032bab4c92e, 0x477988ed3f6779c2, 0x1b61ba7a64a8c7bf])),
        field_new!(Fq,BigInteger([0x3e33a63cfb09ca, 0xfb61f5a83d134e83, 0x211239006249d702, 0x2d7dade2ce5d0565])),
        field_new!(Fq,BigInteger([0x555a37b931de5bc5, 0xc02064d37014a2a1, 0x90545ed28ddbca4a, 0x1f5114ea4a5617ec])),
        field_new!(Fq,BigInteger([0x47f3cb38565dfe9b, 0xc362f62cd78a37a3, 0xc34c00c0aab60d2, 0x2a146baa8ff6bcd8])),
        field_new!(Fq,BigInteger([0xc5f883460b17ddf1, 0x99f874dbff0a1f60, 0xf8ff08f802c642d8, 0x272540131350d75f])),
        field_new!(Fq,BigInteger([0x7b99da8b4284e040, 0x836396cfbd8f5561, 0xa51a9fd03b8fffad, 0x240676ca8fbb894a])),
        field_new!(Fq,BigInteger([0xd500933f3e43c198, 0xb3ceb2da4b9c4c98, 0xcb377ed1172c948a, 0x1befa4c00cf0bd0c])),
        field_new!(Fq,BigInteger([0xdf5f6cf09b6d8117, 0x5b1d662707f57e6d, 0x5ba6bea3e1ff4f1a, 0x11830374f185631])),
        field_new!(Fq,BigInteger([0x7dc3c615c7a59bf, 0xf492d9726eab590, 0x844a17a856f06138, 0x42ec0c73cebf328])),
        field_new!(Fq,BigInteger([0xca262a60f67cba02, 0xb6b465d3ce9e44fe, 0x2d36c9028f367b27, 0x3fd51089b54a71c8])),
        field_new!(Fq,BigInteger([0xbad023e4a05ee1e0, 0x9cfba0eb8a6c4373, 0x87ad0af9fd6da545, 0x249aafab9d38de52])),
        field_new!(Fq,BigInteger([0xaf742ae2652d7a02, 0x801b9685d8f437b5, 0x474349b88d655466, 0x290a0843e130526e])),
        field_new!(Fq,BigInteger([0xcbc52b5681c4e9a0, 0xd4b8624f8ad2622f, 0x9cc2415def14079f, 0x34d11464439db5f])),
        field_new!(Fq,BigInteger([0x982e9600684d487, 0x698287282b9464b4, 0xba818103a608fac1, 0xfe2a8b421af4e9c])),
        field_new!(Fq,BigInteger([0xdd63e860047ea0be, 0x3e3e3d39cf2493cb, 0x54c476c328a16f10, 0x13643cddddbcd890])),
        field_new!(Fq,BigInteger([0x59dcdb21da0ba5fb, 0x41f2e3da66cac476, 0x7c9e45cda6d9f4b7, 0x97b30ef4ca916a7])),
        field_new!(Fq,BigInteger([0xac01bfa8a86eb37a, 0xea9c9803bc74b64e, 0xd2b435cb7fe77506, 0x3cb8ea215814e432])),
        field_new!(Fq,BigInteger([0xdee85b44a9080d71, 0xff6fa29ac10e049f, 0x99b92867cb750bbb, 0x202743de827b898b])),
        field_new!(Fq,BigInteger([0xba43ea6d70dc822a, 0x6751703b5f2eab57, 0x9c0afa1b4224146f, 0x3296a08b6197b989])),
        field_new!(Fq,BigInteger([0x973eba5a1cb7fcd3, 0xe1585747fbaf4400, 0x10d440e51feac8c6, 0x201b3bf0c9fef355])),
        field_new!(Fq,BigInteger([0xa3487ed23240d027, 0xfb9d4d21976db8ab, 0xc8dc947a14b3b9e, 0x1bf8e31dd78e6de4])),
        field_new!(Fq,BigInteger([0xf6385973aee992c8, 0x69bb76b8b2a0f999, 0x51a965aededaaa78, 0xeb6c91440e63322])),
        field_new!(Fq,BigInteger([0x576ce040d62b1f3b, 0x91667c8e98e06d89, 0x7d9ef93e9aa049d1, 0x2f1f894113228ef3])),
        field_new!(Fq,BigInteger([0x7fc219537ae1d440, 0x65e50086c8567fa, 0x3a7ba01831411ce7, 0x3e4275f58c2f22c6])),
        field_new!(Fq,BigInteger([0x7b6dbf21ec26ed39, 0xdeb78e6955b5c50, 0x476ebbe175601226, 0x2e7c55c3f49de3a8])),
        field_new!(Fq,BigInteger([0x318f4e473cb676bf, 0xc074a0082682a900, 0xbdf53045d2c7fd8a, 0x1cc918dd0d879ec8])),
        field_new!(Fq,BigInteger([0x1999dbc6b408875a, 0x18c171f5b7279dac, 0x3c5db097a0e81e6a, 0x8c75fa67a174246])),
        field_new!(Fq,BigInteger([0xc67854ee519a33c2, 0xe72ecd1c3ad0537f, 0x46599220340338f0, 0x22f7538e68570c91])),
        field_new!(Fq,BigInteger([0x48caa0da6b6f5550, 0x8156e674742f2aab, 0xacd35f83e273d346, 0x298e6ebc0b81b7d8])),
        field_new!(Fq,BigInteger([0x82f18c846d012f0b, 0xfc0bf89b040ea3fb, 0x76d53d2726cb82b0, 0x5a0b439ed3ad3e8])),
        field_new!(Fq,BigInteger([0xcd0f709a29f1028f, 0xe205302048f33097, 0x8255fdd3d5172d0d, 0x241b7c4dc48014ed])),
        field_new!(Fq,BigInteger([0x604a0aa77f14daf3, 0xe7b208bee2f09d41, 0xef1b363ac6a172c6, 0x222d513a3787272b])),
        field_new!(Fq,BigInteger([0x52adb5eebd64b6ce, 0x65205376a621ea8d, 0x953e498454bc252a, 0x394deee2de9502d5])),
        field_new!(Fq,BigInteger([0xba5344655aded649, 0xf9647b488a4b84d2, 0x2a4b3ddba4273765, 0x3ee571b547815c5e])),
        field_new!(Fq,BigInteger([0xb67cb695743bcf2c, 0xc67c66ea57cde2d5, 0x63fdeb8ca090c60c, 0x3bf6dac3264ef3a6])),
        field_new!(Fq,BigInteger([0xe9768ec22ee24d1a, 0x5e4cc9486f74e6f9, 0x7e1cf421212ae83f, 0x2f3d7c79e305d31e])),
        field_new!(Fq,BigInteger([0xe483aaf471284e31, 0x3cc0f4ef0d54f91f, 0x948bd4e66d5206b8, 0x2d435d65bcc3016])),
        field_new!(Fq,BigInteger([0xe1998c7efe966381, 0x31f9244fd2e277db, 0xd582ebe0724789c, 0x3038c00d222cffc0])),
        field_new!(Fq,BigInteger([0xf2a2393d019ffb50, 0x97f6daf52c9bea21, 0x744f86aedc173f6b, 0x375ec74dd9035cd2])),
        field_new!(Fq,BigInteger([0x555673da2792dfd0, 0xeb0f58a3a8aede9f, 0x43c00f24d3126625, 0x1b9b25f88fe1d0f4])),
        field_new!(Fq,BigInteger([0xb2032b7c279123a6, 0x92f6976f487764ae, 0x7345cb5636a8b854, 0x320e74c4236f456c])),
        field_new!(Fq,BigInteger([0x453918b0f7bacbe8, 0x1b7a8bd4f3511976, 0x322d99edd93c4ff5, 0x12ccbec785bee23e])),
        field_new!(Fq,BigInteger([0x484edc89843fb120, 0xde2151d1dcf90181, 0x4efe680f7b6a95a9, 0x2f6112ab9d0141a1])),
        field_new!(Fq,BigInteger([0xf2c3ad788dd78f71, 0x61f9659f82dc3f75, 0xcc896ce1eed5f962, 0x195101d4a3c8a0fd])),
        field_new!(Fq,BigInteger([0x530234a2eb06c55b, 0x3edf997bfa7d75a7, 0x1d9e686b0241cd0, 0xa2aa7e82b11f39e])),
        field_new!(Fq,BigInteger([0x39a428c3edfa7eee, 0x4f56c859f8ed7833, 0x4f5bfcbd5e2e4761, 0x373d5289531d33a4])),
        field_new!(Fq,BigInteger([0x66afa4b156a1165e, 0x700e4d197bf45bbf, 0x1493aff0e3db766d, 0x323711852989a3b7])),
        field_new!(Fq,BigInteger([0x3c9480fc5cef809d, 0x65952840bc86619d, 0xd8a0e217d30a2187, 0x33002107f52027a1])),
        field_new!(Fq,BigInteger([0xbcbe6aac57bcdaa0, 0x2fdc10d23cbb64ea, 0x3f15bd1ba2f9bb04, 0x256cd7038713927c])),
        field_new!(Fq,BigInteger([0xc98878f6b69bdf5, 0xec7773f50b9fa914, 0x432d1245813865d0, 0x22c0544be2887f69])),
        field_new!(Fq,BigInteger([0x4c7d842d9b5ee1de, 0xbaf1bd508bade6d6, 0xe8c54294f7ecff25, 0x2a7f6f00ec1e9889])),
        field_new!(Fq,BigInteger([0xbbac92fefe78d48c, 0x16f864154bcc65fd, 0x874b996cef268dbe, 0x26946cf23541ade1])),
        field_new!(Fq,BigInteger([0x1713b41ac19fe34c, 0x947f909f2feef8c0, 0xa3c6e0ca6c77d65a, 0x3af83cd8cb4a19b2])),
        field_new!(Fq,BigInteger([0xcfcb2dd93eefc6ac, 0xf22798955eb309f9, 0xcb4f17e46d236e0e, 0x34ff10f676da9353])),
        field_new!(Fq,BigInteger([0xe3b2af4da1d58c7c, 0x788f421e1f14de73, 0x2cc57a01171dd45c, 0xf26930f23c52f60])),
        field_new!(Fq,BigInteger([0x8710e7095d5df46a, 0x33f81827249600b0, 0x56165ce5f11c369b, 0x2da2f8c748a30b96])),
        field_new!(Fq,BigInteger([0x57307e915176a1fd, 0xc583c12ff0a390cd, 0xaa2f6b841f57b8e5, 0x1edc04ee505bc98c])),
        field_new!(Fq,BigInteger([0x32ca688722080946, 0xfd2fc52df6252d9b, 0x6efa75fe92f8a6f2, 0x1504971d34a71e18])),
        field_new!(Fq,BigInteger([0x12d319cbc0e9f013, 0xe1c8efac8137f491, 0xdca22c12ee4d8ebe, 0x3e6efa267631e35f])),
        field_new!(Fq,BigInteger([0x7a1eebd3f1d11e0f, 0x57d287e211e06c2a, 0xa3022927cb4df13d, 0x3d94efdc2bbbcfbf])),
        field_new!(Fq,BigInteger([0xd8848ce53c4fb4e, 0x81cb34f2e9cb141a, 0x54531e93d4d6834a, 0x2a9a163db160803f])),
        field_new!(Fq,BigInteger([0x48bfc74a5a1c4a5, 0xa754f30f2a718a02, 0xabd479cf007ae654, 0x162dcc3ba53b9a18])),
        field_new!(Fq,BigInteger([0xbc3d49120e44710a, 0x9fc870509495ac5d, 0x41e08a1ab5c794c3, 0x24f40d2c18fdd6a7])),
        field_new!(Fq,BigInteger([0xc19b4ffb5254a627, 0x7f94d3fe0c271b16, 0x53ac4e4b9e410709, 0x3ef5ae1bf950bde3])),
        field_new!(Fq,BigInteger([0x6448a1f5d6fcbd66, 0xc11bf0b68d610dd7, 0x9c87bf34bd3c5738, 0xa92e7cb3561ac13])),
        field_new!(Fq,BigInteger([0x4445dd4ef1185fee, 0x4b549bc9ccc3d1fe, 0xdceb14adbb12dd7b, 0xf58965d6aaf097a])),
        field_new!(Fq,BigInteger([0x56b49bf9c66a98ea, 0x5c9e8b9b912ec7bd, 0x88420748876337bf, 0x289710f8e75acc71])),
        field_new!(Fq,BigInteger([0x340eaf913331a2a6, 0xb442248283adac5a, 0xd7cacedcc64de3c4, 0x3aedb6c568f6279])),
        field_new!(Fq,BigInteger([0xd704db2dbc165890, 0xbba3218ee51c7e83, 0xb171eecb56b9c9f9, 0x2bde5f0b6b2122fe])),
        field_new!(Fq,BigInteger([0x129300c649adb744, 0xf3b883b819f1caa0, 0xcbb2e8d0f8a1eb49, 0xfbd7f1e7cdfd1db])),
        field_new!(Fq,BigInteger([0xe75374548b266074, 0x9c1040e62090efa9, 0x5e572c38ed7e826c, 0x20da8ced7563ac8e])),
        field_new!(Fq,BigInteger([0x96b116bb627fd06c, 0x3e25baefd405e542, 0x4e34433759e56865, 0x119a4f1598139c26])),
        field_new!(Fq,BigInteger([0x7fd864c14ff297a2, 0x9e0b7dc6644c4714, 0x3d30c020eb23f913, 0x1f4c0c8f92d1fbbf])),
        field_new!(Fq,BigInteger([0x1a9006783c864ad2, 0xe52fd0a29235bf3, 0x11d9e2e01c6d5344, 0x37a1d08255c5536a])),
        field_new!(Fq,BigInteger([0x4cb662a9da81d5ce, 0xe09569073d3b6c2, 0x55d6f3bdf439e088, 0x18be7e660b1db848])),
        field_new!(Fq,BigInteger([0x49b576e569a65790, 0x981a6ca4b78db651, 0x6557c237cbe8f144, 0x2d62f70b80d17a3f])),
        field_new!(Fq,BigInteger([0xb4f47f4e02c852d1, 0xb644698928274503, 0x89c3245307162fe8, 0x3be82c2c87bb41ec])),
        field_new!(Fq,BigInteger([0x845ddac3bba36f57, 0xe2d1e5b76643ce40, 0x2336d60885aca968, 0x30d3cbf9a0ef5798])),
        field_new!(Fq,BigInteger([0xdc1557ca63dfb23c, 0xa7452115adfec5d8, 0xd9ac8a93624affff, 0x788e6b5064467bb])),
        field_new!(Fq,BigInteger([0x626f1b51943907db, 0xf379552928bb5a97, 0x2645339d91d29326, 0x10da812ab184d1ba])),
        field_new!(Fq,BigInteger([0xa27ffc9a876c47ed, 0x99e8b5300193217b, 0xa0373f1281dbd6bb, 0x27c0bec64168b789])),
        field_new!(Fq,BigInteger([0x73aa46a45ec0c411, 0xbe99dc1df00c4779, 0xde531e1f19746b35, 0x39df13ee067ab7ee])),
        field_new!(Fq,BigInteger([0x3109897b91a8ef68, 0xda79962dd6a7d722, 0xb2ef8b466c8e94d2, 0x3ccca83b7a491cf1])),
        field_new!(Fq,BigInteger([0x79754f10ae6de24a, 0x1b888a3750ed4b44, 0xa98fa40d5a14fbc5, 0xde2d8d9f2737756])),
        field_new!(Fq,BigInteger([0x901d4ad175ff9322, 0x3a2308e07efab339, 0x3a340a586ec89195, 0x12ee4afe28473a32])),
        field_new!(Fq,BigInteger([0x36a91936232238f, 0x9493b340eacd1883, 0xad676c9cf335ec77, 0xd695f57882ac71a])),
        field_new!(Fq,BigInteger([0xd252ec763d1c65a0, 0x7d1b6d03475d8860, 0x94482d161b7641fe, 0x28a2817a2a4ab2a8])),
        field_new!(Fq,BigInteger([0x500da3a590f83f3a, 0xdf141eb26993f9ab, 0x530844f7dbcaeff3, 0x13fdb1f8792913a])),
        field_new!(Fq,BigInteger([0xd5bb03a593df55d0, 0xcec4f993d5c47924, 0x2695f531adccd7c4, 0x9518ee54d1321de])),
        field_new!(Fq,BigInteger([0xa08b9f3156a0df4b, 0xd76b9674143a5329, 0xef55f20bae531371, 0x14a4e60c4e197a5a])),
        field_new!(Fq,BigInteger([0x8a2f0430a6ef314, 0x803f31ccb9d1cab1, 0x47cc3f2af6d22d6f, 0x28dc20b3d7b71add])),
        field_new!(Fq,BigInteger([0x17c9f68e6957c78e, 0x281d561e5bb675d6, 0x10838444e6d4093e, 0x3fb5e224545e7af2])),
        field_new!(Fq,BigInteger([0x417c4f03f8a86b68, 0xcc328bb7687269e0, 0xb4593db0507b4370, 0xf79d82ffc82e0bb])),
        field_new!(Fq,BigInteger([0xe26aa2b282138b3d, 0xfc77883d52011bef, 0x6c81f5841ac9e4d0, 0x2e92606d8272c335])),
        field_new!(Fq,BigInteger([0x1f422d3abe24adde, 0xbe1c10059203da4e, 0x16efbecd8e58b030, 0x2a86c879607060ea])),
        field_new!(Fq,BigInteger([0x2b6165888638186f, 0x33f09ad78466ef85, 0x4b849d02ce9482bb, 0x13dfcbc234b70df5])),
        field_new!(Fq,BigInteger([0x5883885e9ab03eb3, 0x6f7de94bca428e10, 0x9ec2b32a48a41dc3, 0x23c19925b4a3363])),
        field_new!(Fq,BigInteger([0xd18622ff39fd18ff, 0x10c2ddbd54282040, 0x5f0e9949a3b47218, 0x1ae9fa54d74f5c6b])),
        field_new!(Fq,BigInteger([0x4190922b87a61d44, 0x6c52b9bbc4e8ef00, 0x5292d3b6d2a0a23f, 0x32c5855a1a14ca62])),
        field_new!(Fq,BigInteger([0x23fee448bc3bd58a, 0xaddccbd2151eec04, 0x366415cef9790c1c, 0x26feb57f47284c34])),
        field_new!(Fq,BigInteger([0x7852c98c8541f089, 0xc1415fb1a5916060, 0x3595efca60ece345, 0x29c3c478087acaa7])),
        field_new!(Fq,BigInteger([0x4acd8e459b04be89, 0x2057d5df139d5b83, 0x2455b505d27715a9, 0x917c6824fcd1da3])),
        field_new!(Fq,BigInteger([0xb1673ae31d16780c, 0x2e47075aea2b1bf8, 0xc728e5bb31b1f399, 0x11c7bdea0fdb41de])),
        field_new!(Fq,BigInteger([0xaf4a61dec5222701, 0xde433a6635cd3db0, 0x98e45d524ac93f56, 0x451c0f400202b6e])),
        field_new!(Fq,BigInteger([0x6e7d9a745bf95c4f, 0x6ab15b1a7776d14f, 0xf96cbf62fbbec405, 0x634c9f5a8eb2e5e])),
        field_new!(Fq,BigInteger([0x59b2b36421593a6e, 0x82ec83bf8fc7a5bd, 0xd41d2c483a3fe62d, 0x2de522e77fe2129a])),
        field_new!(Fq,BigInteger([0xac263186458f0f51, 0x148f1fc848a2d163, 0xed8e6a35aa5090fa, 0x13a14d576ec70fc6])),
        field_new!(Fq,BigInteger([0xb7127ba146dca430, 0x21f8c18b76c7c13b, 0x1853c5a7f6652f45, 0x19f8e93fdbf04b25])),
        field_new!(Fq,BigInteger([0x596137eddaf77349, 0xdde8e5b43b97f652, 0x331d508520b2e3cc, 0x95c3957cc056200])),
        field_new!(Fq,BigInteger([0xe6479b4f00aff3b7, 0x1f43a26e75b6be05, 0x734d71add4642899, 0x1fcfa7a698ada48b])),
        field_new!(Fq,BigInteger([0xdc96a78189da4c4a, 0x48841f195e52972b, 0xe07e489a9dbeec18, 0x923d62262e4f7ab])),
        field_new!(Fq,BigInteger([0x4c0a92cd76201b8d, 0x3cc8966a4df13a45, 0xccc4122430fd1f5a, 0x6ac7d7366d1b22c])),
        field_new!(Fq,BigInteger([0xe436717698bdc615, 0xe22b265ea7b734d3, 0x16f15e0cb340397e, 0x1cc85e729961b94a])),
        field_new!(Fq,BigInteger([0x5be53b714fc0af9, 0x8851e8b01febd320, 0x20068ebd4967b25b, 0x609ca12cdd035b4])),
        field_new!(Fq,BigInteger([0x22a03fa9a044a78c, 0xdedc6ee8f7c9e52e, 0x209630871946d3a7, 0x3d6e963640b6db19])),
        field_new!(Fq,BigInteger([0xaaffbba51defd448, 0x6f29ff0296f0f3b, 0x59d6173f1c6ff060, 0x16d506abba876d48])),
        field_new!(Fq,BigInteger([0x2597c36e76e7686e, 0xb14aec1624ab1368, 0xacb895b0e4fd20aa, 0x332145524e56bb01])),
        field_new!(Fq,BigInteger([0x2703c4ffc1137df6, 0xac96cd8948ffe4dd, 0x7238a7e1fb784a4b, 0x30c2212db3b6e301])),
        field_new!(Fq,BigInteger([0x7c6546ead240ac68, 0xd25f0483542dd5ad, 0xbb7f246c634bdf5, 0x1452629aa5b39182])),
        field_new!(Fq,BigInteger([0x8e5ad8d93613babd, 0xb9c05bba93cccbca, 0x788187c21a14442b, 0x12bcd4d179dcbabd])),
        field_new!(Fq,BigInteger([0xb129e62dbcb9a91b, 0x8fdbcf833ee6fd1, 0xab9ede88f2f06864, 0x9f38ec4936c6bde])),
        field_new!(Fq,BigInteger([0x9a00e749941e0d8c, 0x9a169d3e9ce0c139, 0x2312a1caf47d0176, 0xc0b91e85804302e])),
        field_new!(Fq,BigInteger([0xfae859804e61a0d5, 0xc263be077a3b68, 0x722a4cfccf939f28, 0x15bccb87c0766fd7])),
        field_new!(Fq,BigInteger([0xb0bee321d1e0250f, 0x125112684afd028, 0xd7878403e58b28f5, 0x2ac74133da661ea1])),
        field_new!(Fq,BigInteger([0x72f97858bd11f1a, 0xcd060ca95a1fbac7, 0x3a6cd258dc017065, 0xf602fd0d942898c])),
        field_new!(Fq,BigInteger([0x4f7b3a88d733b3d7, 0xfbb50091ae99eba3, 0xbdcf47452942fc35, 0x188e6d12f0ececac])),
        field_new!(Fq,BigInteger([0x6e58e3126bded79c, 0x72246f77c89c304d, 0xcb4f15d868f17327, 0x117dc182f0d839e6])),
        field_new!(Fq,BigInteger([0xbac9ce245243bbc5, 0xc02325940dab7145, 0xc9b3712e9c294935, 0x34ba755149ad8327])),
        field_new!(Fq,BigInteger([0x600cb41c54dbc0f2, 0xb5fb4ef2384d61d8, 0x7346ca6ae6d845c, 0x331fe71959b76375])),
        field_new!(Fq,BigInteger([0x7464ef960600747d, 0x6f40aed8cd8c2bb7, 0x877fff6dc61c246c, 0x3d7183fa10977d47])),
        field_new!(Fq,BigInteger([0x6ad71a86c358f6eb, 0x9e49529b3bc53e0c, 0x44cc4eea0fa0f744, 0x252f512074943963])),
        field_new!(Fq,BigInteger([0x3765ebc0e56eabb4, 0xd0c679c9e0c290d0, 0xb1146e2f8854e195, 0x273183affcde552a])),
        field_new!(Fq,BigInteger([0xbe13b9ea010db0c6, 0x84cbfc5c22097720, 0x8e5a00f9eb7d9887, 0xc7dd41188fb0faf])),
        field_new!(Fq,BigInteger([0xf91af7d5ad596e45, 0x5875334f2fd1a838, 0xcc829f3e4b839c48, 0x100dedfd1baecd8e])),
        field_new!(Fq,BigInteger([0x12439287ca69f742, 0x995f6534ca39c850, 0xdbfe45ce4c939f43, 0x1a2d6ade69f3d829])),
        field_new!(Fq,BigInteger([0x7d944cd7d391fb8a, 0x68dd32bc134507fa, 0xe8270eb80acb53a1, 0x2c91ee4c8942baaf])),
        field_new!(Fq,BigInteger([0xe3b1ad4263d24db7, 0xed87d6d5577b582e, 0x15da7d94fb2bf7fa, 0x246120ff6a47f57])),
        field_new!(Fq,BigInteger([0x5a6f12b77bd26f0c, 0xe618de6ba6cbf4ba, 0x79c64ae1402ac70e, 0x24f4d4e9b3fb5bab])),
        field_new!(Fq,BigInteger([0xf1b3ccfd5709d1a5, 0xeea324986e63bd24, 0x5eb3d12db87a96f2, 0x342b9ece5f9ff30f])),
        field_new!(Fq,BigInteger([0x2d65a87ab5280a76, 0x43a00ba359f820cd, 0x115ddf62abdb9527, 0x2d4378faaa603e20])),
        field_new!(Fq,BigInteger([0xf67f6aa46f8c8cce, 0xa47d94e1c43d555b, 0x87b3aa9cfe1c4e45, 0x5c32c9bbc1bf5c8])),
        field_new!(Fq,BigInteger([0x94b477bc97c4788d, 0x7b447e7c187f162e, 0xc3fa78647d1ff0ea, 0x3e300479330c40c3])),
        field_new!(Fq,BigInteger([0x885a137c8a104363, 0xd7fb83d9bd6bed52, 0xffaca15407701b1f, 0x12970c8d4cf62b53])),
        field_new!(Fq,BigInteger([0x575514a214546ef5, 0x8f1bd4bc29bb82f, 0x3720efc0848755be, 0x3a101aec2b33885])),
        field_new!(Fq,BigInteger([0x2ed8ab47952def78, 0x228e5a3b7aaf3392, 0x5978e38e91363d57, 0xc4f896022c767b8])),
        field_new!(Fq,BigInteger([0x4a136f844a1c7912, 0xb3519938197cbad5, 0x865f34275dabc7f2, 0x3c9caf853618bcb1])),
        field_new!(Fq,BigInteger([0xd3f35d88649ceadb, 0x6d2a3e0eec52c431, 0xbadd444fa9e62051, 0x2285816e9956d22a])),
        field_new!(Fq,BigInteger([0x1c082cf6d7b38838, 0x90c438d75116833c, 0xc419264150c90f0a, 0x2183f34500fa12aa])),
        field_new!(Fq,BigInteger([0xd28019ae88159d70, 0x820e1ea7149f558c, 0x89a177b677080fa8, 0x21763d85b5a52828])),
        field_new!(Fq,BigInteger([0xd527c80c19698028, 0x3bb04753a336e0e3, 0x431ec8ecced29219, 0x15d8c24aed6b23d7])),
        field_new!(Fq,BigInteger([0x6a3637ec41b3f5f4, 0xf374468f1c584235, 0x1959eb782e141d2e, 0x12e7d4ab87a5aa2b])),
        field_new!(Fq,BigInteger([0x1663428c6b627b0b, 0xf012b0931e56be93, 0x7395f71b2cad41b, 0xd60b8e33a904f0f])),
        field_new!(Fq,BigInteger([0xe3fe46faaf5b1ff0, 0xc8843ce9b13a0620, 0x9cca732243c481c4, 0x14a4038d3717a409])),
        field_new!(Fq,BigInteger([0x1b3610a8d98a23a, 0x1f4e37772be8a4e9, 0x4b70cffd4569a525, 0x164efcd4c43ba78c])),
        field_new!(Fq,BigInteger([0x9860dcd07a941fff, 0x2c0e00636a9419ce, 0x604f16280e26d19a, 0x1a0bb916eae56218])),
        field_new!(Fq,BigInteger([0xb136b61001c3185b, 0x4dbf8d8ae0f59eee, 0x33a0070a4f277c00, 0x19e7726a3e6982df])),
        field_new!(Fq,BigInteger([0x2cdca51d46e9f964, 0x4b804c1826038be8, 0x847a24d5bcefbab2, 0x17ff932a8b75358b])),
        field_new!(Fq,BigInteger([0x932b12398d3fe8e2, 0x98c29df2b4a39053, 0x4fefcfb6f8bc8c13, 0x126e3aa52f6c7470])),
        field_new!(Fq,BigInteger([0xc69da8084379e90f, 0x827ba880103ef228, 0x360c9f6ea5c39527, 0x1526e30de7ca9ebc])),
        field_new!(Fq,BigInteger([0x46199352011bce99, 0xce9074fabffe0d9a, 0x6e8615b78a40d2fe, 0x2358246f5de1adad])),
        field_new!(Fq,BigInteger([0xf8fe3f6492c0cd3b, 0x9edb501d30391354, 0x8b642cbfde6bc761, 0x2f4857c63e62f5b6])),
        field_new!(Fq,BigInteger([0x3e4b5cf92a34b084, 0xa38976062bfb1267, 0xcb0ce38afa78ae65, 0x19fd3eee90082403])),
        field_new!(Fq,BigInteger([0x268c7a8bd767ee6c, 0x1853d8a5b854a62f, 0x6e878c3e70db847f, 0x37152814b6947317])),
        field_new!(Fq,BigInteger([0x46e9eb9777a38662, 0x3eb7262b64ea301e, 0xe40cbbc7371dd6b6, 0x361dca667f332d9f])),
        field_new!(Fq,BigInteger([0xae87b6fceca32196, 0x98f3794b67813131, 0xf95160756c56abeb, 0x34a03e3a277adacc])),
        field_new!(Fq,BigInteger([0xa3fe1133b7d8ce48, 0xe8fc23f53f92e858, 0x460844ae61efbbf5, 0x30f6779730781ea9])),
        field_new!(Fq,BigInteger([0x40e150733d0f3fd1, 0xef41de3206bbffa2, 0xe670b6837bd77272, 0x2b8897dc3dd2a56c])),
        field_new!(Fq,BigInteger([0x18c0e2da234492e, 0x511372eafa7d007b, 0x46c58ddfab6e9008, 0x3ffc2becb19591f0])),
        field_new!(Fq,BigInteger([0x90aeacca790d0f11, 0xff1c1d52d8812f9d, 0x2e3a8ebfc5deb7f6, 0xd606c4cb061db5e])),
        field_new!(Fq,BigInteger([0x94f2d3bde8489823, 0xda901ac8f507c6c6, 0xd38e9d61559ad6de, 0x3c88376fdc506764])),
        field_new!(Fq,BigInteger([0xe44a0072e89dcfe0, 0x14e5079dcc979bcf, 0x5f6997ae94b5382, 0x1d4093c3e01d7d6f])),
        field_new!(Fq,BigInteger([0xe7ec48d052d482dd, 0x3db9f331bca3503e, 0x58bd73a55ca8cdaf, 0x30791dea442d7861])),
        field_new!(Fq,BigInteger([0xb9366145fa10f54d, 0x1d815fe76662ed04, 0x1a3dbe9c87beb120, 0x123e0b78fa8073e8])),
        field_new!(Fq,BigInteger([0xc03797e1453b83c1, 0x99b2de19d8c20893, 0xc1e1cc9970436615, 0x3b02d0e5efec7711])),
        field_new!(Fq,BigInteger([0x3d34b1927b75bc8c, 0xdee0debec72a30aa, 0x45961b90995135fe, 0x75b4694548e6ed])),
        field_new!(Fq,BigInteger([0x547cf37b29a0cc79, 0x969ccad45141e105, 0xee955dee7efead3b, 0x19f78e17050c696b])),
        field_new!(Fq,BigInteger([0x1e49017dfa4c686, 0x2850fa62dd741324, 0x7f4a94f01d3ae02d, 0x2713caa17c68f69e])),
        field_new!(Fq,BigInteger([0x561ce3d42ecf9b9e, 0xc81404ca6b99e904, 0x2913922f264840e0, 0x3ac3154b225227f4])),
        field_new!(Fq,BigInteger([0xf4a2e30fe1d64830, 0x9da528cba5add844, 0xb361c5bf9d8b63eb, 0xc1a760e07305808])),
        field_new!(Fq,BigInteger([0x9bbd8b64a8f4e96d, 0x92804f3fc1acb2c3, 0x3976b1d74bacab98, 0x4eeb45f9938ab2c])),
        field_new!(Fq,BigInteger([0xd5738917a4cac30c, 0xe8fe369f6a8fa67c, 0x5da1fc771087250b, 0x38135706f2791fae])),
        field_new!(Fq,BigInteger([0xc09988dd302471be, 0x7beb6e0aead4c4b6, 0x3a263709ae5aedf0, 0x3419bebd05eb99cb])),
        field_new!(Fq,BigInteger([0xd308a27acf6021b0, 0x9ebc3947a2cbe81e, 0x3eeb919b4f2550e4, 0x2110955caf926462])),
        field_new!(Fq,BigInteger([0x9e717877527a5fb2, 0x9712a53b9a37770b, 0x3357a9e169810e2f, 0x3c036e3ee03c8542])),
        field_new!(Fq,BigInteger([0x9eade5b84d50d25c, 0x4058a7bfeca96142, 0x1b673888bda18af8, 0x16cdd86faf56ff66])),
        field_new!(Fq,BigInteger([0xbd56d521201ca244, 0xd55878f139cb6277, 0x5e42f9b670cf28bd, 0x3583f72970167a32])),
        field_new!(Fq,BigInteger([0x618f808dd19095b6, 0x476ae84ba0589c32, 0x34ada966df9db652, 0x3fcfdc414cc67951])),
        field_new!(Fq,BigInteger([0xc9ae023b5cfe502d, 0x64a02e7469acfe69, 0xf2b0dec5804788aa, 0xdeb3f0cfa8ff42e])),
        field_new!(Fq,BigInteger([0x88b8abfc2346cb75, 0xc10e31b70140582f, 0x2d44e1671473e74b, 0xbf0e099bab9597c])),
        field_new!(Fq,BigInteger([0x5828892deedd9938, 0x2b197e22331663e6, 0xa9b48a99bdba9086, 0x38dbd8c43572c645])),
        field_new!(Fq,BigInteger([0x77da34cc6e33d806, 0xf7d7546ddda63c6f, 0xd18c6df80efc18b0, 0x31b1845bc8efabd1])),
        field_new!(Fq,BigInteger([0x4ccab249c164efd1, 0x1624025cc9e31516, 0x374ea5b5343b3456, 0x5882ce8f1b65933])),

    ];

    // The MDS matrix constants
    const MDS_CST: &'static [Fq] = &[
        // Constants in Montgomery representation 
        field_new!(Fq,BigInteger([0x139a14725017926c, 0x3dbc5dd06e12e2f6, 0x27d5413206033874, 0x3da6a2e8b4718a26])),
        field_new!(Fq,BigInteger([0xd8c3238d5671f1e, 0x2854b226105da9d8, 0x6959f138101c0d28, 0x17a06fe8aec59a37])),
        field_new!(Fq,BigInteger([0x9d0a838937ae621e, 0x9d26fa798a9eef9d, 0xc56d137dd07d9385, 0x676b826a04e312e])),
        field_new!(Fq,BigInteger([0x2b292e56344e1e00, 0x37723cf73b435564, 0xe9b9bd04fbf22c17, 0x188f24c829545ab6])),
        field_new!(Fq,BigInteger([0xb8cc6d7994033456, 0xb42279c9b422448, 0x66e68272120391e4, 0x1042f1ddd8c9fc10])),
        field_new!(Fq,BigInteger([0xad65dabbf17dd08a, 0x5c99493d45eceed5, 0x27258d64087cef3e, 0x30b89dda2caf24a7])),
        field_new!(Fq,BigInteger([0x58e55bef9dbf3918, 0xf0bcc8a49c36c57b, 0x927ccd2041e71f29, 0x1afc658173c9cd65])),
        field_new!(Fq,BigInteger([0xde184719c8d05a93, 0x84b9d6f74be1bc0, 0x94258de2771b4a58, 0x10562bb91046d8be])),
        field_new!(Fq,BigInteger([0xa7e184139e2aef63, 0x54de0f418695ba1c, 0x1a00369a1f8def52, 0x25bba3782054c3c9])),
    ];
}

pub type FqQuinticSbox = PoseidonQuinticSBox<Fq, FqPoseidonParameters>;
pub type FqPoseidonHash = PoseidonHash<Fq, FqPoseidonParameters, FqQuinticSbox>;
pub type FqBatchPoseidonHash = PoseidonBatchHash<Fq, FqPoseidonParameters, FqQuinticSbox>;