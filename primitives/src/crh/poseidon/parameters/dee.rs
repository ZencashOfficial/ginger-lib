use crate::crh::{
    PoseidonParameters,
    FieldBasedHashParameters, PoseidonHash, batched_crh::PoseidonBatchHash,
    PoseidonQuinticSBox,
};
use algebra::fields::tweedle::*;

use algebra::biginteger::BigInteger256 as BigInteger;
use algebra::field_new;

#[derive(Clone)]
/// x^5-POSEIDON-128 parameters for scalar field of the Tweedle Dee (= Fr).
/// The number of rounds are from their [paper](https://eprint.iacr.org/2019/458.pdf), the 
/// round constants and MDS matrix are generated using the sage script from 
/// [IAIK](https://extgit.iaik.tugraz.at/krypto/hadeshash), commit 7ecf9a7d4f37e777ea27e4c4d379443151270563. 
pub struct FrPoseidonParameters;

impl FieldBasedHashParameters for FrPoseidonParameters {
    type Fr = Fr;
    // Number of partial rounds
    const R: usize = 2;  // The rate of the hash function
}

// x^5-POSEIDON-128 parameters for the base field of the Tweedle Dum curve
impl PoseidonParameters for FrPoseidonParameters {

    const T: usize = 3; // Size of the internal state (in field elements)
    const R_F: i32 = 4; // Half number of full rounds (the R_f in the paper) 
    const R_P: i32 = 57; // Number of partial rounds

    // The zero element of the field
    const ZERO: Fr = field_new!(Fr, BigInteger([0x0, 0x0, 0x0, 0x0]));
    
    // The constant 3 to add to the position corresponding to the capacity (Montgomery rep.)
    // (we don't use the domain separator for now)
    const C2: Fr = field_new!(Fr, 
        BigInteger([
            0x123baa49fffffff5,
            0xd90b134e59472506,
            0xffffffffffffffff,
            0x3fffffffffffffff
        ])
    );

    // State vector after permutation of zero state vector (Montgomery rep.)
	const AFTER_ZERO_PERM: &'static [Fr] = &[
        field_new!(Fr,BigInteger([0x24474df3c64b6467,	0xc0c46698dadb69fb,	0x89479fdaa2d4a3d3,	0x1038e2503c010653])),
        field_new!(Fr,BigInteger([0xdbd0f46c618cbe3f,	0xd4c1353c7911ca30,	0xc7d9c14b932bea56,	0x3e8a462618bedbfc])),
        field_new!(Fr,BigInteger([0x24022ec1dbfc566b,	0xed45883272d605ba,	0x147a173707407b8b,	0xebc6f1edf47ed55])),
    ];

    // Array of round constants
    const ROUND_CST: &'static [Fr] = &[
        // Constants converted to Montgomery representation.
        // For rounds 4 + 57 + 4 = 65

        field_new!(Fr,BigInteger([0x96e49caa308157a9,0x21d4aac20df4ad70,0x5b9562555ec3ee9b,0x286e579e0cd44f71,])),
        field_new!(Fr,BigInteger([0x8691bb3441e3ccaa,0xbacb8e29165c8512,0x13476958dd2b0cba,0x2075c2c33a95bcf2,])),
        field_new!(Fr,BigInteger([0xedd3ab8aa3f6feea,0x84a682bd70dd3cde,0x8a991374ca8e51b7,0x3baa33e337193ebe,])),
        field_new!(Fr,BigInteger([0x2b70f1442552f235,0x971124f514cbbdfd,0xcd08c12001147d7c,0x1e8053682decb7fa,])),
        field_new!(Fr,BigInteger([0xa5d6b9d41de543ab,0x3f1bcb7206821903,0x29d5a582ea0304a3,0x25354b8ddbe6a074,])),
        field_new!(Fr,BigInteger([0xa9ebb39291e0c568,0x1418bbe311fd3534,0x492aba9b52b211b0,0x3f873b33a7918a1c,])),
        field_new!(Fr,BigInteger([0xbda9e7fc0a3098b4,0xb540794b4c450369,0xfabbe7491782408d,0x2bb53e236f2d2ce,])),
        field_new!(Fr,BigInteger([0x68bfba1545afed1c,0xc06d77120183cb06,0x657c32aa75dd1b85,0x339918b1bc346224,])),
        field_new!(Fr,BigInteger([0x9da841366cb6faa0,0x370c7bb69a87bafd,0x1d5517e4c392d1e7,0x3bf3521efc6ec45b,])),
        field_new!(Fr,BigInteger([0xef9f2f3969a8020d,0x7d20145243a59f7e,0x647b2f00e156d1aa,0x21b19e486cf0fadc,])),
        field_new!(Fr,BigInteger([0xf1cee9ce4f64cd1a,0x71e5d9cdb65202ba,0x47c7e723b990c2c8,0x122d8528b120c062,])),
        field_new!(Fr,BigInteger([0x191831c3a036efbe,0xb932d101e1e1b2b6,0xed839db6af04efc4,0x18898a60dfdb6e12,])),
        field_new!(Fr,BigInteger([0x1145b9b3b45f4b1e,0x4fe953f4e6a6403b,0xa9228be7b7078e4e,0x8b9dd5b754706ad,])),
        field_new!(Fr,BigInteger([0x404a331006c5410,0xada3015cbda68076,0xedd9092f3e7d3d6c,0x161604a28d333fe5,])),
        field_new!(Fr,BigInteger([0xe5ac55598707f3b3,0xcf97bf3f99b62394,0x182fd2e47be7abc8,0x209869291edb34b5,])),
        field_new!(Fr,BigInteger([0xc01347c4e578c6b2,0x1a9b2b8bc5ceb24,0x5f477f62b1224579,0x2876067a7f8f17f9,])),
        field_new!(Fr,BigInteger([0xb968124b4e2ee39b,0xf153cb6f82a92d14,0x92b9a472a4aa8635,0x398014eac5bd153c,])),
        field_new!(Fr,BigInteger([0x4c3424600afce50b,0xbb499350cd73769,0xa5aec8760597be85,0x25c0a6de06b98847,])),
        field_new!(Fr,BigInteger([0xc136d3ccd53c0c77,0x215e2499ebfa8287,0x131b6a6c67b6d34a,0x19ba276a822ba18e,])),
        field_new!(Fr,BigInteger([0xbcdfec243df8a880,0x69221f39732003,0x42575b3e08c64247,0x25a0df7c4b423cd1,])),
        field_new!(Fr,BigInteger([0x66e08079b80d3054,0xde4017fc4e2e3b13,0x197143d79dac0f2b,0x290cf56942717117,])),
        field_new!(Fr,BigInteger([0xb9d729b0aaa2ccda,0x13519e2ad9a7df52,0xd2df92c2d2da3692,0xa22764dd82d4402,])),
        field_new!(Fr,BigInteger([0x7d921ea4ed82ca74,0x369135e23f4f7bc6,0xf6c453b06fc2855b,0x1157000ae5a56b4f,])),
        field_new!(Fr,BigInteger([0x6a8cc7162724cab9,0xb6c2a377b7810ecf,0xcf1e38dd536e53e3,0xb3e75a6c26811d5,])),
        field_new!(Fr,BigInteger([0x1c34cc95619573bd,0xce4393fe0b40d306,0xa6b2a3a7b21cc93d,0x102fb0984de470a7,])),
        field_new!(Fr,BigInteger([0xa029cadc02dbee07,0x9685292bc549313e,0x39e7d4bc0a087cb5,0x33af821a108edd0f,])),
        field_new!(Fr,BigInteger([0x7ef8496c6960f53b,0x85244382d43b0296,0x108a5dedf31de1d7,0x2938595a354f3676,])),
        field_new!(Fr,BigInteger([0x5d26d35b2a167d7f,0xa3f116e23ba3c47,0x63ba716d8997040b,0x2d33e722e6c557fc,])),
        field_new!(Fr,BigInteger([0xf0010dd8ddbdfb90,0x11842ad3f8001f49,0xcfc955f574fe8eff,0x16a01e51c798fea0,])),
        field_new!(Fr,BigInteger([0xd29c0b3c67341f5d,0xde83cf468dab9fee,0x6e324b41959dc450,0xf7aef1e842ddb0e,])),
        field_new!(Fr,BigInteger([0xd75ca4547d514158,0xad7d2fba6b4ee87e,0x1ad69bb72a10b91e,0x10a3f1fbeead0098,])),
        field_new!(Fr,BigInteger([0xf38788a6a5ea48bc,0x7ec38ab17f9869f6,0x3b1e036ec22d9d88,0x2ccaaa619723aef6,])),
        field_new!(Fr,BigInteger([0xda1d874a92767272,0x222ea3352b315e9e,0x489a5e6e9d22b306,0x19f8adb0fb91466c,])),
        field_new!(Fr,BigInteger([0x8eb49dbaa74dda7d,0xecd2190329fe7c1c,0x471542a14378c306,0x197d940ac7c64169,])),
        field_new!(Fr,BigInteger([0xb796fcbed8e00da3,0x25b0f6b3c8cae312,0x3e82bbb0e61ed0b6,0x378f7066efe5df0f,])),
        field_new!(Fr,BigInteger([0x613ec5e01af5e481,0x5beba25e481e4647,0xac58c6699925efb2,0x38837de1266b70ad,])),
        field_new!(Fr,BigInteger([0xfb9213ce73dcced,0xfaf3a1fde67b84e8,0x55bb397b1fb9a303,0x1b2234020061d2c,])),
        field_new!(Fr,BigInteger([0x900ae19b639ae6fd,0x6ab8122a65aab019,0x527f0d2ea7034f0c,0x67137cd83dce5cd,])),
        field_new!(Fr,BigInteger([0xbd37fbac1201d052,0x9f63e6a4cc3f45fb,0xc7b3284494825589,0x3d0f06c82fe77b06,])),
        field_new!(Fr,BigInteger([0x3a9b075c4193ce4,0x7ec82421279c5249,0x25f1c0d11fcaa5d9,0x2efc6924decb5f2d,])),
        field_new!(Fr,BigInteger([0x2844d27b5d8de129,0xdf97e91391f14995,0x721fa02d5c3fd556,0x3b8a922fef852d4b,])),
        field_new!(Fr,BigInteger([0x1e44662a11f5ed5e,0xab3f6b3a94e3dd0a,0x2e91a457a6372784,0x357c328c34079547,])),
        field_new!(Fr,BigInteger([0xd413ad3a62fbbe64,0xb24307f80e9f452f,0x3b39d61dc4f6f7a2,0x3935b46dd131d22b,])),
        field_new!(Fr,BigInteger([0xb654772250350d11,0x375def9e13e75c7a,0x91a6259ef218944d,0x33849aa954c3a6a4,])),
        field_new!(Fr,BigInteger([0x555d3a2178b88d8,0xdf7ba7718237d4a6,0x117ee1c2e626feac,0x1d2b319d422c870e,])),
        field_new!(Fr,BigInteger([0x5c6619e0641c43c,0x6be9ccc5980464a,0x7232e1c064581be3,0x2120e0260a50c100,])),
        field_new!(Fr,BigInteger([0x105b891d0666ff14,0x605a0b332757b116,0x670f916438abbfde,0x241bbb4cb29c9611,])),
        field_new!(Fr,BigInteger([0xa42b5d2e72a4484,0x9a6904d3b7ef5e9b,0xeaa27186952c7100,0x232b8a7bd39cefca,])),
        field_new!(Fr,BigInteger([0xd2905a61a60e6041,0x899c33ec722f1492,0xf81d61172feb70f8,0x1adee8117b10487d,])),
        field_new!(Fr,BigInteger([0x1a767a4dfdb69578,0x2452a7893acfdb7f,0x2619366249cd5e98,0x1d0d1d5bf6280806,])),
        field_new!(Fr,BigInteger([0xb874658c0679b807,0xce00c33bd8c081c1,0xaa915801307affa4,0x2172f4b2cb61048d,])),
        field_new!(Fr,BigInteger([0xd24f87803febcd40,0x68b92da968f99600,0xd31f2cf8fd7a75ad,0xa7a57084a2256a9,])),
        field_new!(Fr,BigInteger([0xa3700e7bd137a1fe,0xb0c8ec5e5ebab893,0x4e415fa7835612dc,0x3f45a089f2cd1cc0,])),
        field_new!(Fr,BigInteger([0x56cb6f4891119f9b,0xb46e5f22e6e03b11,0x9f8e4199b2768762,0x741ac09033cf9bf,])),
        field_new!(Fr,BigInteger([0xb7ae92c42aabe136,0xd9ba5d820c085518,0xf8500c4b1587f2e5,0x3d0598e241e3faad,])),
        field_new!(Fr,BigInteger([0xafbeada7315e85ee,0xd17a95549ddc1a18,0xea61cdbe754ca3e9,0x27ddf4e85fbcca97,])),
        field_new!(Fr,BigInteger([0xad669b4ba8736840,0x530b30bccf439701,0xec116bdb1a3ab110,0x2e4cf052e809f53d,])),
        field_new!(Fr,BigInteger([0x55b4ab60c18c767e,0xb5c859060673ea39,0x69befcb83de80557,0x30b714a19bc19bb8,])),
        field_new!(Fr,BigInteger([0xdad50aa97c47bd65,0x567079a3ae6f3f7c,0xe61f26f49e92e610,0x35dacbd31703a3fe,])),
        field_new!(Fr,BigInteger([0x855e1755d971da61,0x6099a981e212b7f4,0x56cce7802d7719e8,0x1bd44f8356d4bd22,])),
        field_new!(Fr,BigInteger([0xdaa90fe0502851dc,0xf42dcd92b207147b,0x930fd0a2d55560ad,0x13d198adb65273c1,])),
        field_new!(Fr,BigInteger([0x2dbd2a02376fd004,0x5cdb36b191a00553,0x4241ba8364ebace,0x530a01feb911ca6,])),
        field_new!(Fr,BigInteger([0x130b4d244fdd1c37,0x500fbf60b9c856bb,0x28b11d1e2c32ca6,0xc0a6e6e36bf0a52,])),
        field_new!(Fr,BigInteger([0xaa3126dc4fc5274b,0xa62e7521570accb4,0x788d1f3ae57ef742,0x265b258fe90eedf1,])),
        field_new!(Fr,BigInteger([0x32036232210d984e,0xb9b2a8b58cd418c3,0x64dde61411c6cb07,0x39d0dde15138bd9e,])),
        field_new!(Fr,BigInteger([0xa1f6d09c7ae4fb46,0x5ca6b22f1c806c7,0x949e1fa886cd8f05,0x2de77c1f70569fa0,])),
        field_new!(Fr,BigInteger([0xc79c4ea9bf7df739,0xd7811101216dbc5f,0x6c5e0c2a16edd1b6,0x3f92d1a311740053,])),
        field_new!(Fr,BigInteger([0xa1239b9b6db0c5aa,0x9faa7ee688bf8d8,0xc72d4b2c98a3eca,0x6175c335ecb9411,])),
        field_new!(Fr,BigInteger([0x5a0806f3489ac665,0x3bf555c6f809326f,0x5dd894b982cccb3a,0x10c26b42184e2417,])),
        field_new!(Fr,BigInteger([0x49051c91d5d92cfc,0x2e3456653e338cf2,0xd51d669611799fff,0x96c20d9cded33cd,])),
        field_new!(Fr,BigInteger([0x1b0433c427e25ffb,0xda75da197898d2e7,0xa2bd435922b68c4d,0x1cb108347217f2d9,])),
        field_new!(Fr,BigInteger([0xf962ec7c80f243d6,0x4c526529b261357a,0x1b3d5b14924bc726,0x300e66d265163ca6,])),
        field_new!(Fr,BigInteger([0x695f0ca90ea0a921,0x3f3f144cc7273c9e,0x260fce9ffd26d074,0x81ba669eca7712b,])),
        field_new!(Fr,BigInteger([0x896bec269a095d02,0x40bb30acbe77311d,0x55ae996bb2b0b029,0x204882949c911712,])),
        field_new!(Fr,BigInteger([0x86d77b734624aca5,0x8a23044f00217dba,0x10eaf46f96c8e958,0xda0b77899936fcd,])),
        field_new!(Fr,BigInteger([0x4eeeed5416efc2b2,0x3fd96716a0651a14,0xd1582b838e00d405,0x238590a005a9261,])),
        field_new!(Fr,BigInteger([0xcda9def81a13694b,0xdc546ee0e40713aa,0x4061f5c54bcf4b3e,0x336fa4750bfc5144,])),
        field_new!(Fr,BigInteger([0xe21c95b18d495c15,0x33549ba45e7e8efd,0xe079e72c67b1e20d,0xa3c2dd5d8d0f8ea,])),
        field_new!(Fr,BigInteger([0x5c9b3ced1d4c6c0,0x472cbd2d64e5ed6d,0xc6ac7e1f9d20db9,0x2b9a24fdee775645,])),
        field_new!(Fr,BigInteger([0x368ec971c25526c7,0x82cea63463209096,0x9c93447d03605ba5,0x3db0d90439db352a,])),
        field_new!(Fr,BigInteger([0xa63883743b3fbc5d,0xf70843b35138358a,0x4207ee441e849f86,0x1976d24dae2f4402,])),
        field_new!(Fr,BigInteger([0x8931cb85aa9fbfad,0xaa06af2ae21b43e6,0xb88c735d5e459af,0x136fc39c52bac357,])),
        field_new!(Fr,BigInteger([0x3e4327133752d934,0xba77e0a3a0c9f17b,0x443322fd0825fa7a,0x3256f57fc43050bf,])),
        field_new!(Fr,BigInteger([0xd06dc3235761dd91,0xaf5a449dcfda30cd,0x55da6c31d05dbb62,0x3e8609d68aeab5d5,])),
        field_new!(Fr,BigInteger([0xb27f5076e9d516b4,0xeb8e90f978d1cb00,0xd5a81b105b6e4b73,0x22ecf60bfff720fb,])),
        field_new!(Fr,BigInteger([0xeee08420632f05ce,0x14adc2e637461f72,0x1a90be11356e17a1,0xd39db7c450214ca,])),
        field_new!(Fr,BigInteger([0x93ac880085863c4e,0xfae0316b7b12633,0xd2dd354b640a4178,0x39eab423b3a47b00,])),
        field_new!(Fr,BigInteger([0x101149e573669740,0x2a746ac13c1af657,0x5c2342b823dbf172,0x1b82a4d85630fd93,])),
        field_new!(Fr,BigInteger([0x1e031aaa5c6ea175,0xad92739acddd18d6,0x2a0a29c32e527994,0x1fad14f82f0e5345,])),
        field_new!(Fr,BigInteger([0x80cd44b05d072568,0xb2148dbcf15ecdde,0xd3094a3608d648f0,0xb1246becd42e7fd,])),
        field_new!(Fr,BigInteger([0xd07a2e2aa3d33846,0xb3b3fc9ad13e13c1,0x9c4364b647b23648,0x2a414e9a6da982e7,])),
        field_new!(Fr,BigInteger([0x8c5292db87bc1b,0x1b649a433757254b,0xf985f4efe7f634aa,0x11f4e20ebc378dd0,])),
        field_new!(Fr,BigInteger([0xfbbbd4b03a68b6d3,0xe886bf17ad7f1239,0x48defa4aad5e4c04,0x23e60b7a2e1170b,])),
        field_new!(Fr,BigInteger([0xdbbc56302b8ca0cd,0x1b7e59b8d8a73284,0x5a0ccff7360b3cce,0x82230ed52f4fa10,])),
        field_new!(Fr,BigInteger([0x67cfd64043031894,0x137b46d4968fccb6,0xe4f37461f5d1069,0x11b694bfe564cd0d,])),
        field_new!(Fr,BigInteger([0x6755d0d450c39ad6,0x1714d1889f9d2786,0xdcd1ecff8f7e422,0x27ebc7dd627ed1ce,])),
        field_new!(Fr,BigInteger([0xaf7d237edcb60760,0x1a00f166e283ff8c,0x19eddfe995403334,0x375fd0144fbb607b,])),
        field_new!(Fr,BigInteger([0xafda41207449ffeb,0x8fbf2d942e1cd5e6,0x6227f4e6af3f2286,0x2e5759c1b4aadd01,])),
        field_new!(Fr,BigInteger([0xd5c84258dcfae245,0xfc975a9ab0045a0,0xc73dea0cd4da9119,0x3b39f258da612976,])),
        field_new!(Fr,BigInteger([0x5e081bd9a682e19e,0x66a636b4c0123230,0xae33c049f0362f35,0x1f688804dae84786,])),
        field_new!(Fr,BigInteger([0x5ea84d6ac6e027f4,0x2ea4a6f786906232,0x4625e748b9fc017c,0x20b5aaf8bea72843,])),
        field_new!(Fr,BigInteger([0xcc10fa6370ffacb6,0xc6b48ad178f98813,0x521ce75207dc5f25,0x29fe5cb47289a1a5,])),
        field_new!(Fr,BigInteger([0xab78b36e5f109719,0x331f9de80dd059c1,0x6410238154b8f53e,0x33a9382ca67d6b3f,])),
        field_new!(Fr,BigInteger([0x3aa7aed95c6acf12,0xce82d756bdc55654,0x8314ece42b451cc8,0x39ea252bffe4e726,])),
        field_new!(Fr,BigInteger([0x5c42595b83a52a07,0x6764f443b127b061,0xb6a8ce661c7d79d2,0x196ba8fb7bc1d357,])),
        field_new!(Fr,BigInteger([0x44f990b1c283f601,0xb855dbb120223869,0xa2f87305ba4cc315,0x37e61e1c7498b704,])),
        field_new!(Fr,BigInteger([0x4981cafdbdbb0504,0x7928f57629ae380e,0x958bbaafc84e7f2f,0xccf700c62759cca,])),
        field_new!(Fr,BigInteger([0x4504b88f61180e5c,0x3df383c13a15c799,0x23a36dd5b19825e7,0x2c519343f0208c8f,])),
        field_new!(Fr,BigInteger([0x3d447daf8aa99e64,0x51f19b77b7da1240,0x9e7a3809940974a7,0x31c7652e9fedf9f2,])),
        field_new!(Fr,BigInteger([0xab6c4e1a360600ad,0x6615d0f0d3093416,0xfa2c6503ef5462c5,0x452bf0b082a2c1f,])),
        field_new!(Fr,BigInteger([0x2501ff61bb8b354f,0xce69a275a3fe2995,0x8f6508f81aee2aa1,0x38c3cc66db640b2a,])),
        field_new!(Fr,BigInteger([0x1011733bd844831c,0x5191bc470a159382,0xecb4671c62ea27ce,0x1ef18abfc0624f99,])),
        field_new!(Fr,BigInteger([0xce92cea0b5cc865b,0xb21e2625565ec764,0x960fb408279c1c59,0x3a6dbe54faca500a,])),
        field_new!(Fr,BigInteger([0xf6bf668ee5bbeb9e,0xf26bae754d2c82e5,0xfb8aff87ae3e7577,0x143a73249eb4ef97,])),
        field_new!(Fr,BigInteger([0xb3ee56543b5911f7,0x9416987baec46d04,0x827b1f7bec722736,0x23d95b715f02e4e3,])),
        field_new!(Fr,BigInteger([0x17cd1bb3c3cd71bc,0x34a4edfdb49664d5,0x1dfe05d354ff8e1,0x2b5666dbf888b8e4,])),
        field_new!(Fr,BigInteger([0x78691a6d1cdead8b,0x6ba2d5a4181e7324,0x3f42a001d144ed50,0x8b0b83cdb86fd14,])),
        field_new!(Fr,BigInteger([0xc6609eeb12f476f4,0x272fea94200707f1,0xee82f99f18b24de1,0x2d14a145f8cd2b3d,])),
        field_new!(Fr,BigInteger([0x877c81bb6e9acfff,0xed2a378b04e12161,0x20bf8568ab4afaf1,0x3e51738a9f35acaa,])),
        field_new!(Fr,BigInteger([0x3abace9514d62136,0x608721fc26647308,0xdceb23572b0326d9,0x176eefec1312a6e1,])),
        field_new!(Fr,BigInteger([0x55e9c5f0b6ebd40d,0x5613bc78dee421f1,0x726cc6c1a4c151f,0x44387b22f78125d,])),
        field_new!(Fr,BigInteger([0x8cf21c8c11c50a19,0xf38d3ed86c63ea92,0x9722d6a7e715fa5c,0x1cefb36d57333b0b,])),
        field_new!(Fr,BigInteger([0xf8514c5e99ee3478,0x136e39b0938b3fa9,0x2170db8f605584b,0x27f83ab788346fa4,])),
        field_new!(Fr,BigInteger([0xb67103e504bc8076,0x440531e18ce55446,0x73821c90ac72f6e7,0x14086e161728fcae,])),
        field_new!(Fr,BigInteger([0x4f7e22aae0fb3c0e,0x43611fc35a37f326,0x7cbe75be1e5adf08,0x26d2e93ee3620392,])),
        field_new!(Fr,BigInteger([0x960150a67b785caf,0x7dd6dcecf70783c8,0x8424edb16e1d8982,0x39dff2850ab656b3,])),
        field_new!(Fr,BigInteger([0xfd8c3571bfd2b8b7,0x72a251629e22c773,0x98258b4e6ab903d3,0x5a4da31df24c07e,])),
        field_new!(Fr,BigInteger([0x8100d081575fc412,0x9d945ce8df780d46,0xf8dffc51d21ba378,0x17653e02879a1c1a,])),
        field_new!(Fr,BigInteger([0x1cc6afbecd72eaec,0xd914c92bb4645d68,0xb567342d0b6e8abf,0x145b8302f3e6926f,])),
        field_new!(Fr,BigInteger([0x3ee4ae3b94d7d356,0xf720ad0bad443fdf,0x1f53bda3fcdb9d9e,0x3ff8c84933f9061b,])),
        field_new!(Fr,BigInteger([0x62d740d826a802a6,0xa23fadc95e66099d,0xb68a95ba48e2fe9f,0x1125c3cf1346a4a1,])),
        field_new!(Fr,BigInteger([0xeffbc5e4e68fa473,0x44aeb59551b40f22,0x2ad7773bd96d486d,0xa465d85529558aa,])),
        field_new!(Fr,BigInteger([0xa4ffb0da10092123,0x7e4664bbeefed955,0x2da4df97d7326cf,0x12284e4005acecf9,])),
        field_new!(Fr,BigInteger([0xe157045353d95aeb,0xebea343016236184,0xc40492d045edba28,0x3670ee52183c9f65,])),
        field_new!(Fr,BigInteger([0x56422f9dfe955f91,0xc8970385fd393f68,0x9e6c78fd0a0718a0,0xc9f876cad5b02d5,])),
        field_new!(Fr,BigInteger([0x5d983dbba66c6ab2,0xc605eb42e6228c49,0x44fddc327c07cd04,0x58626badc33f377,])),
        field_new!(Fr,BigInteger([0xd73325db296f0e10,0x6b048acaaa6945c6,0xbb70b8b860592cbd,0x376df75f1ad7a5d5,])),
        field_new!(Fr,BigInteger([0xcab5a5555d7407a6,0x978fbd3f5b6207ee,0xeb6f3dfe65efbd9c,0x1ace4baf623e7682,])),
        field_new!(Fr,BigInteger([0x92c2795998881a1d,0x7df6e9ecb27062ca,0x213dd0b3a1912a88,0x189add18730daf9f,])),
        field_new!(Fr,BigInteger([0x374649ecf6479719,0xb00c0d5e2a4832c0,0x8511cc4419b6d0b7,0x183e23d1c5775536,])),
        field_new!(Fr,BigInteger([0x34c43c68869522cc,0xbb9dca70287ce883,0x487bf4974c64a93a,0xb29634965016941,])),
        field_new!(Fr,BigInteger([0xbe20faacc419ee61,0xcb694e2620f38eb,0xa65db0e16a8fd201,0x3b54cee9b1496996,])),
        field_new!(Fr,BigInteger([0x207ebc3baba00ea9,0x5a294ff37a030241,0xeb94abdbea6fc5a9,0x2cd11336fcc76c56,])),
        field_new!(Fr,BigInteger([0xac7a7ca129b4fbcf,0x66042747e91d30d2,0x9f278a07ce224e74,0xd2d74d202a2d6be,])),
        field_new!(Fr,BigInteger([0x4431e3fd948bab0c,0x34aaaaa3e4f1156c,0x5d3c91752fe8b731,0x3c25783382f4ee7,])),
        field_new!(Fr,BigInteger([0x5ed8a0069c1d5512,0x7adad54826caeda,0xa4a4780d618b0f6d,0x23577a49569a2dfe,])),
        field_new!(Fr,BigInteger([0xaf0b6c314d870420,0x4bd784dcd8ad568e,0xf3c467f2ae003a42,0x363506fbe2bd0742,])),
        field_new!(Fr,BigInteger([0x602e5e50e2d11832,0x6e644029eb6e31c5,0x6dcd40893c47394a,0x6f722ac3c4b1a1d,])),
        field_new!(Fr,BigInteger([0x1c8a95e482ce33ab,0x526712e384d258,0xcb6ba047618df75,0x230805ec67f3cdc2,])),
        field_new!(Fr,BigInteger([0xad5460a73f3167d,0x8fae9d189d4fe84c,0x16f3824fa6c243ff,0x3a6c1c9ea7edca1f,])),
        field_new!(Fr,BigInteger([0x3e8f48039d677fa5,0xe356ebaba639943f,0x8d7a67ef2eaef811,0x37402b2c8c4986c7,])),
        field_new!(Fr,BigInteger([0xd838b890d6b98fa0,0xaa9b690ecc6eec54,0x7de2d92e79039631,0xb10fd90bf81c6dc,])),
        field_new!(Fr,BigInteger([0xedb5fa5b4f66b76,0x2773e55eba73c4bd,0x15be6f68ba792ff,0x2ef1fe836507366,])),
        field_new!(Fr,BigInteger([0xccf519d3508ac590,0x62fac0ba2e1ed280,0xa41a647e6e9d3ea6,0x27def5b2b6cf50a4,])),
        field_new!(Fr,BigInteger([0x5582eedf724c06eb,0x3843b4c7d396b8a9,0xffea221e843d91b8,0xb4762567f45a87d,])),
        field_new!(Fr,BigInteger([0xf885ae06faaf9648,0xf8ce25b0ff85f44d,0x680aeffb1bba63fe,0x5e3bb3c50c128fb,])),
        field_new!(Fr,BigInteger([0x6c3a2f5f796aa8b,0x8a007c7f67aeb28b,0x35b5e098dba3501c,0x146d766947af2758,])),
        field_new!(Fr,BigInteger([0xf27781d007dba27f,0xb7aacc02602ce2b,0xfd46e418c869f0b7,0x1acb4e3c024d450d,])),
        field_new!(Fr,BigInteger([0x1c2bdd990aaa89c3,0xd1e957570dabf355,0xaf41ffff2d077df,0x382702f73e8ca069,])),
        field_new!(Fr,BigInteger([0x64519cffed69bfa,0xac4e33014da13e00,0xba9dd8c3cb7a8172,0x46677549158c4e6,])),
        field_new!(Fr,BigInteger([0x1f77bfde28853e75,0x7337b4678d97e349,0xeecc22a8d76d2157,0x3d1a1cfc83a82a9c,])),
        field_new!(Fr,BigInteger([0x980a14f78bbb1f4c,0x7ef4dd2287474ee6,0x4a1079953b1645a0,0x300a620f54a8c0c1,])),
        field_new!(Fr,BigInteger([0x2d188648e9e3e5b6,0x5cf33207d679671d,0x40883d89b1253195,0x10ee49a953f56dc7,])),
        field_new!(Fr,BigInteger([0x3adadd2101d8bd63,0x33df0eb31797df7e,0xfcf86a39b9a45d9f,0x4d291c2cf0e0491,])),
        field_new!(Fr,BigInteger([0x61eea72d80ef6eaf,0x6f481ab2a00dab2d,0x33ecf8632d4671eb,0x263695fc8374dec9,])),
        field_new!(Fr,BigInteger([0x7bde72466887140f,0x772cbfdc28930ef0,0xc7f5244132bad5f8,0x337b03e80dfc68b4,])),
        field_new!(Fr,BigInteger([0x8d27d67f27a7f00c,0x3c11f7b6b49c0ae5,0x84309da69ff47c92,0x3435609d126dc6b7,])),
        field_new!(Fr,BigInteger([0xa837402105bbac05,0x4f8dec27d57822a,0x75845e3f571d59fd,0x7217a73f6d0b87f,])),
        field_new!(Fr,BigInteger([0x4928496226507efe,0x69ee20bec63aab0d,0x406ece9ebaced00,0x12c4af067c95dd18,])),
        field_new!(Fr,BigInteger([0x22780460c27d721c,0x96e7e3e52a3edc51,0xb4cf2754800b8ed9,0x186ce2256b772536,])),
        field_new!(Fr,BigInteger([0xa4dc31ccdd8acb7a,0xf854811e84e5dd10,0x523a8a9c60e9adc3,0x2723d921b8da0202,])),
        field_new!(Fr,BigInteger([0x1fd8951398b93603,0x21516dc03db96d8a,0x693c6b81b1b87981,0x205502e47f244d0a,])),
        field_new!(Fr,BigInteger([0x16f89f034dcec8c,0x1ae30f848f59a2b9,0xe40c4d57d16f7f8d,0x2c3feed6d9cfc99d,])),
        field_new!(Fr,BigInteger([0xff5010bf79fef123,0xd5f6a0543cff07cb,0xf1d1843b9ba397f3,0x1c53734dd27d1d19,])),
        field_new!(Fr,BigInteger([0x4678c6d2f5240880,0x82e87e24f7928ae4,0xf9d239f64491808d,0x2827c2883e78289f,])),
        field_new!(Fr,BigInteger([0x3c16f6844e77ba39,0xca6627dde714f4f1,0x3a7e05abeefcd616,0xc1db7a7053d3d11,])),
        field_new!(Fr,BigInteger([0x940fbb0d8f15fc89,0xa83ef0a29146e284,0xf334033e00da252e,0x2a753e83505367b4,])),
        field_new!(Fr,BigInteger([0x579ae7dbb4b5543d,0xc11e5d69c1f37c5f,0x108825316f0c826d,0x29550573a6a38c40,])),
        field_new!(Fr,BigInteger([0x72384a0e94db0bc9,0x56ad9e5623143221,0x4f507a6bb04490e9,0x1cfd98357916ba9a,])),
        field_new!(Fr,BigInteger([0xf5753093010cd994,0x8bbfe602479d124,0x95e1d83b2df4ca23,0x962b457c960f8,])),
        field_new!(Fr,BigInteger([0x2cb94df89ef7190,0xdec1268a034738d3,0x111aa8f0cf21138f,0x1357cd70bf53613,])),
        field_new!(Fr,BigInteger([0x649d8110f1bca7c4,0x19c8c23243063c88,0x7fb6592a1dfb3844,0x1734f127e83eaa59,])),
        field_new!(Fr,BigInteger([0x7fa139528dfe1917,0x3e6386d912812445,0xa0a814eeb76159d4,0x210f408101e0cdd5,])),
        field_new!(Fr,BigInteger([0xda9314c1877d416e,0x19190aed09448753,0x7b29dc10a8119a6a,0x2c2ded01c439129d,])),
        field_new!(Fr,BigInteger([0x50dc934c3e37a923,0xf23ecf566619ee89,0x3e9e384d163fc9b0,0x413a4b63cc56c1b,])),
        field_new!(Fr,BigInteger([0x9c9f08ca3f89903a,0xfedd50fbe2535acf,0xcfe0d443be3bc02a,0x1cd9a4eeff2084cd,])),
        field_new!(Fr,BigInteger([0x92307e46f2746b4c,0x8a57fcb7d070928c,0x2b88c02deab71c2f,0xb47ca4c3aa8b0c0,])),
        field_new!(Fr,BigInteger([0x6c0361dc084910e8,0xfe62117ae04f9db3,0x31c79368b87af81f,0x8619840e876ffb2,])),
        field_new!(Fr,BigInteger([0x3054afa48624292b,0xa7de75dd7f7bf8ab,0x5a9e215f8d2416,0x3fd3bd5b55b562bb,])),
        field_new!(Fr,BigInteger([0xa55e01b5a3143758,0x8851da79e5a81b82,0xd927c1f66c043cc1,0x245f2b5cb696f1cd,])),
        field_new!(Fr,BigInteger([0x2fb07f7df91845b1,0xdc65b7c3ceb2fdd3,0x2919b5d3bd0cde3f,0x32c6f61fbd7faadc,])),
        field_new!(Fr,BigInteger([0xdc9f298ce93858d5,0xa31ab8dc15b782e1,0xc2deed052375743,0x955ccce101132f4,])),
        field_new!(Fr,BigInteger([0xe8ee7e239d4f3c8a,0xe66900e9a576d4c7,0x2ffa52d554a07a06,0x31caad3af941507,])),
        field_new!(Fr,BigInteger([0x90fabb1e77dd7d96,0x7e8e130e429f6c92,0xec9834846267bc6b,0x1cce425bc75f31d5,])),
        field_new!(Fr,BigInteger([0xd0936b8bc86a9512,0x231fb65967f4990,0x35535748824cb029,0x2ace0a4fb270b8ff,])),
    ];

    // The MDS matrix constants
    const MDS_CST: &'static [Fr] = &[
        // Constants in Montgomery representation 
        field_new!(Fr,BigInteger([0xc5eb49398ac1be20,0x3f4083c5acbc9e35,0x130ea244d8cce92c,0x1d906cae17b5d067,])),
        field_new!(Fr,BigInteger([0x64a9256c5fd5c30c,0x70f6b8a30e416cd9,0xee873749a224c657,0x1c30ec990bae03ec,])),
        field_new!(Fr,BigInteger([0x4984b2df88757eb1,0x9a0fe48abd25137a,0xa55ed28b545527f5,0x11553250cedfd33d,])),
        field_new!(Fr,BigInteger([0x76bcc431379ac589,0xf76ec3d9150c6172,0x5348edf1bde3498a,0x1e173d4e4ffcd67d,])),
        field_new!(Fr,BigInteger([0x316e6cb528fdb7b2,0x4c1d99e919912b23,0x62c55f4f280478,0x18e7b300db5b6a06,])),
        field_new!(Fr,BigInteger([0x85b2566725a5269,0xe6afd4bdf8c46e92,0xdadcd5554f21ff18,0x24c0542db0effe16,])),
        field_new!(Fr,BigInteger([0x4a28d193d74ca6d1,0x6245d63b2da60541,0xd0c4a6f0a331650e,0x2ec987a777c27fab,])),
        field_new!(Fr,BigInteger([0x75c3d5740cc466bc,0x3c77a0736a1a5783,0xfe921de3ef52dcde,0x317a58d045c1dc1a,])),
        field_new!(Fr,BigInteger([0xc4e5cf8a59da9913,0xab22581d5f570f7e,0x95951970d2bf677a,0x1d596246dacf0119,])),
    ];
}

pub type FrQuinticSbox = PoseidonQuinticSBox<Fr, FrPoseidonParameters>;
pub type FrPoseidonHash = PoseidonHash<Fr, FrPoseidonParameters, FrQuinticSbox>;
pub type FrBatchPoseidonHash = PoseidonBatchHash<Fr, FrPoseidonParameters, FrQuinticSbox>;
