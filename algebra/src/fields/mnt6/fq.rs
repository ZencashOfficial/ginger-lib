//! Base field Fq for the MNT6-298, a 298 bit prime field with duadicity 34.

use crate::{
    biginteger::BigInteger320 as BigInteger,
    fields::{Fp320, Fp320Parameters, FpParameters},
};

pub type Fq = Fp320<FqParameters>;

pub struct FqParameters;

impl Fp320Parameters for FqParameters {}
impl FpParameters for FqParameters {
    type BigInt = BigInteger;

    /// MODULUS =
    /// 4759222861692613257533492496530484515451248785528235155532677357391646\
    /// 47307408490559963137
    const MODULUS: BigInteger = BigInteger([
        0xbb4334a400000001,
        0xfb494c07925d6ad3,
        0xcaeec9635cf44194,
        0xa266249da7b0548e,
        0x3bcf7bcd473,
    ]);

    const MODULUS_BITS: u32 = 298;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 22;

    /// Montgomery constant =
    /// 2233646483262814149388017053592230295549237255497924206830512748722002\
    /// 60503540791531766876
    const R: BigInteger = BigInteger([
        0xc3177aefffbb845c,
        0x9b80c702f9961788,
        0xc5df8dcdac70a85a,
        0x29184098647b5197,
        0x1c1223d33c3,
    ]);

    /// Montgomery constant squared =
    /// 1639831447225064468267151243689723805258943971272055777812343054963258\
    /// 61831001705438796139
    const R2: BigInteger = BigInteger([
        0x465a743c68e0596b,
        0x34f9102adb68371,
        0x4bbd6dcf1e3a8386,
        0x2ff00dced8e4b6d,
        0x149bb44a342,
    ]);

    const INV: u64 = 0xbb4334a3ffffffff;

    /// generator = 10
    const GENERATOR: BigInteger = BigInteger([
        0xb1ddfacffd532b94,
        0x25e295ff76674008,
        0x8f00647b48958d36,
        0x1159f37d4e0fddb2,
        0x2977770b3d1,
    ]);

    const TWO_ADICITY: u32 = 34;

    /// 2^34-th root of unity =
    /// 1206388178269131734587688294856900998453770080308916180101097729373635\
    /// 54409782252579816313
    const ROOT_OF_UNITY: BigInteger = BigInteger([
        0x818b361df1af7be4,
        0x2ae2750d46a53957,
        0x5784a8fe792c5f8a,
        0xf9bd39c0cdcf1bb6,
        0x6a24a0f8a8,
    ]);

    /// (q-1)/2 =
    /// 2379611430846306628766746248265242257725624392764117577766338678695823\
    /// 23653704245279981568
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xdda19a5200000000,
        0x7da4a603c92eb569,
        0x657764b1ae7a20ca,
        0xd133124ed3d82a47,
        0x1de7bde6a39,
    ]);

    /// T = (q-1)/2^duadicity =
    /// 2770232305450256248897344628657729199302411164115319933935928482906687\
    /// 1159442729
    const T: BigInteger = BigInteger([
        0xe4975ab4eed0cd29,
        0xd73d10653ed25301,
        0x69ec1523b2bbb258,
        0x3def351ce8998927,
        0xef,
    ]);

    /// (T-1)/2 =
    /// 1385116152725128124448672314328864599651205582057659966967964241453343\
    /// 5579721364
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0xf24bad5a77686694,
        0x6b9e88329f692980,
        0xb4f60a91d95dd92c,
        0x9ef79a8e744cc493,
        0x77,
    ]);
}
