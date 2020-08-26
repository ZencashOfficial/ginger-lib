//! Exponent field of the MNT6-298, a 298 bit prime field with duadicity 17.

use crate::{
    biginteger::BigInteger320 as BigInteger,
    fields::{Fp320, Fp320Parameters, FpParameters},
};

pub type Fr = Fp320<FrParameters>;

pub struct FrParameters;

impl Fp320Parameters for FrParameters {}
impl FpParameters for FrParameters {
    type BigInt = BigInteger;

    /// MODULUS =
    /// 4759222861692613257533492496530484515451248792426947253955551285762102\
    /// 62817955800483758081
    const MODULUS: BigInteger = BigInteger([
        14487189785281953793u64,
        4731562877756902930u64,
        14622846468719063274u64,
        11702080941310629006u64,
        4110145082483u64,
    ]);

    const MODULUS_BITS: u32 = 298;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 22;

    /// Montgomery constant =
    /// 2233646483262814149388017053592230295518275145728992506352593304452273\
    /// 73121871807443600476
    const R: BigInteger = BigInteger([
        1784298994435064924u64,
        16852041090100268533u64,
        14258261760832875328u64,
        2961187778261111191u64,
        1929014752195u64,
    ]);

    /// Montgomery constant squared =
    /// 2730004785232377209109816556011608606400831266272357197129806122962639\
    /// 66512828033847775776
    const R2: BigInteger = BigInteger([
        28619103704175136u64,
        11702218449377544339u64,
        7403203599591297249u64,
        2248105543421449339u64,
        2357678148148u64,
    ]);

    const INV: u64 = 12714121028002250751u64;

    /// generator = 17
    const GENERATOR: BigInteger = BigInteger([
        2709730703260633621u64,
        13556085429182073539u64,
        10903316137158576359u64,
        5319113788683590444u64,
        4022235209932u64,
    ]);

    const TWO_ADICITY: u32 = 17;

    /// 2^17-th root of unity =
    /// 2647062505718000807580693023696543055301256755212639760340548780175809\
    /// 02343339784464690243
    const ROOT_OF_UNITY: BigInteger = BigInteger([
        9821480371597472441u64,
        9468346035609379175u64,
        9963748368231707135u64,
        14865337659602750405u64,
        3984815592673u64,
    ]);

    /// T = (r-1)/2^duadicity =
    /// 3630998887399759870554727551674258816109656366292531779446068791017229\
    /// 177993437198515
    const T: BigInteger = BigInteger([
        0x70964866b2d38b3,
        0x987520d4f1af2890,
        0x2a47657764b1ae89,
        0x6a39d133124ed3d8,
        0x1de7bde,
    ]);

    /// (T-1)/2 =
    /// 1815499443699879935277363775837129408054828183146265889723034395508614\
    /// 588996718599257
    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x384b24335969c59,
        0xcc3a906a78d79448,
        0x1523b2bbb258d744,
        0x351ce899892769ec,
        0xef3def,
    ]);

    /// (r-1)/2 =
    /// 2379611430846306628766746248265242257725624396213473626977775642881051\
    /// 31408977900241879040
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        0x64866b2d38b30000,
        0x20d4f1af28900709,
        0x657764b1ae899875,
        0xd133124ed3d82a47,
        0x1de7bde6a39,
    ]);
}
