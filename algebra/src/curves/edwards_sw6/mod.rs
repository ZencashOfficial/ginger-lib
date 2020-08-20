//! [Zexe's second inner curve](https://eprint.iacr.org/2018/962.pdf) for the sw6, a twisted
//! Edwards curve with a prime order subgroup of 374 bit.
//!
//! Its security level is about 190 bit.

use crate::field_new;
use crate::{
    biginteger::BigInteger384 as BigInteger,
    curves::{
        models::{ModelParameters, TEModelParameters, MontgomeryModelParameters},
        twisted_edwards_extended::{GroupAffine, GroupProjective},
    },
    fields::edwards_sw6::{fq::Fq, fr::Fr},
};
use std::str::FromStr;

#[cfg(test)]
mod tests;

pub type EdwardsAffine = GroupAffine<EdwardsParameters>;
pub type EdwardsProjective = GroupProjective<EdwardsParameters>;

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct EdwardsParameters;

impl ModelParameters for EdwardsParameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl TEModelParameters for EdwardsParameters {
    /// COEFF_A = -1 =
    /// 2586644260129690940106527336948935335363935127549146605398842626667204\
    /// 68348340822774968888139573360124440321458176
    const COEFF_A: Fq = field_new!(Fq, BigInteger([
        9384023879812382873,
        14252412606051516495,
        9184438906438551565,
        11444845376683159689,
        8738795276227363922,
        81297770384137296,
    ]));

    /// COEFF_D = 79743
    const COEFF_D: Fq = field_new!(Fq, BigInteger([
        0x4669ffffff46a638,
        0xa56bbe0a7f9fae05,
        0x403b425466a710b4,
        0xf6648db6ea4e988b,
        0x74d51b5923d35a8d,
        0xf8ed90b17fe903,
    ]));

    /// COFACTOR = 8
    const COFACTOR: &'static [u64] = &[8];

    /// COFACTOR^(-1) mod r =
    /// 1212489496935792628174934689194813438451844591038662471278843170572544\
    /// 1736421489799867521238554906438478484045560
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger([
        7353538464571651976,
        2030910049503177537,
        16726103313845754033,
        1110650741117127777,
        5304838729792721053,
        4975067790294675,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) = (GENERATOR_X, GENERATOR_Y);

    type MontgomeryModelParameters = EdwardsParameters;

    /// Multiplication by `a` is just negation.
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        -*elem
    }
}

impl MontgomeryModelParameters for EdwardsParameters {
    /// COEFF_A =
    /// 9008362308427189103711687048774306798471008020953914968541414705532906\
    /// 3590616489392386084989619674926965747987765
    const COEFF_A: Fq = field_new!(Fq, BigInteger([
        7594254284108454966u64,
        14287343397973578077u64,
        6490358977072726023u64,
        8023375322051995268u64,
        8242802613686040715u64,
        100541941146122331u64,
    ]));
    /// COEFF_B =
    /// 1685808029286972029735358632071504655516834325453755108544701156113914\
    /// 04757724333382582803149953685197474573470410
    const COEFF_B: Fq = field_new!(Fq, BigInteger([
        11173793475516310780u64,
        14217481814129454913u64,
        11878518835804377107u64,
        14866315431314324110u64,
        9234787938768687129u64,
        62053599622152261u64,
    ]));

    type TEModelParameters = EdwardsParameters;
}

impl FromStr for EdwardsAffine {
    type Err = ();

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        s = s.trim();
        if s.is_empty() {
            return Err(());
        }
        if s.len() < 3 {
            return Err(());
        }
        if !(s.starts_with('(') && s.ends_with(')')) {
            return Err(());
        }
        let mut point = Vec::new();
        for substr in s.split(|c| c == '(' || c == ')' || c == ',' || c == ' ') {
            if !substr.is_empty() {
                point.push(Fq::from_str(substr)?);
            }
        }
        if point.len() != 2 {
            return Err(());
        }
        let point = EdwardsAffine::new(point[0], point[1]);

        if !point.is_on_curve() {
            Err(())
        } else {
            Ok(point)
        }
    }
}

/// GENERATOR_X =
/// 1747017723244855069416909035124235519982943529688336599600423627426848\
/// 69862495746426366187462669992073196420267127
const GENERATOR_X: Fq = field_new!(Fq, BigInteger([
    3737364149926089590,
    13002967008679663837,
    9954144214462864555,
    3365719140389487049,
    8643066672427471196,
    120355578793479865,
]));

/// GENERATOR_Y =
/// 2084872000522588454953403744515407754454084396549301913240116355601425\
/// 23886549663106522691296420655144190624954833
const GENERATOR_Y: Fq = field_new!(Fq, BigInteger([
    6027299446526298157,
    12854429557810467099,
    11207279014226687864,
    17040621363687352702,
    6112671509202865855,
    44040319652922447,
]));
