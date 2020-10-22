mod nvidia;
mod utils;

use algebra::{PrimeField, FpParameters};
use itertools::*;

static COMMON_SRC: &str = include_str!("cl/common.cl");
static FIELD_SRC: &str = include_str!("cl/field.cl");

pub trait Limb: Sized + Clone + Copy {
    type LimbType: Clone + std::fmt::Display;
    fn zero() -> Self;
    fn new(val: Self::LimbType) -> Self;
    fn value(&self) -> Self::LimbType;
    fn bits() -> usize;
    fn ptx_info() -> (&'static str, &'static str);
    fn opencl_type() -> &'static str;
    fn limbs_of<T>(value: T) -> Vec<Self> {
        utils::limbs_of::<T, Self::LimbType>(value)
            .into_iter()
            .map(|l| Self::new(l))
            .collect()
    }
    // /// Calculate the `INV` parameter of Montgomery reduction algorithm for 32/64bit limbs
    // /// * `a` - Is the first limb of modulus
    // fn calc_inv(a: Self) -> Self;
    // fn calculate_r2<F: PrimeField>() -> Vec<Self>;
}

#[derive(Clone, Copy)]
pub struct Limb32(u32);
impl Limb for Limb32 {
    type LimbType = u32;
    fn zero() -> Self {
        Self(0)
    }
    fn new(val: Self::LimbType) -> Self {
        Self(val)
    }
    fn value(&self) -> Self::LimbType {
        self.0
    }
    fn bits() -> usize {
        32
    }
    fn ptx_info() -> (&'static str, &'static str) {
        ("u32", "r")
    }
    fn opencl_type() -> &'static str {
        "uint"
    }
    // fn calc_inv(a: Self) -> Self {
    //     let mut inv = 1u32;
    //     for _ in 0..31 {
    //         inv = inv.wrapping_mul(inv);
    //         inv = inv.wrapping_mul(a.value());
    //     }
    //     Self(inv.wrapping_neg())
    // }
    // fn calculate_r2<F: PrimeField>() -> Vec<Self> {
    //     calculate_r2::<F>()
    //         .into_iter()
    //         .map(|l| Self::new(l))
    //         .collect()
    // }
}

#[derive(Clone, Copy)]
pub struct Limb64(u64);
impl Limb for Limb64 {
    type LimbType = u64;
    fn zero() -> Self {
        Self(0)
    }
    fn new(val: Self::LimbType) -> Self {
        Self(val)
    }
    fn value(&self) -> Self::LimbType {
        self.0
    }
    fn bits() -> usize {
        64
    }
    fn ptx_info() -> (&'static str, &'static str) {
        ("u64", "l")
    }
    fn opencl_type() -> &'static str {
        "ulong"
    }
    // fn calc_inv(a: Self) -> Self {
    //     let mut inv = 1u64;
    //     for _ in 0..63 {
    //         inv = inv.wrapping_mul(inv);
    //         inv = inv.wrapping_mul(a.value());
    //     }
    //     Self(inv.wrapping_neg())
    // }
    // fn calculate_r2<F: PrimeField>() -> Vec<Self> {
    //     calculate_r2::<F>()
    //         .into_iter()
    //         .tuples()
    //         .map(|(lo, hi)| Self::new(((hi as u64) << 32) + (lo as u64)))
    //         .collect()
    // }
}

fn define_field<L: Limb>(name: &str, limbs: Vec<L>) -> String {
    format!(
        "#define {} ((FIELD){{ {{ {} }} }})",
        name,
        join(limbs.iter().map(|l| l.value()), ", ")
    )
}

/// Calculates `R ^ 2 mod P` and returns the result as a vector of 32bit limbs
// fn calculate_r2<F: PrimeField>() -> Vec<u32> {
//     // R ^ 2 mod P
//     BigInteger::new(utils::limbs_of::<_, u32>(F::one()))
//         .modpow(
//             &BigUint::from_slice(&[2]),                          // ^ 2
//             &BigUint::new(utils::limbs_of::<_, u32>(F::char())), // mod P
//         )
//         .to_u32_digits()
// }

/// Generates OpenCL constants and type definitions of prime-field `F`
fn params<F, L: Limb>() -> String
where
    F: PrimeField,
{
    let one = L::limbs_of(F::one()); // Get Montgomery form of F::one()
    let p = L::limbs_of(F::Params::MODULUS); // Get regular form of field modulus
    let r2 = F::Params::R2;
    let limbs = one.len(); // Number of limbs
    let inv = F::Params::INV;
    let limb_def = format!("#define FIELD_limb {}", L::opencl_type());
    let limbs_def = format!("#define FIELD_LIMBS {}", limbs);
    let limb_bits_def = format!("#define FIELD_LIMB_BITS {}", L::bits());
    let p_def = define_field("FIELD_P", p);
    let r2_def = define_field("FIELD_R2", L::limbs_of(r2));
    let one_def = define_field("FIELD_ONE", one);
    let zero_def = define_field("FIELD_ZERO", vec![L::zero(); limbs]);
    let inv_def = format!("#define FIELD_INV {}", inv);
    let typedef = format!("typedef struct {{ FIELD_limb val[FIELD_LIMBS]; }} FIELD;");
    join(
        &[
            limb_def,
            limbs_def,
            limb_bits_def,
            one_def,
            p_def,
            r2_def,
            zero_def,
            inv_def,
            typedef,
        ],
        "\n",
    )
}

/// Returns OpenCL source-code of a ff::PrimeField with name `name`
/// Find details in README.md
pub fn field<F, L: Limb>(name: &str) -> String
where
    F: PrimeField,
{
    join(
        &[
            COMMON_SRC.to_string(),
            params::<F, L>(),
            nvidia::field_add_sub_nvidia::<F, L>(),
            String::from(FIELD_SRC),
        ],
        "\n",
    )
    .replace("FIELD", name)
}