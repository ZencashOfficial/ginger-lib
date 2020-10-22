use algebra::{PrimeField, FpParameters, PairingEngine};
use itertools::*;

use super::nvidia;
use super::utils;

// Instead of having a very large OpenCL program written for a specific curve, with a lot of
// rudandant codes (As OpenCL doesn't have generic types or templates), this module will dynamically
// generate OpenCL codes given different PrimeFields and curves.

static FFT_SRC: &str = include_str!("cl/fft.cl");
static FIELD2_SRC: &str = include_str!("cl/field2.cl");
static EC_SRC: &str = include_str!("cl/ec.cl");
static MULTIEXP_SRC: &str = include_str!("cl/multiexp.cl");
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
}

fn define_field<L: Limb>(name: &str, limbs: Vec<L>) -> String {
    format!(
        "#define {} ((FIELD){{ {{ {} }} }})",
        name,
        join(limbs.iter().map(|l| l.value()), ", ")
    )
}

/// Generates OpenCL constants and type definitions of prime-field `F`
fn params<F, L: Limb>() -> String
where
    F: PrimeField
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
fn field<F, L: Limb>(name: &str) -> String
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

fn field2(field2: &str, field: &str) -> String {
    String::from(FIELD2_SRC)
        .replace("FIELD2", field2)
        .replace("FIELD", field)
}

fn fft(field: &str) -> String {
    String::from(FFT_SRC).replace("FIELD", field)
}

fn ec(field: &str, point: &str) -> String {
    String::from(EC_SRC)
        .replace("FIELD", field)
        .replace("POINT", point)
}

fn multiexp(point: &str, exp: &str) -> String {
    String::from(MULTIEXP_SRC)
        .replace("POINT", point)
        .replace("EXPONENT", exp)
}

pub fn kernel_multiexp<E>(limb64: bool) -> String
where
    E: PairingEngine
{
    vec![
        if limb64 {
            field::<E::Fr, Limb64>("Fp")
        } else {
            field::<E::Fr, Limb32>("Fp")
        },
        if limb64 {
            field::<E::Fq, Limb64>("Fq")
        } else {
            field::<E::Fq, Limb32>("Fq")
        },
        ec("Fq", "G1"),
        multiexp("G1", "Fp"),
        field2("Fq2", "Fq"),
        ec("Fq2", "G2"),
        multiexp("G2", "Fp"),
    ]
    .join("\n\n")
}

pub fn kernel_fft<E: PairingEngine>(limb64: bool) -> String 
{
    vec![
        if limb64 {
            field::<E::Fr, Limb64>("Fp")
        } else {
            field::<E::Fr, Limb32>("Fp")
        },
        fft("Fp"),
    ]
    .join("\n\n")
}
