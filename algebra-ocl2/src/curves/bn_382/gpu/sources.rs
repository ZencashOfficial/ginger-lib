use algebra::fields::bn_382::{Fq, Fr};

use crate::ffgen;
use crate::gpu::sources::{ec, field2, fft, multiexp};

pub fn kernel(limb64: bool) -> String {

    vec![
        if limb64 {
            ffgen::field::<Fr, ffgen::Limb64>("Fr")
        } else {
            ffgen::field::<Fr, ffgen::Limb32>("Fr")
        },
        fft("Fp"),
        if limb64 {
            ffgen::field::<Fq, ffgen::Limb64>("Fq")
        } else {
            ffgen::field::<Fq, ffgen::Limb32>("Fq")
        },
        ec("Fq", "G1"),
        multiexp("G1", "Fp"),
        field2("Fq2", "Fq"),
        ec("Fq2", "G2"),
        multiexp("G2", "Fp"),
    ]
    .join("\n\n")
}
