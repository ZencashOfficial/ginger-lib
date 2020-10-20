#[macro_use]
extern crate criterion;

use criterion::Criterion;
//use criterion::{Criterion, BatchSize};

use algebra::{
    fields::bn_382::Fr,
    curves::bn_382::{G1Projective, G1Affine},
    msm::VariableBaseMSM,
    BigInteger384, PrimeField, UniformRand, ProjectiveCurve,
    ToBytes, FromBytes,
};

use rand_xorshift::XorShiftRng;
use rand::SeedableRng;

use std::fs::File;

//TODO: Maybe "macroize" and use the same bench framework as the other tests in the directory ?

// ***************************************************************************************
// FAST METHOD
// ***************************************************************************************

// fn variable_msm_affine_fast_4(c: &mut Criterion) {
//
//     const PARAM_C: usize = 4;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=4", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_5(c: &mut Criterion) {
//
//     const PARAM_C: usize = 5;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=5", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_6(c: &mut Criterion) {
//
//     const PARAM_C: usize = 6;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=6", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_7(c: &mut Criterion) {
//
//     const PARAM_C: usize = 7;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=7", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_8(c: &mut Criterion) {
//
//     const PARAM_C: usize = 8;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=8", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_9(c: &mut Criterion) {
//
//     const PARAM_C: usize = 9;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=9", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_10(c: &mut Criterion) {
//
//     const PARAM_C: usize = 10;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=10", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_11(c: &mut Criterion) {
//
//     const PARAM_C: usize = 11;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=11", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_12(c: &mut Criterion) {
//
//     const PARAM_C: usize = 12;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=12", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_13(c: &mut Criterion) {
//
//     const PARAM_C: usize = 13;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=13", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_14(c: &mut Criterion) {
//
//     const PARAM_C: usize = 14;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=14", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_15(c: &mut Criterion) {
//
//     const PARAM_C: usize = 15;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=15", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_16(c: &mut Criterion) {
//
//     const PARAM_C: usize = 16;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=16", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_17(c: &mut Criterion) {
//
//     const PARAM_C: usize = 17;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=17", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_18(c: &mut Criterion) {
//
//     const PARAM_C: usize = 18;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=18", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_19(c: &mut Criterion) {
//
//     const PARAM_C: usize = 19;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=19", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_20(c: &mut Criterion) {
//
//     const PARAM_C: usize = 20;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=20", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_21(c: &mut Criterion) {
//
//     const PARAM_C: usize = 21;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=21", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_22(c: &mut Criterion) {
//
//     const PARAM_C: usize = 22;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=22", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_fast_23(c: &mut Criterion) {
//
//     const PARAM_C: usize = 23;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine_fast c=23", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }

// ***************************************************************************************
// AFFINE
// ***************************************************************************************

// fn variable_msm_affine_4(c: &mut Criterion) {
//
//     const PARAM_C: usize = 4;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=4", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_5(c: &mut Criterion) {
//
//     const PARAM_C: usize = 5;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=5", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_6(c: &mut Criterion) {
//
//     const PARAM_C: usize = 6;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=6", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_7(c: &mut Criterion) {
//
//     const PARAM_C: usize = 7;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=7", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_8(c: &mut Criterion) {
//
//     const PARAM_C: usize = 8;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=8", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_9(c: &mut Criterion) {
//
//     const PARAM_C: usize = 9;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=9", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_10(c: &mut Criterion) {
//
//     const PARAM_C: usize = 10;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=10", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_11(c: &mut Criterion) {
//
//     const PARAM_C: usize = 11;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=11", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_12(c: &mut Criterion) {
//
//     const PARAM_C: usize = 12;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=12", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_13(c: &mut Criterion) {
//
//     const PARAM_C: usize = 13;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=13", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_14(c: &mut Criterion) {
//
//     const PARAM_C: usize = 14;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=14", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_15(c: &mut Criterion) {
//
//     const PARAM_C: usize = 15;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=15", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_16(c: &mut Criterion) {
//
//     const PARAM_C: usize = 16;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=16", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_17(c: &mut Criterion) {
//
//     const PARAM_C: usize = 17;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=17", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_18(c: &mut Criterion) {
//
//     const PARAM_C: usize = 18;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=18", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_19(c: &mut Criterion) {
//
//     const PARAM_C: usize = 19;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=19", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_20(c: &mut Criterion) {
//
//     const PARAM_C: usize = 20;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=20", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_21(c: &mut Criterion) {
//
//     const PARAM_C: usize = 21;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=21", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_22(c: &mut Criterion) {
//
//     const PARAM_C: usize = 22;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=22", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }
// fn variable_msm_affine_23(c: &mut Criterion) {
//
//     const PARAM_C: usize = 23;
//
//     let (v, g) = load_data();
//
//     c.bench_function("Variable MSM with affine c=23", move |b| {
//         b.iter(|| {
//             VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(),PARAM_C);
//         }
//         )
//     });
// }

// ***************************************************************************************
// AFFINE SD
// ***************************************************************************************

fn variable_msm_affine_sd_4(c: &mut Criterion) {

    const PARAM_C: usize = 4;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=4", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_5(c: &mut Criterion) {

    const PARAM_C: usize = 5;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=5", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_6(c: &mut Criterion) {

    const PARAM_C: usize = 6;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=6", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_7(c: &mut Criterion) {

    const PARAM_C: usize = 7;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=7", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_8(c: &mut Criterion) {

    const PARAM_C: usize = 8;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=8", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_9(c: &mut Criterion) {

    const PARAM_C: usize = 9;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=9", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_10(c: &mut Criterion) {

    const PARAM_C: usize = 10;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=10", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_11(c: &mut Criterion) {

    const PARAM_C: usize = 11;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=11", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_12(c: &mut Criterion) {

    const PARAM_C: usize = 12;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=12", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_13(c: &mut Criterion) {

    const PARAM_C: usize = 13;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=13", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_14(c: &mut Criterion) {

    const PARAM_C: usize = 14;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=14", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_15(c: &mut Criterion) {

    const PARAM_C: usize = 15;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=15", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_16(c: &mut Criterion) {

    const PARAM_C: usize = 16;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=16", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_17(c: &mut Criterion) {

    const PARAM_C: usize = 17;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=17", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_18(c: &mut Criterion) {

    const PARAM_C: usize = 18;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=18", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_19(c: &mut Criterion) {

    const PARAM_C: usize = 19;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=19", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_20(c: &mut Criterion) {

    const PARAM_C: usize = 20;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=20", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_21(c: &mut Criterion) {

    const PARAM_C: usize = 21;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=21", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_22(c: &mut Criterion) {

    const PARAM_C: usize = 22;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=22", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}
fn variable_msm_affine_sd_23(c: &mut Criterion) {

    const PARAM_C: usize = 23;

    let (v, g) = load_data();

    c.bench_function("Variable MSM with affine_sd c=23", move |b| {
        b.iter(|| {
            VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(),PARAM_C);
        }
        )
    });
}

const SAMPLES: usize = 1<<23;

fn load_data() -> (Vec<BigInteger384>,Vec<G1Affine>) {

    let mut fs = File::open("./scalars_bases").unwrap();
    let mut v = Vec::with_capacity(SAMPLES);
    let mut g = Vec::with_capacity(SAMPLES);

    for _i in 0..SAMPLES {
        let elem = BigInteger384::read(&mut fs).unwrap();
        v.push(elem);
        let elem = G1Affine::read(&mut fs).unwrap();
        g.push(elem);
    }
    (v,g)
}

#[allow(dead_code)]
fn generate_data() {

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let mut scalar_fs = File::create("./scalars_bases").unwrap();

    for _i in 0..SAMPLES {
        let elem = Fr::rand(&mut rng).into_repr();
        elem.write(&mut scalar_fs).unwrap();
        let elem = G1Projective::rand(&mut rng).into_affine();
        elem.write(&mut scalar_fs).unwrap();
    }
}

criterion_group! {
    name = variable_msm_eval;
    config = Criterion::default().sample_size(10);
    //targets = variable_msm_affine_fast_4,variable_msm_affine_fast_5,variable_msm_affine_fast_6,variable_msm_affine_fast_7,variable_msm_affine_fast_8,variable_msm_affine_fast_9,variable_msm_affine_fast_10,variable_msm_affine_fast_11,variable_msm_affine_fast_12,variable_msm_affine_fast_13,variable_msm_affine_fast_14,variable_msm_affine_fast_15,variable_msm_affine_fast_16,variable_msm_affine_fast_17,variable_msm_affine_fast_18,variable_msm_affine_fast_19,variable_msm_affine_fast_20,variable_msm_affine_fast_21,variable_msm_affine_fast_22,variable_msm_affine_fast_23
    //targets = variable_msm_affine_4,variable_msm_affine_5,variable_msm_affine_6,variable_msm_affine_7,variable_msm_affine_8,variable_msm_affine_9,variable_msm_affine_10,variable_msm_affine_11,variable_msm_affine_12,variable_msm_affine_13,variable_msm_affine_14,variable_msm_affine_15,variable_msm_affine_16,variable_msm_affine_17,variable_msm_affine_18,variable_msm_affine_19,variable_msm_affine_20,variable_msm_affine_21,variable_msm_affine_22,variable_msm_affine_23
    targets = variable_msm_affine_sd_4,variable_msm_affine_sd_5,variable_msm_affine_sd_6,variable_msm_affine_sd_7,variable_msm_affine_sd_8,variable_msm_affine_sd_9,variable_msm_affine_sd_10,variable_msm_affine_sd_11,variable_msm_affine_sd_12,variable_msm_affine_sd_13,variable_msm_affine_sd_14,variable_msm_affine_sd_15,variable_msm_affine_sd_16,variable_msm_affine_sd_17,variable_msm_affine_sd_18,variable_msm_affine_sd_19,variable_msm_affine_sd_20,variable_msm_affine_sd_21,variable_msm_affine_sd_22,variable_msm_affine_sd_23

}

criterion_main!(
    variable_msm_eval
);