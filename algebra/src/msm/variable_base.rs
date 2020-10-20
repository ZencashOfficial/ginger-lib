use crate::{AffineCurve, BigInteger, FpParameters, Field, PrimeField, ProjectiveCurve};
use rayon::prelude::*;

pub struct VariableBaseMSM;

impl VariableBaseMSM {

    // Function that recodes the scalars into SD numbers
    // The output is a vector
    pub fn recode_sd<G: AffineCurve>(
        scalar: &<G::ScalarField as PrimeField>::BigInt,
        c: usize
    ) -> Vec<i64> {

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;

        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        let mut vec_coeff = Vec::new();

        window_starts.iter().rev().for_each(|x| {
            let mut scal = (*scalar).clone();
            scal.divn(*x as u32);
            // We mod the remaining bits by the window size.
            let a = scal.as_ref()[0] % (1 << c);
            vec_coeff.push(a as i64);
        });

        for idx in (0..vec_coeff.len()).rev() {
            if vec_coeff[idx] >= (1 << (c-1)) {
                vec_coeff[idx] -= 1 << c;
                if idx!= 0 {
                    vec_coeff[idx-1] += 1;
                }
            }
        }

        // last element is the least significant digit
        return vec_coeff;

    }

    // Recodes in SD the portions of the scalar divided in c bits
    // Takes as input the carry_in and generates the vector of SD digits corresponding to
    // the processing chunk and the carry_out
    // For example: for 4 cores, it takes the carry_in and
    // generates carry_out, a_(i+3), a_(i+2), a_(i+1), a_i
    #[inline]
    pub fn recode_sd_chunk<G:AffineCurve>(
        scalar: &<G::ScalarField as PrimeField>::BigInt,
        c: usize,
        chunk_pos: usize,
        num_cores: usize,
        vec_coeff: &mut Vec<i64>,
        carry_in: &mut i64
    )   {

        // starting bit position of the chunk
        let start_w_bit = chunk_pos * (c * num_cores);

        for i in 0..num_cores {
            let windows_starts = start_w_bit + i * c;
            let mut scal = (*scalar).clone();
            scal.divn(windows_starts as u32);
            let mut a = (scal.as_ref()[0] % (1 << c)) as i64 + *carry_in;
            if a >= (1 << (c-1)) {
                a -= 1 << c;
                *carry_in = 1;
            } else {
                *carry_in = 0;
            }
            (*vec_coeff)[i] = a;
        }
    }

    pub fn multi_scalar_mul_affine_sd<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt]
    ) -> G::Projective {

        //In the case of SD recoding, we can use one more bit for c for the same amount of memory usage
        let c = if scalars.len() < 32 {
            3 + 1
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize + 1
        };

        let cc = 1 << c;

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;

        let zero = G::zero().into_projective();

        // The number of windows of the scalars in groups of c bits
        let num_w = (num_bits as f64 / c as f64).ceil() as usize;
        // The number of cores present
        let cpus = rayon::current_num_threads();
        // The number of chunks of cpus digits
        let num_chunks = (num_w as f64/cpus as f64).floor() as usize;
        // The remaining digits to process sequentially
        let remaining_digits = num_w - (num_chunks * cpus);
        // println!("num_w = {}", num_w);
        // println!("cpus = {}", cpus);
        // println!("c = {}", c);
        // println!("remaining_digits = {}", remaining_digits);
        // println!("num_chunks = {}", num_chunks);

        let mut window_sums =  Vec::new();
        let mut small_window_sums:Vec<_>;

        let mut carry_vec = vec![0; scalars.len()];
        let mut vec_coeff = vec![vec![0; cpus]; scalars.len()];

        for i in 0..num_chunks {
            // index of the digit within the small window of cpus number of digits
            let idx: Vec<_> = (0..cpus).rev().collect();

            vec_coeff.par_iter_mut().zip(carry_vec.par_iter_mut()).enumerate().for_each( |(l, (v1, c1))| {
                Self::recode_sd_chunk::<G>(&scalars[l], c, i, cpus, v1, c1);
            });

            small_window_sums = idx
                .into_par_iter()
                .map(|w_idx| {
                    // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                    let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc / 2];
                    for i in 0..scalars.len() {
                        if !scalars[i].is_zero() {

                            let scalar = vec_coeff[i][w_idx];

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 && bases[i].is_zero() == false {
                                if scalar < 0 {
                                    buckets[(-scalar - 1) as usize].push(-(bases[i]));
                                } else {
                                    buckets[(scalar - 1) as usize].push(bases[i]);
                                }
                            }
                        }
                    }
                    G::add_points(&mut buckets);
                    let mut res = zero.clone();

                    let mut running_sum = zero.clone();
                    for b in buckets[0..cc / 2].iter_mut().rev() {
                        if b.len() != 0 && b[0].is_zero() == false {
                            running_sum.add_assign_mixed(&b[0])
                        }
                        res += &running_sum;
                    }

                    res
                })
                .collect();

            small_window_sums.iter().rev().for_each(|x| {
                window_sums.push(x.clone());
            });
        }

        if remaining_digits != 0 {

            let idx:Vec<_> = (0..remaining_digits).rev().collect();

            vec_coeff.par_iter_mut().zip(carry_vec.par_iter_mut()).enumerate().for_each( |(l, (v1, c1))| {
                Self::recode_sd_chunk::<G>(&scalars[l], c, num_chunks, cpus, v1, c1);
            });

            // Each window is of size `c`.
            // We divide up the bits 0..num_bits into windows of size `c`, and
            // in parallel process each such window.
            small_window_sums = idx
                .into_par_iter()
                .map(|w_idx| {
                    // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                    let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc / 2];
                    for i in 0..scalars.len() {
                        if !scalars[i].is_zero() {

                            let scalar = vec_coeff[i][w_idx];

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 && bases[i].is_zero() == false {
                                if scalar < 0 {
                                    buckets[(-scalar - 1) as usize].push(-(bases[i]));
                                } else {
                                    buckets[(scalar - 1) as usize].push(bases[i]);
                                }
                            }
                        }
                    }
                    G::add_points(&mut buckets);
                    let mut res = zero.clone();

                    let mut running_sum = zero.clone();
                    for b in buckets[0..cc / 2].iter_mut().rev() {
                        if b.len() != 0 && b[0].is_zero() == false {
                            running_sum.add_assign_mixed(&b[0])
                        }
                        res += &running_sum;
                    }
                    res
                })
                .collect();

            small_window_sums.iter().rev().for_each(|x| {
                window_sums.push(x.clone());
            });
        }

        // We store the sum for the lowest window.
        let lowest = window_sums.first().unwrap();

        // We're traversing windows from high to low.
        window_sums[1..].iter().rev().fold(zero, |mut total, &sum_i| {
            total += &sum_i;
            for _ in 0..c {
                total.double_in_place();
            }
            total
        }) + lowest
    }

    pub fn multi_scalar_mul_affine<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt]
    ) -> G::Projective {
        let c = if scalars.len() < 32 {
            3
        } else {
            (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
        };
        let cc = 1 << c;

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = G::ScalarField::one().into_repr();

        let zero = G::zero().into_projective();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = window_starts
            .into_par_iter()
            .map(|w_start| {
                // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                let mut buckets = vec![Vec::with_capacity(bases.len()/cc*2); cc];
                scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero()).for_each(|(&scalar, base)|  {
                    if scalar == fr_one {
                        // We only process unit scalars once in the first window.
                        if w_start == 0 && base.is_zero() == false {
                            buckets[cc-1].push(*base);
                        }
                    } else {
                        let mut scalar = scalar;

                        // We right-shift by w_start, thus getting rid of the
                        // lower bits.
                        scalar.divn(w_start as u32);

                        // We mod the remaining bits by the window size.
                        let scalar = scalar.as_ref()[0] % (1 << c);

                        // If the scalar is non-zero, we update the corresponding
                        // bucket.
                        // (Recall that `buckets` doesn't have a zero bucket.)
                        if scalar != 0 && base.is_zero() == false {
                            buckets[(scalar - 1) as usize].push(*base);
                        }
                    }
                });
                G::add_points(&mut buckets);
                let mut res = if buckets[cc-1].len() == 0 {zero} else {buckets[cc-1][0].into_projective()};

                let mut running_sum = zero;
                for b in buckets[0..cc-1].iter_mut().rev() {
                    if b.len() != 0 && b[0].is_zero() == false {
                        running_sum.add_assign_mixed(&b[0])
                    }
                    res += &running_sum;
                }
                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = window_sums.first().unwrap();

        // We're traversing windows from high to low.
        window_sums[1..].iter().rev().fold(zero, |mut total, &sum_i| {
            total += &sum_i;
            for _ in 0..c {
                total.double_in_place();
            }
            total
        }) + lowest
    }

    fn msm_inner<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective
        where
            G::Projective: ProjectiveCurve<Affine = G>,
    {
        let c = if scalars.len() < 32 {
            3
        } else {
            super::ln_without_floats(scalars.len()) + 2
        };

        let num_bits = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = G::ScalarField::one().into_repr();

        let zero = G::Projective::zero();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        let window_starts_iter = window_starts.into_par_iter();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = window_starts_iter
            .map(|w_start| {
                let mut res = zero;
                // We don't need the "zero" bucket, so we only have 2^c - 1 buckets
                let mut buckets = vec![zero; (1 << c) - 1];
                scalars
                    .iter()
                    .zip(bases)
                    .filter(|(s, _)| !s.is_zero())
                    .for_each(|(&scalar, base)| {
                        if scalar == fr_one {
                            // We only process unit scalars once in the first window.
                            if w_start == 0 {
                                res.add_assign_mixed(base);
                            }
                        } else {
                            let mut scalar = scalar;

                            // We right-shift by w_start, thus getting rid of the
                            // lower bits.
                            scalar.divn(w_start as u32);

                            // We mod the remaining bits by the window size.
                            let scalar = scalar.as_ref()[0] % (1 << c);

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 {
                                buckets[(scalar - 1) as usize].add_assign_mixed(base);
                            }
                        }
                    });
                let buckets = G::Projective::batch_normalization_into_affine(&buckets);

                let mut running_sum = G::Projective::zero();
                for b in buckets.into_iter().rev() {
                    running_sum.add_assign_mixed(&b);
                    res += &running_sum;
                }

                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = *window_sums.first().unwrap();

        // We're traversing windows from high to low.
        window_sums[1..]
            .iter()
            .rev()
            .fold(zero, |mut total, &sum_i| {
                total += &sum_i;
                for _ in 0..c {
                    total.double_in_place();
                }
                total
            }) + &lowest
    }

    pub fn multi_scalar_mul<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        Self::msm_inner(bases, scalars)
    }

    // *********************************************************************************************
    // Parametrized by the window size c
    // *********************************************************************************************

    pub fn multi_scalar_mul_affine_sd_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c:usize
    ) -> G::Projective {

        // In the case of SD recoding, we can use one more bit for c for the same amount of memory usage
        // let c = if scalars.len() < 32 {
        //     3 + 1
        // } else {
        //     (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize + 1
        // };

        let cc = 1 << c;

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;

        let zero = G::zero().into_projective();

        // The number of windows of the scalars in groups of c bits
        let num_w = (num_bits as f64 / c as f64).ceil() as usize;
        // The number of cores present
        let cpus = rayon::current_num_threads();
        // The number of chunks of cpus digits
        let num_chunks = (num_w as f64/cpus as f64).floor() as usize;
        // The remaining digits to process sequentially
        let remaining_digits = num_w - (num_chunks * cpus);
        // println!("num_w = {}", num_w);
        // println!("cpus = {}", cpus);
        // println!("c = {}", c);
        // println!("remaining_digits = {}", remaining_digits);
        // println!("num_chunks = {}", num_chunks);

        let mut window_sums =  Vec::new();
        let mut small_window_sums:Vec<_>;

        let mut carry_vec = vec![0; scalars.len()];
        let mut vec_coeff = vec![vec![0; cpus]; scalars.len()];

        for i in 0..num_chunks {
            // index of the digit within the small window of cpus number of digits
            let idx: Vec<_> = (0..cpus).rev().collect();

            vec_coeff.par_iter_mut().zip(carry_vec.par_iter_mut()).enumerate().for_each( |(l, (v1, c1))| {
                Self::recode_sd_chunk::<G>(&scalars[l], c, i, cpus, v1, c1);
            });

            small_window_sums = idx
                .into_par_iter()
                .map(|w_idx| {
                    // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                    let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc / 2];
                    for i in 0..scalars.len() {
                        if !scalars[i].is_zero() {

                            let scalar = vec_coeff[i][w_idx];

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 && bases[i].is_zero() == false {
                                if scalar < 0 {
                                    buckets[(-scalar - 1) as usize].push(-(bases[i]));
                                } else {
                                    buckets[(scalar - 1) as usize].push(bases[i]);
                                }
                            }
                        }
                    }
                    G::add_points(&mut buckets);
                    let mut res = zero.clone();

                    let mut running_sum = zero.clone();
                    for b in buckets[0..cc / 2].iter_mut().rev() {
                        if b.len() != 0 && b[0].is_zero() == false {
                            running_sum.add_assign_mixed(&b[0])
                        }
                        res += &running_sum;
                    }

                    res
                })
                .collect();

            small_window_sums.iter().rev().for_each(|x| {
                window_sums.push(x.clone());
            });
        }

        if remaining_digits != 0 {

            let idx:Vec<_> = (0..remaining_digits).rev().collect();

            vec_coeff.par_iter_mut().zip(carry_vec.par_iter_mut()).enumerate().for_each( |(l, (v1, c1))| {
                Self::recode_sd_chunk::<G>(&scalars[l], c, num_chunks, cpus, v1, c1);
            });

            // Each window is of size `c`.
            // We divide up the bits 0..num_bits into windows of size `c`, and
            // in parallel process each such window.
            small_window_sums = idx
                .into_par_iter()
                .map(|w_idx| {
                    // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                    let mut buckets = vec![Vec::with_capacity(bases.len() / cc * 2); cc / 2];
                    for i in 0..scalars.len() {
                        if !scalars[i].is_zero() {

                            let scalar = vec_coeff[i][w_idx];

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 && bases[i].is_zero() == false {
                                if scalar < 0 {
                                    buckets[(-scalar - 1) as usize].push(-(bases[i]));
                                } else {
                                    buckets[(scalar - 1) as usize].push(bases[i]);
                                }
                            }
                        }
                    }
                    G::add_points(&mut buckets);
                    let mut res = zero.clone();

                    let mut running_sum = zero.clone();
                    for b in buckets[0..cc / 2].iter_mut().rev() {
                        if b.len() != 0 && b[0].is_zero() == false {
                            running_sum.add_assign_mixed(&b[0])
                        }
                        res += &running_sum;
                    }
                    res
                })
                .collect();

            small_window_sums.iter().rev().for_each(|x| {
                window_sums.push(x.clone());
            });
        }

        // We store the sum for the lowest window.
        let lowest = window_sums.first().unwrap();

        // We're traversing windows from high to low.
        window_sums[1..].iter().rev().fold(zero, |mut total, &sum_i| {
            total += &sum_i;
            for _ in 0..c {
                total.double_in_place();
            }
            total
        }) + lowest
    }

    pub fn multi_scalar_mul_affine_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c: usize
    ) -> G::Projective {
        // let c = if scalars.len() < 32 {
        //     3
        // } else {
        //     (2.0 / 3.0 * (f64::from(scalars.len() as u32)).log2() - 2.0).ceil() as usize
        // };
        let cc = 1 << c;

        let num_bits =
            <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = G::ScalarField::one().into_repr();

        let zero = G::zero().into_projective();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = window_starts
            .into_par_iter()
            .map(|w_start| {
                // We don't need the "zero" bucket, we use 2^c-1 bucket for units
                let mut buckets = vec![Vec::with_capacity(bases.len()/cc*2); cc];
                scalars.iter().zip(bases).filter(|(s, _)| !s.is_zero()).for_each(|(&scalar, base)|  {
                    if scalar == fr_one {
                        // We only process unit scalars once in the first window.
                        if w_start == 0 && base.is_zero() == false {
                            buckets[cc-1].push(*base);
                        }
                    } else {
                        let mut scalar = scalar;

                        // We right-shift by w_start, thus getting rid of the
                        // lower bits.
                        scalar.divn(w_start as u32);

                        // We mod the remaining bits by the window size.
                        let scalar = scalar.as_ref()[0] % (1 << c);

                        // If the scalar is non-zero, we update the corresponding
                        // bucket.
                        // (Recall that `buckets` doesn't have a zero bucket.)
                        if scalar != 0 && base.is_zero() == false {
                            buckets[(scalar - 1) as usize].push(*base);
                        }
                    }
                });
                G::add_points(&mut buckets);
                let mut res = if buckets[cc-1].len() == 0 {zero} else {buckets[cc-1][0].into_projective()};

                let mut running_sum = zero;
                for b in buckets[0..cc-1].iter_mut().rev() {
                    if b.len() != 0 && b[0].is_zero() == false {
                        running_sum.add_assign_mixed(&b[0])
                    }
                    res += &running_sum;
                }
                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = window_sums.first().unwrap();

        // We're traversing windows from high to low.
        window_sums[1..].iter().rev().fold(zero, |mut total, &sum_i| {
            total += &sum_i;
            for _ in 0..c {
                total.double_in_place();
            }
            total
        }) + lowest
    }

    fn msm_inner_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c:usize
    ) -> G::Projective
        where
            G::Projective: ProjectiveCurve<Affine = G>,
    {
        // let c = if scalars.len() < 32 {
        //     3
        // } else {
        //     super::ln_without_floats(scalars.len()) + 2
        // };

        let num_bits = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as usize;
        let fr_one = G::ScalarField::one().into_repr();

        let zero = G::Projective::zero();
        let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

        let window_starts_iter = window_starts.into_par_iter();

        // Each window is of size `c`.
        // We divide up the bits 0..num_bits into windows of size `c`, and
        // in parallel process each such window.
        let window_sums: Vec<_> = window_starts_iter
            .map(|w_start| {
                let mut res = zero;
                // We don't need the "zero" bucket, so we only have 2^c - 1 buckets
                let mut buckets = vec![zero; (1 << c) - 1];
                scalars
                    .iter()
                    .zip(bases)
                    .filter(|(s, _)| !s.is_zero())
                    .for_each(|(&scalar, base)| {
                        if scalar == fr_one {
                            // We only process unit scalars once in the first window.
                            if w_start == 0 {
                                res.add_assign_mixed(base);
                            }
                        } else {
                            let mut scalar = scalar;

                            // We right-shift by w_start, thus getting rid of the
                            // lower bits.
                            scalar.divn(w_start as u32);

                            // We mod the remaining bits by the window size.
                            let scalar = scalar.as_ref()[0] % (1 << c);

                            // If the scalar is non-zero, we update the corresponding
                            // bucket.
                            // (Recall that `buckets` doesn't have a zero bucket.)
                            if scalar != 0 {
                                buckets[(scalar - 1) as usize].add_assign_mixed(base);
                            }
                        }
                    });
                let buckets = G::Projective::batch_normalization_into_affine(&buckets);

                let mut running_sum = G::Projective::zero();
                for b in buckets.into_iter().rev() {
                    running_sum.add_assign_mixed(&b);
                    res += &running_sum;
                }

                res
            })
            .collect();

        // We store the sum for the lowest window.
        let lowest = *window_sums.first().unwrap();

        // We're traversing windows from high to low.
        window_sums[1..]
            .iter()
            .rev()
            .fold(zero, |mut total, &sum_i| {
                total += &sum_i;
                for _ in 0..c {
                    total.double_in_place();
                }
                total
            }) + &lowest
    }

    pub fn multi_scalar_mul_c<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
        c:usize
    ) -> G::Projective {
        Self::msm_inner_c(bases, scalars, c)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::time::Instant;
    use crate::{
        fields::bn_382::Fr,
        curves::bn_382::{
            G1Projective, G1Affine,
        },
        UniformRand,
        FixedBaseMSM,
    };
    use rand::{
        Rng, SeedableRng
    };
    use rand_xorshift::XorShiftRng;

    fn naive_var_base_msm<G: AffineCurve>(
        bases: &[G],
        scalars: &[<G::ScalarField as PrimeField>::BigInt],
    ) -> G::Projective {
        let mut acc = <G::Projective as ProjectiveCurve>::zero();

        for (base, scalar) in bases.iter().zip(scalars.iter()) {
            acc += &base.mul(*scalar);
        }
        acc
    }

    #[test]
    fn test_with_bn_382_c() {

        const SAMPLES: usize = 1 << 23;

        let mut rng = XorShiftRng::seed_from_u64(234872845u64);

        println!("Generating scalars...");
        let v = (0..SAMPLES)
            .map(|_| Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        println!("Generating bases...");
        let g = (0..SAMPLES)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        for c in 4..25 {

            let now_2 = Instant::now();
            let fast = VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(), c);
            let now_3 = Instant::now();

            let now_4 = Instant::now();
            let affine = VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(), c);
            let now_5 = Instant::now();

            let now_6 = Instant::now();
            let affine_sd = VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(), c);
            let now_7 = Instant::now();

            println!("==========================================================================");
            println!("c = {}", c);
            println!("Duration fast      = {:?}", now_3.duration_since(now_2).as_micros());
            println!("Duration affine    = {:?}", now_5.duration_since(now_4).as_micros());
            println!("Duration affine_sd = {:?}", now_7.duration_since(now_6).as_micros());

            assert_eq!(fast, affine);
            assert_eq!(affine, affine_sd);
        }

    }

    #[test]
    fn test_all() {

        const SAMPLES: usize = 1 << 10;

        let mut rng = XorShiftRng::seed_from_u64(234872845u64);

        println!("Generating scalars...");
        let v = (0..SAMPLES)
            .map(|_| Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        println!("Generating bases...");
        let g = (0..SAMPLES)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        println!("Starting naive...");
        let now_0 = Instant::now();
        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let now_1 = Instant::now();

        println!("Starting fast...");
        let now_2 = Instant::now();
        let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());
        let now_3 = Instant::now();

        println!("Starting affine...");
        let now_4 = Instant::now();
        let affine = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());
        let now_5 = Instant::now();

        println!("Starting affine_sd...");
        let now_6 = Instant::now();
        let affine_sd = VariableBaseMSM::multi_scalar_mul_affine_sd(g.as_slice(), v.as_slice());
        let now_7 = Instant::now();

        println!("Duration naive     = {:?}", now_1.duration_since(now_0).as_micros());
        println!("Duration fast      = {:?}", now_3.duration_since(now_2).as_micros());
        println!("Duration affine    = {:?}", now_5.duration_since(now_4).as_micros());
        println!("Duration affine_sd = {:?}", now_7.duration_since(now_6).as_micros());

        assert_eq!(naive, fast);
        assert_eq!(fast, affine);
        assert_eq!(affine, affine_sd);

    }

    #[test]
    fn test_with_unequal_numbers() {
        const SAMPLES: usize = 1 << 10;

        let mut rng = XorShiftRng::seed_from_u64(234872845u64);

        let v = (0..SAMPLES-1)
            .map(|_| Fr::rand(&mut rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..SAMPLES)
            .map(|_| G1Projective::rand(&mut rng).into_affine())
            .collect::<Vec<_>>();

        let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
        let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());
        let affine = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());

        assert_eq!(naive, fast);
        assert_eq!(naive, affine)
    }

    #[test]
    fn batch_addition()
    {
        let mut length = 1000000;
        let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

        let size_in_bits = <Fr as PrimeField>::size_in_bits();
        let window_size = FixedBaseMSM::get_mul_window_size(length);
        let mut v = FixedBaseMSM::multi_scalar_mul::<G1Projective> (
            size_in_bits,
            window_size,
            &FixedBaseMSM::get_window_table
                (
                    size_in_bits,
                    window_size,
                    G1Projective::prime_subgroup_generator()
                ),
            &(0..length).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>(),
        );
        ProjectiveCurve::batch_normalization(&mut v);

        let vector = v.iter().map(|e| e.into_affine()).collect::<Vec<_>>();

        println!();
        println!("{}", "*****BENCHMARKING OPTIMISED AFFINE VS SERIAL JACOBIAN MIXED ADDITION*****");
        loop {
            let mut vectors = (0..1000000/length).map(|_| vector[rng.gen_range(0, length/2)..rng.gen_range(length/2, length)].to_vec()).collect::<Vec<_>>();

            //println!("{}{:?}", "Lengths: ".magenta(), vectors.iter().map(|v| v.len()).collect::<Vec<_>>());
            println!();
            println!("{}{:?}", "Length: ", length);
            println!("{}{:?}", "Buckets: ", 1000000/length);

            let start = Instant::now();

            let sum = vectors.iter().map
            (
                |v|
                    {
                        let mut sum = G1Projective::zero();
                        for point in v.iter()
                            {
                                sum.add_assign_mixed(point);
                            }
                        sum.into_affine()
                    }
            ).collect::<Vec<_>>();

            let serial = start.elapsed();
            println!("     {}{:?}", "serial: ", serial);
            let start = Instant::now();

            AffineCurve::add_points(&mut vectors);

            let batch = start.elapsed();
            println!("     {}{:?}", "batch: ", batch);
            println!("     {}{:?}", "batch/serial: ", batch.as_secs_f32()/serial.as_secs_f32());

            assert_eq!(vectors.iter().map(|v| if v.len() == 0 {G1Affine::zero()} else {v[0]}).collect::<Vec<_>>().iter().eq(sum.iter()), true);

            length = length/10;
            if length == 1 {break}
        }
    }

    #[test]
    fn multiexp() {
        let mut length = 1000000;
        let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

        let size_in_bits = <Fr as PrimeField>::size_in_bits();
        let window_size = FixedBaseMSM::get_mul_window_size(length);
        let mut v = FixedBaseMSM::multi_scalar_mul::<G1Projective> (
            size_in_bits,
            window_size,
            &FixedBaseMSM::get_window_table
                (
                    size_in_bits,
                    window_size,
                    G1Projective::prime_subgroup_generator()
                ),
            &(0..length).map(|_| Fr::rand(rng)).collect::<Vec<Fr>>(),
        );
        ProjectiveCurve::batch_normalization(&mut v);

        let bases = v.iter().map(|e| e.into_affine()).collect::<Vec<_>>();
        let scalars = (0..length).map(|_| Fr::rand(rng).into_repr()).collect::<Vec<_>>();

        println!();
        println!("{}", "*****BENCHMARKING MULTIEXP WITH OPTIMISED AFFINE VS*****");
        println!("{}", "******MULTIEXP WITH SERIAL JACOBIAN MIXED ADDITION******");
        length = 1000000;
        loop {
            let base = bases[0..length].to_vec();
            let scalar = scalars[0..length].to_vec();

            println!("{}{:?}", "Length: ", length);

            let start = Instant::now();

            let s1 = VariableBaseMSM::multi_scalar_mul_affine(&base, &scalar);

            let batch = start.elapsed();
            println!("     {}{:?}", "batch: ", batch);
            let start = Instant::now();

            let s2 = VariableBaseMSM::multi_scalar_mul(&base, &scalar);

            let serial = start.elapsed();
            println!("     {}{:?}", "serial: ", serial);
            println!("     {}{:?}", "batch/serial: ", batch.as_secs_f32()/serial.as_secs_f32());

            assert_eq!(s1, s2);
            if length == 1 {break}
            length = length/10;
        }
    }
}
