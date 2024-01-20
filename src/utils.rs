use crate::{domain::Domain, polynomial::Polynomial};

use ff::Field;

/// Computes the quotient polynomial for a kzg proof
///
/// The state being proved is p(z) = y
/// Where:
/// - `z` is the point being passed as input
pub fn compute(
    poly: &Polynomial,
    input_point: blstrs::Scalar,
    output_point: blstrs::Scalar,
    domain: &Domain,
) -> Polynomial {
    match domain.find(&input_point) {
        Some(index_in_domain) => {
            compute_quotient_in_domain(poly, index_in_domain, output_point, domain)
        }
        None => compute_quotient_outside_domain(poly, input_point, output_point, domain),
    }
}

fn compute_quotient_in_domain(
    poly: &Polynomial,
    index_in_domain: usize,
    output_point: blstrs::Scalar,
    domain: &Domain,
) -> Polynomial {
    let polynomial_shifted: Vec<_> = poly
        .evaluations
        .iter()
        .map(|evaluation| evaluation - output_point)
        .collect();

    let input_point = domain[index_in_domain];
    let mut denominator_poly: Vec<_> = domain
        .roots()
        .iter()
        .map(|root| root - input_point)
        .collect();
    denominator_poly[index_in_domain] = blstrs::Scalar::one();
    serial_batch_inversion(&mut denominator_poly);

    let mut quotient_poly = vec![blstrs::Scalar::zero(); domain.size()];
    for i in 0..domain.size() {
        if i == index_in_domain {
            quotient_poly[i] =
                compute_quotient_eval_within_domain(poly, index_in_domain, output_point, domain)
        } else {
            quotient_poly[i] = polynomial_shifted[i] * denominator_poly[i]
        }
    }

    Polynomial::new(quotient_poly)
}

fn compute_quotient_eval_within_domain(
    poly: &Polynomial,
    index_in_domain: usize,
    output_point: blstrs::Scalar,
    domain: &Domain,
) -> blstrs::Scalar {
    // TODO Assumes that index_in_domain is in the domain
    // TODO: should we use a special Index struct/enum to encode this?
    let input_point = domain[index_in_domain];

    // TODO: optimize with batch_inverse
    let mut result = blstrs::Scalar::zero();
    for (index, root) in domain.roots().iter().enumerate() {
        if index == index_in_domain {
            continue;
        }

        let f_i = poly[index] - output_point;
        let numerator = f_i * root;
        let denominator = input_point * (input_point - root);
        result += numerator * denominator.invert().unwrap()
    }

    todo!()
}
fn compute_quotient_outside_domain(
    poly: &Polynomial,
    input_point: blstrs::Scalar,
    output_point: blstrs::Scalar,
    domain: &Domain,
) -> Polynomial {
    // Compute the denominator and store it in the quotient vector, to avoid re-allocation
    let mut quotient: Vec<_> = domain
        .roots()
        .iter()
        .map(|domain_element| *domain_element - input_point)
        .collect();
    // This should not panic, since we assume `input_point` is not in the domain
    serial_batch_inversion(&mut quotient);

    // Compute the numerator polynomial and multiply it by the quotient which holds the
    // denominator
    quotient
        .iter_mut()
        .zip(&poly.evaluations)
        .for_each(|(quotient_i, eval_i)| *quotient_i = (*eval_i - output_point) * *quotient_i);

    // Simple way to do this
    // let domain_size = domain.len();
    // let mut quotient = vec![Fr::zero(); domain_size];
    // for i in 0..domain_size {
    // let denominator = inverse(domain[i] - point);
    //     quotient[i] = (poly.evaluations[i] - output) * denominator
    // }
    // quotient

    Polynomial::new(quotient)
}


use std::ops::MulAssign;

/// Given a vector of field elements {v_i}, compute the vector {coeff * v_i^(-1)}
/// This method is explicitly single core.
pub fn serial_batch_inversion(v: &mut [blstrs::Scalar]) {
    
    // Montgomery’s Trick and Fast Implementation of Masked AES
    // Genelle, Prouff and Quisquater
    // Section 3.2
    // but with an optimization to multiply every element in the returned vector by coeff

    // First pass: compute [a, ab, abc, ...]
    let mut prod: Vec<_> = Vec::with_capacity(v.len());
    let mut tmp = blstrs::Scalar::one();
    for f in v.iter().filter(|f| !f.is_zero_vartime()) {
        tmp.mul_assign(f);
        prod.push(tmp);
    }

    assert_eq!(prod.len(), v.len(), "inversion by zero is not allowed");

    // Invert `tmp`.
    tmp = tmp.invert().unwrap(); // Guaranteed to be nonzero.

    // Second pass: iterate backwards to compute inverses
    for (f, s) in v
        .iter_mut()
        // Backwards
        .rev()
        // Ignore normalized elements
        .filter(|f| !f.is_zero_vartime())
        // Backwards, skip last element, fill in one for last term.
        .zip(prod.into_iter().rev().skip(1).chain(Some(blstrs::Scalar::one())))
    {
        // tmp := tmp * f; f := tmp * s = 1/f
        let new_tmp = tmp * *f;
        *f = tmp * &s;
        tmp = new_tmp;
    }
}
