use crate::Scalar;
use ff::Field;

// Batch inversion of multiple elements
// This method will panic if one of the elements is zero
pub fn batch_inverse(elements: &mut [Scalar]) {
    batch_inversion(elements)
}

fn batch_inversion(v: &mut [Scalar]) {
    serial_batch_inversion(v);
}

/// Given a vector of field elements {v_i}, compute the vector {coeff * v_i^(-1)}
/// This method is explicitly single core.
fn serial_batch_inversion(v: &mut [Scalar]) {
    use std::ops::MulAssign;

    // Montgomery’s Trick and Fast Implementation of Masked AES
    // Genelle, Prouff and Quisquater
    // Section 3.2
    // but with an optimization to multiply every element in the returned vector by coeff

    // First pass: compute [a, ab, abc, ...]
    let mut prod = Vec::with_capacity(v.len());
    let mut tmp = Scalar::one();
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
        .zip(prod.into_iter().rev().skip(1).chain(Some(Scalar::one())))
    {
        // tmp := tmp * f; f := tmp * s = 1/f
        let new_tmp = tmp * *f;
        *f = tmp * &s;
        tmp = new_tmp;
    }
}
