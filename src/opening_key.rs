
use pairing_lib::{group::Group, MillerLoopResult, MultiMillerLoop};
use blstrs::{Bls12, G2Prepared};

/// Opening Key is used to verify opening proofs made about a committed polynomial.
#[derive(Clone, Debug)]
pub struct OpeningKey {
    /// The generator of G1 used in the setup
    pub g1_gen: blstrs::G1Affine,
    /// The generator of G2 used in the setup
    pub g2_gen: blstrs::G2Affine,
    /// \tau times the generator of G2.
    pub tau_g2_gen: blstrs::G2Affine,
    /// The generator of G2, prepared for use in pairings.
    pub prepared_g2: G2Prepared,
    /// \tau times the above generator of G2, prepared for use in pairings.
    pub prepared_beta_g2: G2Prepared,
}

impl OpeningKey {
    pub fn new(g1_gen: blstrs::G1Affine, g2_gen: blstrs::G2Affine, tau_g2_gen: blstrs::G2Affine) -> OpeningKey {
        // Store cached elements for verifying multiple proofs.
        let prepared_g2 = G2Prepared::from(g2_gen);
        let prepared_beta_g2 = G2Prepared::from(tau_g2_gen);
        OpeningKey { g1_gen, g2_gen, tau_g2_gen, prepared_g2, prepared_beta_g2, }
    }

    /// Checks that a polynomial `p` was evaluated at a point `z` and returned the value specified `y`.
    /// ie. y = p(z).
    pub fn verify(
        &self,
        input_point: blstrs::Scalar,
        output_point: blstrs::Scalar,
        poly_comm: blstrs::G1Affine,
        witness_comm: blstrs::G1Affine,
    ) -> bool {
        // TODO: .into is doing an inversion. Check if batch normalization saves anything here
        // codepath : G1Projective::batch_normalize(p, q)
        let inner_a: blstrs::G1Affine = (poly_comm - (self.g1_gen * output_point)).into();
        let inner_b: blstrs::G2Affine = (self.tau_g2_gen - (self.g2_gen * input_point)).into();
        let prepared_inner_b = G2Prepared::from(-inner_b);

        let terms = [(&inner_a, &self.prepared_g2), (&witness_comm, &prepared_inner_b)];
        let pairing = Bls12::multi_miller_loop(&terms).final_exponentiation();

        pairing.is_identity().into()
    }
}
