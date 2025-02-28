use crate::{commit_key::*, opening_key::OpeningKey, domain::Domain};

// This is the SRS in lagrange form.
//
// The lagrange form is used to avoid the need to do an inverse fft to commit to polynomials.
pub struct PublicParameters {
    pub commit_key: CommitKeyLagrange,
    pub opening_key: OpeningKey,
}

impl PublicParameters {
    pub fn from_secret_insecure(tau: u64, domain: &Domain) -> Self {
        use ff::Field;
        use group::prime::PrimeCurveAffine;

        let tau_fr = blstrs::Scalar::from(tau);
        let g1_gen = blstrs::G1Affine::generator();
        let g2_gen = blstrs::G2Affine::generator();
        let tau_g2_gen = (g2_gen * tau_fr).into();

        let powers_of_tau_g1: Vec<blstrs::G1Affine> = (0..domain.size())
            .map(|index| {
                let secret_exp = tau_fr.pow_vartime(&[index as u64]);
                (blstrs::G1Affine::generator() * secret_exp).into()
            })
            .collect();

        let commit_key = CommitKey::new(powers_of_tau_g1).into_lagrange(&domain);
        let opening_key = OpeningKey::new(g1_gen, g2_gen, tau_g2_gen);
        PublicParameters { commit_key, opening_key }
    }
}
