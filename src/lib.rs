
pub mod domain;
pub mod commit_key;
pub mod opening_key;
pub mod polynomial;
pub mod proof;
pub mod params;
pub mod utils;

pub type G1Point = blstrs::G1Affine;
pub type G2Point = blstrs::G2Affine;
pub type Scalar = blstrs::Scalar;
pub type KZGCommitment = G1Point;
pub type G1Projective = blstrs::G1Projective;

// The number of bytes needed to represent a scalar
pub const SCALAR_SERIALIZED_SIZE: usize = 32;
// The number of bytes needed to represent a compressed G1 point
pub const G1_POINT_SERIALIZED_SIZE: usize = 48;
// The number of bytes needed to represent a compressed G2 point
pub const G2_POINT_SERIALIZED_SIZE: usize = 96;
