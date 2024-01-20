
pub mod domain;
pub mod commit_key;
pub mod opening_key;
pub mod polynomial;
pub mod proof;
pub mod params;
pub mod utils;

// The number of bytes needed to represent a scalar
pub const SCALAR_SERIALIZED_SIZE: usize = 32;
// The number of bytes needed to represent a compressed G1 point
pub const G1_POINT_SERIALIZED_SIZE: usize = 48;
// The number of bytes needed to represent a compressed G2 point
pub const G2_POINT_SERIALIZED_SIZE: usize = 96;
