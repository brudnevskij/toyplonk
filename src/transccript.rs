use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use sha2::{Digest, Sha256};

pub fn hash_to_field<F: PrimeField>(label: &str, inputs: &[impl CanonicalSerialize]) -> F {
    let mut hasher = Sha256::new();
    hasher.update(label.as_bytes());

    for input in inputs {
        let mut buf = Vec::new();
        input.serialize_compressed(&mut buf).unwrap();
        hasher.update(&buf);
    }

    let hash_bytes = hasher.finalize();

    F::from_be_bytes_mod_order(&hash_bytes)
}
