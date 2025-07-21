use ark_ec::pairing::Pairing;
use crate::prover::Proof;
use crate::transccript::hash_to_field;

struct Challenges<E: Pairing> {
    alpha: E::ScalarField,
    beta: E::ScalarField,
    gamma: E::ScalarField,
    zeta: E::ScalarField,
    v: E::ScalarField,
    u: E::ScalarField,
}

impl <E: Pairing> Challenges<E> {
    pub fn new(proof: &Proof<E>)->Self{
        let mut commitment_buffer = vec![proof.a, proof.b, proof.c];
        let beta = hash_to_field("beta", &commitment_buffer);
        let gamma = hash_to_field("gamma", &commitment_buffer);
        
        commitment_buffer.push(proof.z);
        let alpha = hash_to_field("alpha", &commitment_buffer);

        commitment_buffer.push(proof.t_lo);
        commitment_buffer.push(proof.t_mid);
        commitment_buffer.push(proof.t_hi);
        let zeta = hash_to_field("zeta", &commitment_buffer);
        let v = hash_to_field("v", &commitment_buffer);
        
        commitment_buffer.push(proof.w_zeta);
        commitment_buffer.push(proof.w_zeta_omega);
        let u = hash_to_field("u", &commitment_buffer);
        
        Challenges{
            alpha,
            beta,
            gamma,
            zeta,
            v,
            u,
        }
    }
}

pub fn verify_kzg_proof<E: Pairing>(proof: &Proof<E>) -> bool{

    true
}