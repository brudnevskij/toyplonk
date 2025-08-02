use crate::circuit::{Circuit, Witness};
use crate::gate::Gate;
use crate::prover::KZGProver;
use crate::verifier::{VerifierPreprocessedInput, verify_kzg_proof};
use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, Group};
use ark_ff::{FftField, Field};
use ark_std::{UniformRand, test_rng};

mod circuit;
mod fft;
mod gate;
mod permutation;
mod prover;
mod verifier;

fn fr(n: u64) -> Fr {
    Fr::from(n)
}

fn get_srs<E: Pairing>(
    g1: E::G1Affine,
    g2: E::G2Affine,
    degree: usize,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let mut rng = test_rng();
    let tau = E::ScalarField::rand(&mut rng);

    let crs_g1 = (0..=degree)
        .map(|i| (g1 * tau.pow(&[i as u64])).into_affine())
        .collect();

    let crs_g2 = (0..=degree)
        .map(|i| (g2 * tau.pow(&[i as u64])).into_affine())
        .collect();

    (crs_g1, crs_g2)
}

fn main() {
    // === 1. Setup ===
    println!("Proving the relation: c = ab + a + b + (a + b)^2");
    println!("Public inputs: a = 2, b = 3, expected c = 36\n");

    let n = 8;
    let omega = Fr::get_root_of_unity(n as u64).unwrap();
    let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();
    let g1 = G1Projective::generator().into_affine();
    let g2 = G2Projective::generator().into_affine();
    let (srs_g1, srs_g2) = get_srs::<Bls12_381>(g1, g2, n + 5);

    // === 2. Circuit Construction ===
    let gates = vec![
        Gate::public_input_gate(fr(1)), // PI a
        Gate::public_input_gate(fr(1)), // PI b
        Gate::public_input_gate(fr(1)), // PI c
        Gate::simple_addition_gate(),   // d = a + b
        Gate::simple_mul_gate(),        // f = a * b
        Gate::simple_mul_gate(),        // e = d * d
        Gate::simple_addition_gate(),   // g = f + d
        Gate::simple_addition_gate(),   // c = g + e
    ];
    println!("✅ Circuit constructed with 8 gates\n");

    let public_inputs = vec![fr(2), fr(3), fr(36)];

    // === 3. Wiring and Witness ===
    let wiring = vec![
        vec![0, 3, 4],       // a
        vec![1, 11, 12],     // b
        vec![2, 23],         // c
        vec![5, 13, 14, 19], // d
        vec![6, 20],         // f
        vec![15, 21],        // e
        vec![7, 22],         // g
    ];

    let witness = Witness {
        // a_wire:   [a=2,   b=3,   c=36,  a=2,   a=2,   d=5,  f=6,   g=11]
        a: vec![fr(2), fr(3), fr(36), fr(2), fr(2), fr(5), fr(6), fr(11)],
        // b_wire:   [0,     0,    0,     b=3,  b=3,   d=5,  d=5,   e=25]
        b: vec![fr(0), fr(0), fr(0), fr(3), fr(3), fr(5), fr(5), fr(25)],
        // c_wire:   [0,     0,    0,     d=5,  f=6,   e=25, g=11,  c=36]
        c: vec![fr(0), fr(0), fr(0), fr(5), fr(6), fr(25), fr(11), fr(36)],
    };

    let circuit = Circuit::new(
        gates,
        witness,
        public_inputs,
        domain.clone(),
        wiring,
        fr(3),
        fr(5),
    );

    // === 4. Prover ===
    let mut rng = test_rng();
    let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();
    let mut prover = KZGProver::<Bls12_381>::new(srs_g1.clone(), domain, g1, false);
    let proof = prover.generate_proof(circuit.clone(), &blinding_scalars);
    println!("📦 Proof successfully generated\n");

    // === 5. Verifier Input ===
    let selector_polys = circuit.get_selector_polynomials();
    let sigma_maps = circuit.permutation.get_sigma_maps();
    let sigma_polys = circuit
        .permutation
        .generate_sigma_polynomials(sigma_maps, &circuit.domain);

    let preprocessed_input = &VerifierPreprocessedInput {
        q_m: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_m, &srs_g1, g1),
        q_l: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_l, &srs_g1, g1),
        q_r: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_r, &srs_g1, g1),
        q_o: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_o, &srs_g1, g1),
        q_c: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_c, &srs_g1, g1),
        sigma_1: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.0, &srs_g1, g1),
        sigma_2: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.1, &srs_g1, g1),
        sigma_3: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.2, &srs_g1, g1),
        x: srs_g2[1],
    };

    println!("🔐 Verifier preprocessed input prepared\n");

    // === 6. Verification ===
    let check = verify_kzg_proof(
        &proof,
        preprocessed_input,
        &circuit.public_inputs,
        &circuit.domain,
        circuit.permutation.k1,
        circuit.permutation.k2,
        g1,
        g2,
    );

    println!(
        "📄 PLONK proof: {}",
        if check {
            "ACCEPTED ✅"
        } else {
            "REJECTED ❌"
        }
    );

    std::process::exit(if check { 0 } else { 1 });
}
