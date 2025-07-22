use crate::circuit::Circuit;
use crate::fft::{inverse_fft, vec_to_poly};
use crate::prover::Proof;
use crate::transccript::hash_to_field;
use ark_ec::AffineRepr;
use ark_ec::bn::G1Affine;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, Zero};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;

struct Challenges<E: Pairing> {
    alpha: E::ScalarField,
    beta: E::ScalarField,
    gamma: E::ScalarField,
    zeta: E::ScalarField,
    v: E::ScalarField,
    u: E::ScalarField,
}

impl<E: Pairing> Challenges<E> {
    pub fn new(proof: &Proof<E>) -> Self {
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

        Challenges {
            alpha,
            beta,
            gamma,
            zeta,
            v,
            u,
        }
    }
}

pub struct VerifierPreprocessedInput<E: Pairing> {
    q_m: E::G1Affine,
    q_l: E::G1Affine,
    q_r: E::G1Affine,
    q_o: E::G1Affine,
    q_c: E::G1Affine,
    sigma_1: E::G1Affine,
    sigma_2: E::G1Affine,
    sigma_3: E::G1Affine,
    x: E::G2Affine,
}

pub fn verify_kzg_proof<E: Pairing>(
    proof: &Proof<E>,
    preprocessed_input: &VerifierPreprocessedInput<E>,
    public_input: &[E::ScalarField],
    domain: &[E::ScalarField],
    k1: E::ScalarField,
    k2: E::ScalarField,
    g1: E::G1Affine,
    g2: E::G2Affine,
) -> bool {
    // 1 - 3, validation is mostly done by type system

    // 4. challenges
    let Challenges {
        alpha,
        beta,
        gamma,
        zeta,
        v,
        u,
    } = Challenges::new(proof);

    // 5. zero polynomial eval
    let vanishing_polynomial = Circuit::vanishing_poly(domain);
    let vanishing_eval = vanishing_polynomial.evaluate(&zeta);

    // 6. Lagrange poly eval
    let lagrange_poly_eval = compute_lagrange_polynomial_eval::<E>(domain.len(), domain[1], zeta);

    // 7. PI evaluation
    let pi = compute_public_input_polynomial::<E>(public_input, domain);
    let pi_eval = pi.evaluate(&zeta);

    // 8. compute constant part of r(X)
    let r0 = compute_r_constant_terms::<E>(
        pi_eval,
        lagrange_poly_eval,
        alpha,
        beta,
        gamma,
        proof.a_bar,
        proof.b_bar,
        proof.c_bar,
        proof.sigma_bar_1,
        proof.sigma_bar_2,
        proof.z_omega_bar,
    );

    // 9. compute first part of batched PC
    let first_part_commitment = compute_first_part_of_batched_poly(
        domain.len(),
        k1,
        k2,
        lagrange_poly_eval,
        vanishing_eval,
        proof,
        preprocessed_input,
    );

    // 10. full batched PC
    let full_batch_commitment = compute_full_batched_polynomial_commitment::<E>(
        v,
        first_part_commitment,
        proof.a,
        proof.b,
        proof.c,
        preprocessed_input.sigma_1,
        preprocessed_input.sigma_2,
    );

    // 11. group encoded batch evaluations
    let group_encoded_batch_evals = compute_group_encoded_batch_evaluations(r0, v, proof, g1);

    // 12. batch validation
    let lhs = E::pairing(proof.w_zeta + proof.w_zeta_omega * u, preprocessed_input.x);
    let rhs = E::pairing(
        proof.w_zeta * zeta + proof.w_zeta_omega * u * zeta * domain[1] + full_batch_commitment - group_encoded_batch_evals,
        g2,
    );

    lhs == rhs
}

fn compute_lagrange_polynomial_eval<E: Pairing>(
    n: usize,
    omega: E::ScalarField,
    zeta: E::ScalarField,
) -> E::ScalarField {
    // w(z^n - 1) / n(zeta-omega)
    let numerator = omega * (zeta.pow(&[n as u64]) - E::ScalarField::one());
    let denominator = E::ScalarField::from(n as u64) * (zeta - omega);
    numerator / denominator
}

pub fn compute_public_input_polynomial<E: Pairing>(
    public_inputs: &[E::ScalarField],
    domain: &[E::ScalarField],
) -> DensePolynomial<E::ScalarField> {
    let mut evaluations = vec![E::ScalarField::zero(); domain.len()];
    for (i, &x) in public_inputs.iter().enumerate() {
        evaluations[i] = -x;
    }

    vec_to_poly(inverse_fft(&evaluations, domain[1]))
}

fn compute_r_constant_terms<E: Pairing>(
    pi: E::ScalarField,
    lagrange_eval: E::ScalarField,
    alpha: E::ScalarField,
    beta: E::ScalarField,
    gamma: E::ScalarField,
    a_bar: E::ScalarField,
    b_bar: E::ScalarField,
    c_bar: E::ScalarField,
    sigma_bar_1: E::ScalarField,
    sigma_bar_2: E::ScalarField,
    z_omega_bar: E::ScalarField,
) -> E::ScalarField {
    pi - (lagrange_eval * alpha.square())
        - alpha
            * (a_bar + beta * sigma_bar_1 + gamma)
            * (b_bar + beta * sigma_bar_2 + gamma)
            * (c_bar + gamma)
            * z_omega_bar
}

fn compute_first_part_of_batched_poly<E: Pairing>(
    n: usize,
    k1: E::ScalarField,
    k2: E::ScalarField,
    lagrange_eval: E::ScalarField,
    vanishing_polynomial_eval: E::ScalarField,
    proof: &Proof<E>,
    verifier_preprocessed_input: &VerifierPreprocessedInput<E>,
) -> E::G1Affine {
    let &Proof {
        z,
        t_lo,
        t_mid,
        t_hi,
        a_bar,
        b_bar,
        c_bar,
        sigma_bar_1,
        sigma_bar_2,
        z_omega_bar,
        ..
    } = proof;
    let &VerifierPreprocessedInput {
        q_m,
        q_l,
        q_r,
        q_o,
        q_c,
        sigma_3,
        ..
    } = verifier_preprocessed_input;
    let Challenges {
        alpha,
        beta,
        gamma,
        zeta,
        u,
        ..
    } = Challenges::new(&proof);

    let constraint_system_summand = q_m.into_group() * a_bar * b_bar
        + q_l.into_group() * a_bar
        + q_r.into_group() * b_bar
        + q_o.into_group() * c_bar
        + q_c.into_group();

    let permutation_summand_1 = z
        * ((a_bar + beta * zeta + gamma)
            * (b_bar + beta * k1 * zeta + gamma)
            * (c_bar + beta * k2 * zeta + gamma)
            * alpha
            + lagrange_eval * alpha.square()
            + u);

    let permutation_summand_2 = sigma_3
        * ((a_bar + beta * sigma_bar_1 + gamma)
            * (b_bar + beta * sigma_bar_2 + gamma)
            * alpha
            * beta
            * z_omega_bar);

    let quotient_summand =
        (t_lo + t_mid * zeta.pow(&[n as u64]) + t_hi * zeta.pow(&[2 * n as u64]))
            * vanishing_polynomial_eval;

    (constraint_system_summand + permutation_summand_1 - permutation_summand_2 - quotient_summand)
        .into()
}

fn compute_full_batched_polynomial_commitment<E: Pairing>(
    v: E::ScalarField,
    d: E::G1Affine,
    a: E::G1Affine,
    b: E::G1Affine,
    c: E::G1Affine,
    sigma_1: E::G1Affine,
    sigma_2: E::G1Affine,
) -> E::G1Affine {
    let powers_of_v = (0..=5).map(|i| v.pow(&[i as u64])).collect::<Vec<_>>();

    (d + a * powers_of_v[1]
        + b * powers_of_v[2]
        + c * powers_of_v[3]
        + sigma_1 * powers_of_v[4]
        + sigma_2 * powers_of_v[5])
        .into()
}

fn compute_group_encoded_batch_evaluations<E: Pairing>(
    r0: E::ScalarField,
    v: E::ScalarField,
    proof: &Proof<E>,
    g1: E::G1Affine,
) -> E::G1Affine {
    let &Proof {
        a_bar,
        b_bar,
        c_bar,
        sigma_bar_1,
        sigma_bar_2,
        ..
    } = proof;
    let powers_of_v = (0..=5).map(|i| v.pow(&[i as u64])).collect::<Vec<_>>();
    (g1 * (-r0
        + powers_of_v[1] * a_bar
        + powers_of_v[2] * b_bar
        + powers_of_v[3] * c_bar
        + powers_of_v[4] * sigma_bar_1
        + powers_of_v[5] * sigma_bar_2))
        .into()
}
