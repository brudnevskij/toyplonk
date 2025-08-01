use crate::circuit::Circuit;
use crate::fft::{compute_lagrange_base, inverse_fft, vec_to_poly};
use crate::prover::Proof;
use crate::transccript::hash_to_field;
use ark_ec::AffineRepr;
use ark_ec::bn::G1Affine;
use ark_ec::pairing::Pairing;
use ark_ff::{Field, One, Zero};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use std::ops::Mul;

#[derive(Debug, Clone)]
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
    // 1-3. validate commitments are in the group

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
    let lagrange_poly_eval = compute_lagrange_base(1, &domain).evaluate(&zeta);

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
    let group_encoded_batch_evals =
        compute_group_encoded_batch_evaluations(r0, v, u, proof.z_omega_bar, proof, g1);

    // 12. batch validation
    let lhs = E::pairing(proof.w_zeta + proof.w_zeta_omega * u, preprocessed_input.x);
    let rhs = E::pairing(
        proof.w_zeta * zeta + proof.w_zeta_omega * u * zeta * domain[1] + full_batch_commitment
            - group_encoded_batch_evals,
        g2,
    );
    lhs == rhs
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
    let lagrange_summand = lagrange_eval.mul(alpha.square());
    let permut_summand = alpha
        * (a_bar + beta * sigma_bar_1 + gamma)
        * (b_bar + beta * sigma_bar_2 + gamma)
        * (c_bar + gamma)
        * z_omega_bar;
    pi - lagrange_summand - permut_summand
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

    let permutation_product = (a_bar + beta * zeta + gamma)
        * (b_bar + beta * k1 * zeta + gamma)
        * (c_bar + beta * k2 * zeta + gamma);
    let permutation_summand_1 =
        z * (permutation_product * alpha + lagrange_eval * alpha.square() + u);

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
    u: E::ScalarField,
    z_omega_bar: E::ScalarField,
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
        + powers_of_v[5] * sigma_bar_2
        + u * z_omega_bar))
        .into()
}

#[cfg(test)]
mod test {
    use crate::circuit::{Circuit, SelectorPolynomials};
    use crate::fft::{compute_lagrange_base, constant_polynomial};
    use crate::gate::Gate;
    use crate::prover::KZGProver;
    use crate::verifier::{
        Challenges, VerifierPreprocessedInput, compute_first_part_of_batched_poly,
        compute_full_batched_polynomial_commitment, compute_group_encoded_batch_evaluations,
        compute_public_input_polynomial, compute_r_constant_terms, verify_kzg_proof,
    };
    use crate::witness::Witness;
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup, Group};
    use ark_ff::{FftField, Field, Zero};
    use ark_poly::Polynomial;
    use ark_poly::univariate::DenseOrSparsePolynomial;
    use ark_std::{UniformRand, test_rng};
    use std::ops::{Add, Mul, Sub};

    fn fr(n: u64) -> Fr {
        Fr::from(n)
    }

    fn dummy_crs<E: Pairing>(
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

    fn dummy_circuit() -> Circuit<Fr> {
        let n = 4;
        let omega = Fr::get_root_of_unity(n as u64).unwrap();
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        let gate0 = Gate::new(fr(1), fr(0), fr(0), fr(0), fr(0)); // a = x₀
        let gate1 = Gate::new(fr(1), fr(0), fr(0), fr(0), fr(0)); // a = x₁
        let gate2 = Gate::new(fr(1), fr(1), fr(1), -fr(1), fr(0)); // a + b + ab = c
        let gate3 = Gate::new(fr(0), fr(0), fr(0), fr(0), fr(0)); // padding (optional)

        let gates = vec![gate0, gate1, gate2, gate3];

        // Witness
        let witness = Witness {
            a: vec![fr(3), fr(4), fr(3), fr(0)], // a[0]=x₀, a[1]=x₁, a[2]=x₀ for gate2
            b: vec![fr(0), fr(0), fr(4), fr(0)], // b[2]=x₁ for gate2
            c: vec![fr(0), fr(0), fr(19), fr(0)], // c[2] = 3 + 4 + 3*4 = 19
        };

        let public_inputs = vec![fr(3), fr(4)];

        Circuit::new(
            gates,
            witness,
            public_inputs,
            domain.clone(),
            vec![vec![1, 6], vec![0, 2]],
            fr(3),
            fr(5),
        )
    }
    #[test]
    fn test_sum_and_product_with_two_public_inputs() {
        let circuit = dummy_circuit();
        let domain = circuit.domain.clone();
        let selector = circuit.get_selector_polynomials();
        let witness = circuit.get_witness_polynomials();
        let pi_poly = circuit.compute_public_input_polynomial();

        // Check constraints
        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(&selector, &witness, &pi_poly));

        let gate_poly = circuit.generate_gate_constraint_polynomial(&selector, &witness, pi_poly);
        let zh = Circuit::vanishing_poly(&domain);
        let gate_dsp = DenseOrSparsePolynomial::from(gate_poly.clone());
        let zh_dsp = DenseOrSparsePolynomial::from(zh);
        let (_, remainder) = gate_dsp.divide_with_q_and_r(&zh_dsp).unwrap();
        assert!(remainder.is_zero());

        let evaluations = domain
            .iter()
            .map(|x| gate_poly.evaluate(x))
            .collect::<Vec<_>>();
        assert!(evaluations.iter().all(|x| x.is_zero()));

        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();
        let (crs_g1, crs_g2) = dummy_crs::<Bls12_381>(g1, g2, circuit.domain.len() + 5);

        let mut prover =
            KZGProver::<Bls12_381>::new(crs_g1.clone(), circuit.domain.clone(), g1, true);

        // Generate blinding scalars (11 required as per the assertion in generate_proof)
        let mut rng = ark_std::test_rng();
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();

        let proof = prover.generate_proof(circuit.clone(), &blinding_scalars);
        //println!("{:?}", proof);

        let ch = Challenges::new(&proof);
        println!("alpha: {}", ch.alpha);
        println!("beta: {}", ch.beta);
        println!("gamma: {}", ch.gamma);
        println!("v: {}", ch.v);
        println!("zeta: {}", ch.zeta);

        let selector_polys = circuit.get_selector_polynomials();
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        let debug = prover.prover_debug_info.unwrap();

        let preprocessed_input = &VerifierPreprocessedInput {
            q_m: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_m, &crs_g1, g1),
            q_l: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_l, &crs_g1, g1),
            q_r: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_r, &crs_g1, g1),
            q_o: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_o, &crs_g1, g1),
            q_c: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_c, &crs_g1, g1),
            sigma_1: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.0, &crs_g1, g1),
            sigma_2: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.1, &crs_g1, g1),
            sigma_3: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.2, &crs_g1, g1),
            x: crs_g2[1],
        };

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

        assert!(check, "verification failed");
    }

    #[test]
    fn test_lineariasation_reconstruction() {
        let circuit = dummy_circuit();
        let selector = circuit.get_selector_polynomials();
        let witness = circuit.get_witness_polynomials();
        let pi_poly = circuit.compute_public_input_polynomial();
        let domain = circuit.domain.clone();

        // Check constraints
        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(&selector, &witness, &pi_poly));

        let gate_poly =
            circuit.generate_gate_constraint_polynomial(&selector, &witness, pi_poly.clone());
        let zh = Circuit::vanishing_poly(&domain);
        let gate_dsp = DenseOrSparsePolynomial::from(gate_poly.clone());
        let zh_dsp = DenseOrSparsePolynomial::from(zh.clone());
        let (_, remainder) = gate_dsp.divide_with_q_and_r(&zh_dsp).unwrap();
        assert!(remainder.is_zero());

        let evaluations = domain
            .iter()
            .map(|x| gate_poly.evaluate(x))
            .collect::<Vec<_>>();
        assert!(evaluations.iter().all(|x| x.is_zero()));

        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();
        let (crs_g1, crs_g2) = dummy_crs::<Bls12_381>(g1, g2, circuit.domain.len() + 5);

        let mut prover =
            KZGProver::<Bls12_381>::new(crs_g1.clone(), circuit.domain.clone(), g1, true);

        // Generate blinding scalars (11 required as per the assertion in generate_proof)
        let mut rng = ark_std::test_rng();
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();

        let proof = prover.generate_proof(circuit.clone(), &blinding_scalars);
        //println!("{:?}", proof);

        let ch = Challenges::new(&proof);
        let selector_polys = circuit.get_selector_polynomials();
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        let debug = prover.prover_debug_info.clone().unwrap();

        let preprocessed_input = &VerifierPreprocessedInput {
            q_m: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_m, &crs_g1, g1),
            q_l: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_l, &crs_g1, g1),
            q_r: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_r, &crs_g1, g1),
            q_o: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_o, &crs_g1, g1),
            q_c: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_c, &crs_g1, g1),
            sigma_1: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.0, &crs_g1, g1),
            sigma_2: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.1, &crs_g1, g1),
            sigma_3: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.2, &crs_g1, g1),
            x: crs_g2[1],
        };

        let lagrange_zeta = compute_lagrange_base(1, &domain).evaluate(&ch.zeta);
        let vanishing_polynomial_eval = Circuit::vanishing_poly(&domain).evaluate(&ch.zeta);
        let r0 = compute_r_constant_terms::<Bls12_381>(
            pi_poly.evaluate(&ch.zeta),
            lagrange_zeta,
            ch.alpha,
            ch.beta,
            ch.gamma,
            proof.a_bar,
            proof.b_bar,
            proof.c_bar,
            proof.sigma_bar_1,
            proof.sigma_bar_2,
            proof.z_omega_bar,
        );

        let d = compute_first_part_of_batched_poly(
            domain.len(),
            circuit.permutation.k1,
            circuit.permutation.k2,
            lagrange_zeta,
            vanishing_polynomial_eval,
            &proof,
            preprocessed_input,
        );
        let f = compute_full_batched_polynomial_commitment::<Bls12_381>(
            ch.v,
            d,
            proof.a,
            proof.b,
            proof.c,
            preprocessed_input.sigma_1,
            preprocessed_input.sigma_2,
        );

        let e =
            compute_group_encoded_batch_evaluations(r0, ch.v, ch.u, proof.z_omega_bar, &proof, g1);

        // check that D + r0 = [r] + u [z]
        assert!(
            debug.linearisation_poly.evaluate(&ch.zeta).is_zero(),
            "r(ζ) must be zero for linearization to work!"
        );
        assert_eq!(
            debug.public_input_poly.evaluate(&ch.zeta),
            compute_public_input_polynomial::<Bls12_381>(&circuit.public_inputs, &domain)
                .evaluate(&ch.zeta)
        );

        assert_eq!(debug.alpha, ch.alpha, "Alpha mismatch!");
        assert_eq!(debug.beta, ch.beta, "Beta mismatch!");
        assert_eq!(debug.gamma, ch.gamma, "Gamma mismatch!");
        assert_eq!(debug.zeta, ch.zeta, "Zeta mismatch!");
        assert_eq!(debug.v, ch.v, "V mismatch!");

        // Recompute D manually to see what's wrong
        let manual_constraint_summand = preprocessed_input.q_m * proof.a_bar * proof.b_bar
            + preprocessed_input.q_l * proof.a_bar
            + preprocessed_input.q_r * proof.b_bar
            + preprocessed_input.q_o * proof.c_bar
            + preprocessed_input.q_c;

        let permutation_product = (proof.a_bar + ch.beta * ch.zeta + ch.gamma)
            * (proof.b_bar + ch.beta * circuit.permutation.k1 * ch.zeta + ch.gamma)
            * (proof.c_bar + ch.beta * circuit.permutation.k2 * ch.zeta + ch.gamma);

        let permutation_summand_1 =
            proof.z * (permutation_product * ch.alpha + lagrange_zeta * ch.alpha.square() + ch.u);

        let permutation_summand_2 = preprocessed_input.sigma_3
            * ((proof.a_bar + ch.beta * proof.sigma_bar_1 + ch.gamma)
                * (proof.b_bar + ch.beta * proof.sigma_bar_2 + ch.gamma)
                * ch.alpha
                * ch.beta
                * proof.z_omega_bar);

        let quotient_summand = (proof.t_lo
            + proof.t_mid * ch.zeta.pow(&[domain.len() as u64])
            + proof.t_hi * ch.zeta.pow(&[2 * domain.len() as u64]))
            * vanishing_polynomial_eval;

        let manual_d = (manual_constraint_summand + permutation_summand_1
            - permutation_summand_2
            - quotient_summand)
            .into_affine();
        assert_eq!(d, manual_d, "D computation mismatch!");

        // individual comparision with r0 components
        let pi = compute_public_input_polynomial::<Bls12_381>(&circuit.public_inputs, &domain);
        assert_eq!(debug.public_input_poly, pi);

        let SelectorPolynomials {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        } = debug.selector_polynomials;
        assert_eq!(
            preprocessed_input.q_l,
            KZGProver::<Bls12_381>::commit_polynomial(&q_l, &crs_g1, g1)
        );
        assert_eq!(
            preprocessed_input.q_r,
            KZGProver::<Bls12_381>::commit_polynomial(&q_r, &crs_g1, g1)
        );
        assert_eq!(
            preprocessed_input.q_m,
            KZGProver::<Bls12_381>::commit_polynomial(&q_m, &crs_g1, g1)
        );
        assert_eq!(
            preprocessed_input.q_o,
            KZGProver::<Bls12_381>::commit_polynomial(&q_o, &crs_g1, g1)
        );
        assert_eq!(
            preprocessed_input.q_c,
            KZGProver::<Bls12_381>::commit_polynomial(&q_c, &crs_g1, g1)
        );

        assert_eq!(proof.a_bar, debug.a_poly_blinded.evaluate(&ch.zeta));
        assert_eq!(proof.b_bar, debug.b_poly_blinded.evaluate(&ch.zeta));
        assert_eq!(proof.c_bar, debug.c_poly_blinded.evaluate(&ch.zeta));

        let test = (preprocessed_input.q_m * proof.a_bar * proof.b_bar)
            + (preprocessed_input.q_l * proof.a_bar)
            + (preprocessed_input.q_r * proof.b_bar)
            + (preprocessed_input.q_o * proof.c_bar)
            + preprocessed_input.q_c
            + (g1 * pi.evaluate(&ch.zeta));

        let lin_constraint_summand = prover.compute_constraint_linearisation_summand(
            proof.a_bar,
            proof.b_bar,
            proof.c_bar,
            &selector_polys,
            &pi,
            ch.zeta,
        );
        let lin_constraint_summand_x =
            KZGProver::<Bls12_381>::commit_polynomial(&lin_constraint_summand, &crs_g1, g1);
        assert_eq!(lin_constraint_summand_x, test.into_affine());

        let lin_quotient_summand = prover.compute_quotient_linearization_summand(
            circuit.gates.len(),
            ch.zeta,
            &zh,
            &debug.t_lo,
            &debug.t_mid,
            &debug.t_hi,
        );
        let lin_quotient_summand_x =
            KZGProver::<Bls12_381>::commit_polynomial(&lin_quotient_summand, &crs_g1, g1);
        assert_eq!(lin_quotient_summand_x, quotient_summand);

        let z_commit = KZGProver::<Bls12_381>::commit_polynomial(&debug.z, &crs_g1, g1);
        assert_eq!(proof.z, z_commit, "Z computation mismatch!");

        let z_omega_zeta_debug = debug.z.evaluate(&(ch.zeta * domain[1]));
        assert_eq!(z_omega_zeta_debug, proof.z_omega_bar, "z(ωζ) mismatch!");

        let r_commit =
            KZGProver::<Bls12_381>::commit_polynomial(&debug.linearisation_poly, &crs_g1, g1);
        let t = d + g1 * r0;
        assert_eq!(d + g1 * r0, r_commit + z_commit * ch.u);

        // CORRECTED: Reconstruct the batched polynomial that should equal [F]₁ - [E]₁
        let batched_poly = debug
            .linearisation_poly // r(X) (full linearized polynomial)
            .add(
                debug
                    .a_poly_blinded
                    .sub(&constant_polynomial(proof.a_bar))
                    .mul(ch.v),
            ) // v·(a(X) - ā)
            .add(
                debug
                    .b_poly_blinded
                    .sub(&constant_polynomial(proof.b_bar))
                    .mul(ch.v.pow(&[2])),
            ) // v²·(b(X) - b̄)
            .add(
                debug
                    .c_poly_blinded
                    .sub(&constant_polynomial(proof.c_bar))
                    .mul(ch.v.pow(&[3])),
            ) // v³·(c(X) - c̄)
            .add(
                debug
                    .sigma_1
                    .sub(&constant_polynomial(proof.sigma_bar_1))
                    .mul(ch.v.pow(&[4])),
            ) // v⁴·(σ₁(X) - σ̄₁)
            .add(
                debug
                    .sigma_2
                    .sub(&constant_polynomial(proof.sigma_bar_2))
                    .mul(ch.v.pow(&[5])),
            ) // v⁵·(σ₂(X) - σ̄₂)
            .add(
                debug
                    .z
                    .sub(&constant_polynomial(proof.z_omega_bar))
                    .mul(ch.u),
            ); // u·(z(X) - z̄ω)

        // Commit to the reconstructed batched polynomial
        let batched_commitment =
            KZGProver::<Bls12_381>::commit_polynomial(&batched_poly, &crs_g1, g1);

        // This should equal [F]₁ - [E]₁
        let f_minus_e = f - e;

        assert_eq!(
            batched_commitment, f_minus_e,
            "Linearization reconstruction failed!"
        );
    }

    #[test]
    fn test_permutation_commitment_reconstruction() {
        let circuit = dummy_circuit();
        let domain = circuit.domain.clone();
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();
        let (crs_g1, crs_g2) = dummy_crs::<Bls12_381>(g1, g2, circuit.domain.len() + 5);

        let mut prover =
            KZGProver::<Bls12_381>::new(crs_g1.clone(), circuit.domain.clone(), g1, true);

        // Generate blinding scalars (11 required as per the assertion in generate_proof)
        let mut rng = ark_std::test_rng();
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();
        let proof = prover.generate_proof(circuit.clone(), &blinding_scalars);

        let ch = Challenges::new(&proof);
        let selector_polys = circuit.get_selector_polynomials();
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        let debug = prover.prover_debug_info.clone().unwrap();

        let preprocessed_input = &VerifierPreprocessedInput::<Bls12_381> {
            q_m: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_m, &crs_g1, g1),
            q_l: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_l, &crs_g1, g1),
            q_r: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_r, &crs_g1, g1),
            q_o: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_o, &crs_g1, g1),
            q_c: KZGProver::<Bls12_381>::commit_polynomial(&selector_polys.q_c, &crs_g1, g1),
            sigma_1: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.0, &crs_g1, g1),
            sigma_2: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.1, &crs_g1, g1),
            sigma_3: KZGProver::<Bls12_381>::commit_polynomial(&sigma_polys.2, &crs_g1, g1),
            x: crs_g2[1],
        };

        let pi =
            compute_public_input_polynomial::<Bls12_381>(&circuit.public_inputs, &circuit.domain);
        let pi_zeta = pi.evaluate(&ch.zeta);
        let lagrange_base_1 = compute_lagrange_base(1, &domain);
        let lagrange_base_1_zeta = lagrange_base_1.evaluate(&ch.zeta);

        // computing the permutation
        let k1 = circuit.permutation.k1.clone();
        let k2 = circuit.permutation.k2.clone();

        let permut_product = (proof.a_bar + ch.beta * ch.zeta + ch.gamma)
            * (proof.b_bar + ch.beta * k1 * ch.zeta + ch.gamma)
            * (proof.c_bar + ch.beta * k2 * ch.zeta + ch.gamma);
        let lin_permut_summand_1 = debug.z.mul(permut_product);
        let lin_permut_summand_1_commit =
            KZGProver::<Bls12_381>::commit_polynomial(&lin_permut_summand_1, &crs_g1, g1);

        let d_permutation_summand_1 =
            proof.z * (permut_product * ch.alpha + lagrange_base_1_zeta * ch.alpha.square() + ch.u);

        assert_eq!(
            lin_permut_summand_1_commit * ch.alpha
                + proof.z * (lagrange_base_1_zeta * ch.alpha.square() + ch.u),
            d_permutation_summand_1
        );

        // (c(z) + beta * S(x) + gamma )
        // because c_bar is constant I add it to gamma
        let permut_2_c = debug.sigma_3.mul(ch.beta) + constant_polynomial(proof.c_bar + ch.gamma);
        let lin_permut_summand_2 = permut_2_c.mul(
            (proof.a_bar + ch.beta * proof.sigma_bar_1 + ch.gamma)
                * (proof.b_bar + ch.beta * proof.sigma_bar_2 + ch.gamma)
                * proof.z_omega_bar,
        );
        let lin_permut_summand_2_commit =
            KZGProver::<Bls12_381>::commit_polynomial(&lin_permut_summand_2, &crs_g1, g1);

        let d_permutation_summand_2 = preprocessed_input.sigma_3
            * ((proof.a_bar + ch.beta * proof.sigma_bar_1 + ch.gamma)
                * (proof.b_bar + ch.beta * proof.sigma_bar_2 + ch.gamma)
                * ch.alpha
                * ch.beta
                * proof.z_omega_bar);
        let r0_2_permut_portion = (proof.a_bar + ch.beta * proof.sigma_bar_1 + ch.gamma)
            * (proof.b_bar + ch.beta * proof.sigma_bar_2 + ch.gamma)
            * (proof.c_bar + ch.gamma)
            * ch.alpha
            * proof.z_omega_bar;

        assert_eq!(
            lin_permut_summand_2_commit * ch.alpha,
            d_permutation_summand_2 + g1 * r0_2_permut_portion
        );

        let lin_init_z_summand = prover.compute_init_z_linearization_summand(
            ch.alpha,
            ch.zeta,
            &debug.z,
            &lagrange_base_1,
        );
        let lin_init_z_summand_commit =
            KZGProver::<Bls12_381>::commit_polynomial(&lin_init_z_summand, &crs_g1, g1);
        assert_eq!(
            proof.z * lagrange_base_1_zeta * ch.alpha.square()
                - g1 * lagrange_base_1_zeta * ch.alpha.square(),
            lin_init_z_summand_commit
        );

        let r0_permut = g1 * (-lagrange_base_1_zeta * ch.alpha.square() - r0_2_permut_portion);
        let full_permutation_lin = (lin_permut_summand_1_commit - lin_permut_summand_2_commit)
            * ch.alpha
            + lin_init_z_summand_commit;
        assert_eq!(
            full_permutation_lin + proof.z * ch.u, // ← This is needed to match the verifier's computation
            d_permutation_summand_1 - d_permutation_summand_2 + r0_permut
        );

        let vanishing_polynomial = Circuit::vanishing_poly(&domain);
        let d_constraint_summand = preprocessed_input.q_m * proof.a_bar * proof.b_bar
            + preprocessed_input.q_l * proof.a_bar
            + preprocessed_input.q_r * proof.b_bar
            + preprocessed_input.q_o * proof.c_bar
            + preprocessed_input.q_c;
        let d_quotient_summand = (proof.t_lo
            + proof.t_mid * ch.zeta.pow(&[domain.len() as u64])
            + proof.t_hi * ch.zeta.pow(&[2 * domain.len() as u64]))
            * vanishing_polynomial.evaluate(&ch.zeta);

        let manual_d = d_constraint_summand + d_permutation_summand_1
            - d_permutation_summand_2
            - d_quotient_summand;
        let d = compute_first_part_of_batched_poly(
            4,
            circuit.permutation.k1,
            circuit.permutation.k2,
            lagrange_base_1.evaluate(&ch.zeta),
            vanishing_polynomial.evaluate(&ch.zeta),
            &proof,
            &preprocessed_input,
        );
        assert_eq!(d, manual_d);

        let lin_constraint = prover.compute_constraint_linearisation_summand(
            proof.a_bar,
            proof.b_bar,
            proof.c_bar,
            &debug.selector_polynomials,
            &pi,
            ch.zeta,
        );
        let lin_quot = prover.compute_quotient_linearization_summand(
            4,
            ch.zeta,
            &vanishing_polynomial,
            &debug.t_lo,
            &debug.t_mid,
            &debug.t_hi,
        );

        let r_manual = (lin_constraint + lin_permut_summand_1.mul(ch.alpha) + lin_init_z_summand)
            .sub(&lin_permut_summand_2.mul(ch.alpha))
            .sub(&lin_quot);
        let r_manual_commit = KZGProver::<Bls12_381>::commit_polynomial(&r_manual, &crs_g1, g1);
        assert_eq!(r_manual, debug.linearisation_poly);

        assert_eq!(
            r_manual_commit + proof.z * ch.u,
            manual_d + r0_permut + g1 * pi_zeta
        );
    }
}
