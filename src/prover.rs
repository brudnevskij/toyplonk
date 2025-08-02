use crate::circuit::{Circuit, SelectorPolynomials};
use crate::fft::{compute_lagrange_base, constant_polynomial, vec_to_poly};
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, One, PrimeField, Zero};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::{DenseUVPolynomial, Polynomial};
use std::ops::{Add, Div, Mul, Sub};
use ark_serialize::CanonicalSerialize;
use sha2::{Digest, Sha256};

/// A PLONK prover using KZG commitments.
/// Handles constraint construction, permutation checks, quotient polynomial,
/// Fiat–Shamir challenge generation, and opening proof construction.
pub struct KZGProver<E: Pairing> {
    pub crs: Vec<E::G1Affine>,
    pub domain: Vec<E::ScalarField>,
    pub g1: E::G1Affine,
    debug_mode: bool,
    pub prover_debug_info: Option<ProverDebugInfo<E>>,
}

/// Central proof structure output by the prover.
/// Contains all commitments and evaluation values needed by the verifier.
#[derive(Clone, Debug)]
pub struct Proof<E: Pairing> {
    pub a: E::G1Affine,
    pub b: E::G1Affine,
    pub c: E::G1Affine,
    pub z: E::G1Affine,
    pub t_lo: E::G1Affine,
    pub t_mid: E::G1Affine,
    pub t_hi: E::G1Affine,
    pub w_zeta: E::G1Affine,
    pub w_zeta_omega: E::G1Affine,
    pub a_bar: E::ScalarField,
    pub b_bar: E::ScalarField,
    pub c_bar: E::ScalarField,
    pub sigma_bar_1: E::ScalarField,
    pub sigma_bar_2: E::ScalarField,
    pub z_omega_bar: E::ScalarField,
}

/// Debug information to inspect prover internals during execution.
#[derive(Debug, Clone)]
pub struct ProverDebugInfo<E: Pairing> {
    pub linearisation_poly: DensePolynomial<E::ScalarField>,

    // quotient poly components
    pub constraint_summand: DensePolynomial<E::ScalarField>,
    pub permutation_summand: DensePolynomial<E::ScalarField>,
    pub init_z_summand: DensePolynomial<E::ScalarField>,
    pub quotient_polynomial: DensePolynomial<E::ScalarField>,
    pub t_lo: DensePolynomial<E::ScalarField>,
    pub t_mid: DensePolynomial<E::ScalarField>,
    pub t_hi: DensePolynomial<E::ScalarField>,

    // opening
    pub opening_poly: DensePolynomial<E::ScalarField>,
    pub opening_omega_poly: DensePolynomial<E::ScalarField>,

    // challenges
    pub alpha: E::ScalarField,
    pub beta: E::ScalarField,
    pub gamma: E::ScalarField,
    pub zeta: E::ScalarField,
    pub v: E::ScalarField,

    // wire polys
    pub a_poly_blinded: DensePolynomial<E::ScalarField>,
    pub b_poly_blinded: DensePolynomial<E::ScalarField>,
    pub c_poly_blinded: DensePolynomial<E::ScalarField>,

    pub selector_polynomials: SelectorPolynomials<E::ScalarField>,

    // permutation
    pub z: DensePolynomial<E::ScalarField>,
    pub sigma_1: DensePolynomial<E::ScalarField>,
    pub sigma_2: DensePolynomial<E::ScalarField>,
    pub sigma_3: DensePolynomial<E::ScalarField>,

    pub public_input_poly: DensePolynomial<E::ScalarField>,
}

impl<E: Pairing> KZGProver<E> {
    pub fn new(
        crs: Vec<E::G1Affine>,
        domain: Vec<E::ScalarField>,
        g1: E::G1Affine,
        debug_mode: bool,
    ) -> Self {
        Self {
            crs,
            domain,
            g1,
            debug_mode,
            prover_debug_info: None,
        }
    }

    /// Generate a PLONK proof using the circuit and witness.
    /// This is the main function responsible for performing all 5 rounds.
    pub fn generate_proof(
        &mut self,
        circuit: Circuit<E::ScalarField>,
        blinding_scalars: &[E::ScalarField],
    ) -> Proof<E> {
        let n = circuit.gates.len();
        assert_eq!(blinding_scalars.len(), 11);
        assert_eq!(self.crs.len(), n + 5 + 1);

        let mut commitment_buffer = Vec::new();
        let vanishing_poly = Circuit::vanishing_poly(&self.domain);
        let witness_polynomial = circuit.get_witness_polynomials();

        // === Round 1: Commit to blinded wire polynomials ===
        let a = witness_polynomial.a.clone();
        let b = witness_polynomial.b.clone();
        let c = witness_polynomial.c.clone();

        let a_poly = Self::compute_wire_coefficients_form(
            blinding_scalars[0],
            blinding_scalars[1],
            &a,
            &vanishing_poly,
        );
        let a_commitment = Self::commit_polynomial(&a_poly, &self.crs, self.g1);
        commitment_buffer.push(a_commitment);

        let b_poly = Self::compute_wire_coefficients_form(
            blinding_scalars[2],
            blinding_scalars[3],
            &b,
            &vanishing_poly,
        );
        let b_commitment = Self::commit_polynomial(&b_poly, &self.crs, self.g1);
        commitment_buffer.push(b_commitment);

        let c_poly = Self::compute_wire_coefficients_form(
            blinding_scalars[4],
            blinding_scalars[5],
            &c,
            &vanishing_poly,
        );
        let c_commitment = Self::commit_polynomial(&c_poly, &self.crs, self.g1);
        commitment_buffer.push(c_commitment);

        // === Round 2: Compute permutation product and commit to Z ===V
        let beta = hash_to_field("beta", &commitment_buffer);
        let gamma = hash_to_field("gamma", &commitment_buffer);
        let rolling_product = circuit
            .permutation
            .get_rolling_product(gamma, beta, &self.domain);
        let z = Self::compute_permutation_polynomial(
            blinding_scalars[6],
            blinding_scalars[7],
            blinding_scalars[8],
            &vanishing_poly,
            &rolling_product,
        );

        let z_commitment = Self::commit_polynomial(&z, &self.crs, self.g1);
        commitment_buffer.push(z_commitment);

        // === Round 3: Compute quotient polynomial and split into parts ===
        // first summand
        let public_input_poly = circuit.compute_public_input_polynomial();
        let selector_polynomials = circuit.get_selector_polynomials();
        let constraint_summand = self.compute_constraint_summand(
            &a_poly,
            &b_poly,
            &c_poly,
            public_input_poly.clone(),
            &vanishing_poly,
            &selector_polynomials,
        );

        // second summand, permutation 2
        let alpha = hash_to_field("alpha", &commitment_buffer);
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polynomials = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &self.domain);

        let lagrange_base_1 = compute_lagrange_base(1, &self.domain);
        let init_z_summand =
            self.compute_init_z_summand(alpha, &vanishing_poly, &z, &lagrange_base_1);

        let permutation_summand = self.compute_permutation_summand(
            &a_poly,
            &b_poly,
            &c_poly,
            &sigma_polynomials,
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            &vanishing_poly,
            alpha,
        );

        let quotient_polynomial =
            constraint_summand.clone() + permutation_summand.clone() + init_z_summand.clone();

        let (t_lo, t_mid, t_hi) = Self::split_quotient_polynomial(
            &quotient_polynomial,
            blinding_scalars[9],
            blinding_scalars[10],
            self.domain.len(),
        );

        let t_lo_commitment = Self::commit_polynomial(&t_lo, &self.crs, self.g1);
        let t_mid_commitment = Self::commit_polynomial(&t_mid, &self.crs, self.g1);
        let t_hi_commitment = Self::commit_polynomial(&t_hi, &self.crs, self.g1);
        commitment_buffer.push(t_lo_commitment);
        commitment_buffer.push(t_mid_commitment);
        commitment_buffer.push(t_hi_commitment);

        // === Round 4: Evaluate wires and sigmas at ζ ===
        let zeta = hash_to_field("zeta", &commitment_buffer);

        let a_bar = a_poly.evaluate(&zeta);
        let b_bar = b_poly.evaluate(&zeta);
        let c_bar = c_poly.evaluate(&zeta);

        let sigma_bar_1 = sigma_polynomials.0.evaluate(&zeta);
        let sigma_bar_2 = sigma_polynomials.1.evaluate(&zeta);

        let z_omega_bar = z.evaluate(&(zeta * self.domain[1]));

        // === Round 5: Construct linearization poly and opening proofs ===
        // previous round's output could be added to transcript
        // however, I will just add v for now.
        let v = hash_to_field("v", &commitment_buffer);

        let linearisation_polynomial = self.compute_linearisation_polynomial(
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
            alpha,
            beta,
            gamma,
            zeta,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &selector_polynomials,
            &public_input_poly,
            &z,
            &sigma_polynomials.2,
            &lagrange_base_1,
            &vanishing_poly,
            &t_lo,
            &t_mid,
            &t_hi,
        );

        let w_zeta = self.compute_opening_proof_polynomial(
            v,
            zeta,
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            &linearisation_polynomial,
            &a_poly,
            &b_poly,
            &c_poly,
            &sigma_polynomials.0,
            &sigma_polynomials.1,
        );

        let w_zeta_commitment = Self::commit_polynomial(&w_zeta, &self.crs, self.g1);
        commitment_buffer.push(w_zeta_commitment);

        let w_zeta_omega =
            self.compute_opening_proof_polynomial_omega(self.domain[1], zeta, z_omega_bar, &z);

        let w_zeta_omega_commitment = Self::commit_polynomial(&w_zeta_omega, &self.crs, self.g1);
        commitment_buffer.push(w_zeta_omega_commitment);

        if self.debug_mode {
            self.prover_debug_info = Some(ProverDebugInfo {
                linearisation_poly: linearisation_polynomial,
                constraint_summand,
                permutation_summand,
                init_z_summand,

                quotient_polynomial,
                t_lo,
                t_mid,
                t_hi,

                opening_poly: w_zeta,
                opening_omega_poly: w_zeta_omega,

                alpha,
                beta,
                gamma,
                zeta,
                v,

                a_poly_blinded: a_poly,
                b_poly_blinded: b_poly,
                c_poly_blinded: c_poly,

                selector_polynomials,

                z,
                sigma_1: sigma_polynomials.0,
                sigma_2: sigma_polynomials.1,
                sigma_3: sigma_polynomials.2,

                public_input_poly,
            });
        }

        Proof {
            a: a_commitment,
            b: b_commitment,
            c: c_commitment,
            z: z_commitment,
            t_lo: t_lo_commitment,
            t_mid: t_mid_commitment,
            t_hi: t_hi_commitment,
            w_zeta: w_zeta_commitment,
            w_zeta_omega: w_zeta_omega_commitment,
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
        }
    }

    /// Constructs the arithmetic constraint polynomial part of the quotient
    ///
    /// Form: q_M·a·b + q_L·a + q_R·b + q_O·c + q_C + PI(X)
    ///
    /// The sum is then divided by the vanishing polynomial Z_H(X)
    fn compute_constraint_summand(
        &self,
        blinded_a: &DensePolynomial<E::ScalarField>,
        blinded_b: &DensePolynomial<E::ScalarField>,
        blinded_c: &DensePolynomial<E::ScalarField>,
        pi: DensePolynomial<E::ScalarField>,
        vanishing_poly: &DensePolynomial<E::ScalarField>,
        selector_polynomials: &SelectorPolynomials<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let SelectorPolynomials {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        } = selector_polynomials;

        let blinded_constraint_poly = blinded_a.naive_mul(&blinded_b).naive_mul(q_m)
            + blinded_a.naive_mul(q_l)
            + blinded_b.naive_mul(q_r)
            + blinded_c.naive_mul(q_o)
            + pi
            + q_c.clone();
        let blinded_constraint_poly = DenseOrSparsePolynomial::from(blinded_constraint_poly);

        let (constraint_summand, r) = blinded_constraint_poly
            .divide_with_q_and_r(&vanishing_poly.into())
            .unwrap();
        assert!(
            r.is_zero(),
            "gate constraint has non zero remainder: {:?}",
            r
        );
        constraint_summand
    }

    /// Evaluate all components at ζ and bundle them into an opening proof.
    /// This is a batch opening at the evaluation point.
    fn compute_opening_proof_polynomial(
        &self,
        v: E::ScalarField,
        zeta: E::ScalarField,
        a_bar: E::ScalarField,
        b_bar: E::ScalarField,
        c_bar: E::ScalarField,
        sigma_bar_1: E::ScalarField,
        sigma_bar_2: E::ScalarField,
        r: &DensePolynomial<E::ScalarField>,
        a: &DensePolynomial<E::ScalarField>,
        b: &DensePolynomial<E::ScalarField>,
        c: &DensePolynomial<E::ScalarField>,
        sigma_1: &DensePolynomial<E::ScalarField>,
        sigma_2: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let powers_of_v = (0..=5).map(|i| v.pow(&[i as u64])).collect::<Vec<_>>();

        (r.clone()
            + (a - &constant_polynomial(a_bar)).mul(powers_of_v[1])
            + (b - &constant_polynomial(b_bar)).mul(powers_of_v[2])
            + (c - &constant_polynomial(c_bar)).mul(powers_of_v[3])
            + (sigma_1 - &constant_polynomial(sigma_bar_1)).mul(powers_of_v[4])
            + (sigma_2 - &constant_polynomial(sigma_bar_2)).mul(powers_of_v[5]))
        // (X-z)
        .div(&vec_to_poly(vec![-zeta, E::ScalarField::one()]))
    }

    /// Evaluate Z at ζw and bundle it into an opening proof.
    fn compute_opening_proof_polynomial_omega(
        &self,
        omega: E::ScalarField,
        zeta: E::ScalarField,
        z_omega_bar: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        // (z(x) - z(ζ)) / X - ζw
        (z - &constant_polynomial(z_omega_bar))
            .div(&vec_to_poly(vec![-zeta * omega, E::ScalarField::one()]))
    }

    /// Constructs the linearization polynomial `r(X)`
    ///
    /// This polynomial combines all constraint evaluations at ζ (zeta), including:
    /// - Gate constraints (q_L · ā + q_R · b̄ + q_M · āb̄ + ...)
    /// - Permutation constraints (Z(X) and σ terms)
    /// - Initialization constraint at L₁(ζ)
    /// - Quotient polynomial evaluation at ζ multiplied by vanishing polynomial
    pub fn compute_linearisation_polynomial(
        &self,
        a_bar: E::ScalarField,
        b_bar: E::ScalarField,
        c_bar: E::ScalarField,
        sigma_bar_1: E::ScalarField,
        sigma_bar_2: E::ScalarField,
        z_omega_bar: E::ScalarField,
        alpha: E::ScalarField,
        beta: E::ScalarField,
        gamma: E::ScalarField,
        zeta: E::ScalarField,
        k1: E::ScalarField,
        k2: E::ScalarField,
        selector_polynomials: &SelectorPolynomials<E::ScalarField>,
        public_input_polynomial: &DensePolynomial<E::ScalarField>,
        z: &DensePolynomial<E::ScalarField>,
        sigma_3: &DensePolynomial<E::ScalarField>,
        lagrange_base_1: &DensePolynomial<E::ScalarField>,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
        t_lo: &DensePolynomial<E::ScalarField>,
        t_mid: &DensePolynomial<E::ScalarField>,
        t_hi: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let constraint_lin_summand = self.compute_constraint_linearisation_summand(
            a_bar,
            b_bar,
            c_bar,
            selector_polynomials,
            public_input_polynomial,
            zeta,
        );

        let permutation_lin_summand = self.compute_permutation_linearization_summand(
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
            alpha,
            beta,
            gamma,
            zeta,
            k1,
            k2,
            z,
            sigma_3,
        );

        let init_z_lin_summand =
            self.compute_init_z_linearization_summand(alpha, zeta, z, lagrange_base_1);

        let quotient_summand = self.compute_quotient_linearization_summand(
            self.domain.len(),
            zeta,
            vanishing_polynomial,
            t_lo,
            t_mid,
            t_hi,
        );

        (constraint_lin_summand + permutation_lin_summand + init_z_lin_summand)
            .sub(&quotient_summand)
    }

    /// Reconstructs t(ζ) * Z_H(ζ), where t = t_lo + ζⁿ t_mid + ζ²ⁿ t_hi
    ///
    /// This term is subtracted in the linearization polynomial to enforce
    /// that the quotient polynomial reconstruction is valid.
    pub fn compute_quotient_linearization_summand(
        &self,
        n: usize,
        zeta: E::ScalarField,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
        t_lo: &DensePolynomial<E::ScalarField>,
        t_mid: &DensePolynomial<E::ScalarField>,
        t_hi: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let vanishing_z = vanishing_polynomial.evaluate(&zeta);
        let zeta_n = zeta.pow(&[n as u64]);
        let zeta_2n = zeta_n.square();

        (t_lo.clone() + t_mid.mul(zeta_n) + t_hi.mul(zeta_2n)).mul(vanishing_z)
    }

    /// Computes the initialization term `α² · (Z(ζ) - 1) · L₁(ζ)`
    ///
    /// Enforces that the permutation product polynomial Z(X)
    /// starts with value 1 at the first position (i.e., Z(ω⁰) = 1)
    pub fn compute_init_z_linearization_summand(
        &self,
        alpha: E::ScalarField,
        zeta: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
        lagrange_base_1: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let lagrange_base_evaluation = lagrange_base_1.evaluate(&zeta);
        // Z(x) - 1
        z.sub(&constant_polynomial(E::ScalarField::one()))
            // (Z(x) - 1)L1(z)
            .mul(lagrange_base_evaluation)
            // a^2 *[(Z(x) - 1)L1(z)]
            .mul(alpha.square())
    }

    /// Constructs the permutation argument part of the linearization polynomial.
    ///
    ///   α · [Z(ζ) · Π(lhs terms) − Z(ζ·ω) · Π(rhs terms)]
    ///
    /// lhs = (ā + β·ζ + γ)(b̄ + β·k₁·ζ + γ)(c̄ + β·k₂·ζ + γ)  
    /// rhs = (ā + β·σ₁(ζ) + γ)(b̄ + β·σ₂(ζ) + γ)(c̄ + β·σ₃(ζ) + γ)
    ///
    /// This linear term is included in the batched polynomial opening proof at ζ.
    pub fn compute_permutation_linearization_summand(
        &self,
        a_bar: E::ScalarField,
        b_bar: E::ScalarField,
        c_bar: E::ScalarField,
        sigma_bar_1: E::ScalarField,
        sigma_bar_2: E::ScalarField,
        z_omega_bar: E::ScalarField,
        alpha: E::ScalarField,
        beta: E::ScalarField,
        gamma: E::ScalarField,
        zeta: E::ScalarField,
        k1: E::ScalarField,
        k2: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
        sigma_3: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let permut_1 = (a_bar + beta * zeta + gamma)
            * (b_bar + beta * k1 * zeta + gamma)
            * (c_bar + beta * k2 * zeta + gamma);
        let permut_summand_1 = z.mul(permut_1);

        // (c(z) + beta * S(x) + gamma )
        // because c_bar is constant I add it to gamma
        let permut_2_c = sigma_3.mul(beta) + constant_polynomial(c_bar + gamma);
        let permut_summand_2 = permut_2_c.mul(
            (a_bar + beta * sigma_bar_1 + gamma)
                * (b_bar + beta * sigma_bar_2 + gamma)
                * z_omega_bar,
        );

        (permut_summand_1.sub(&permut_summand_2)).mul(alpha)
    }

    /// Evaluates arithmetic gate constraints at ζ
    ///
    /// r_gate(ζ) = q_L(ζ) · ā + q_R(ζ) · b̄ + q_M(ζ) · ā·b̄ + q_O(ζ) · c̄ + q_C(ζ) + PI(ζ)
    pub fn compute_constraint_linearisation_summand(
        &self,
        a_bar: E::ScalarField,
        b_bar: E::ScalarField,
        c_bar: E::ScalarField,
        selector_polynomials: &SelectorPolynomials<E::ScalarField>,
        public_polynomial: &DensePolynomial<E::ScalarField>,
        zeta: E::ScalarField,
    ) -> DensePolynomial<E::ScalarField> {
        let pi_eval = constant_polynomial(public_polynomial.evaluate(&zeta));
        let SelectorPolynomials {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        } = selector_polynomials;
        let l_term = q_l.mul(a_bar);
        let r_term = q_r.mul(b_bar);
        let m_term = q_m.mul(a_bar * b_bar);
        let o_term = q_o.mul(c_bar);

        m_term + r_term + l_term + o_term + q_c.clone() + pi_eval
    }

    /// Constructs the permutation part of the quotient polynomial
    ///
    /// Enforces that the copy constraints hold via the grand product argument.
    /// Form: α · (Z(X) · Π(lhs) - Z(ωX) · Π(rhs)) / Z_H(X)
    ///
    /// lhs = (a(X)+β·X+γ)(b(X)+β·k₁X+γ)(c(X)+β·k₂X+γ)  
    /// rhs = (a(X)+β·σ₁(X)+γ)(b(X)+β·σ₂(X)+γ)(c(X)+β·σ₃(X)+γ)
    ///
    /// The result must be divisible by the vanishing polynomial `Z_H(X)`
    fn compute_permutation_summand(
        &self,
        a: &DensePolynomial<E::ScalarField>,
        b: &DensePolynomial<E::ScalarField>,
        c: &DensePolynomial<E::ScalarField>,
        sigma_polynomials: &(
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
        ),
        beta: E::ScalarField,
        gamma: E::ScalarField,
        k1: E::ScalarField,
        k2: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
        vanishing_poly: &DensePolynomial<E::ScalarField>,
        alpha: E::ScalarField,
    ) -> DensePolynomial<E::ScalarField> {
        let permutation_summand_1 = self.compute_first_permutation_summand(
            a.clone(),
            b.clone(),
            c.clone(),
            beta,
            gamma,
            k1,
            k2,
            &z,
        );

        let permutation_summand_2 = self.compute_second_permutation_summand(
            a.clone(),
            b.clone(),
            c.clone(),
            &sigma_polynomials,
            beta,
            gamma,
            &z,
        );

        let permutation_summand = permutation_summand_1.sub(&permutation_summand_2).mul(alpha);

        let (permutation_summand, r) = DenseOrSparsePolynomial::from(permutation_summand)
            .divide_with_q_and_r(&vanishing_poly.into())
            .unwrap();
        assert!(
            r.is_zero(),
            "permutation summand remainder is nonzero: {:?}",
            r
        );

        permutation_summand
    }

    /// Constructs the left-hand side of the permutation argument:
    /// Z(X) · Π (wire_i(X) + β·s + γ)
    ///
    /// where s ∈ {X, k₁·X, k₂·X}
    fn compute_first_permutation_summand(
        &self,
        a: DensePolynomial<E::ScalarField>,
        b: DensePolynomial<E::ScalarField>,
        c: DensePolynomial<E::ScalarField>,
        beta: E::ScalarField,
        gamma: E::ScalarField,
        k1: E::ScalarField,
        k2: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let a_multiply = vec_to_poly(vec![gamma, beta]) + a;
        let b_multiply = vec_to_poly(vec![gamma, beta * k1]) + b;
        let c_multiply = vec_to_poly(vec![gamma, beta * k2]) + c;

        a_multiply
            .naive_mul(&b_multiply)
            .naive_mul(&c_multiply)
            .naive_mul(&z)
    }

    /// Constructs the right-hand side of the permutation argument:
    /// Z(ω·X) · Π (wire_i(X) + β·σ_i(X) + γ)
    ///
    /// Uses σ₁(X), σ₂(X), σ₃(X) and evaluates Z at ω·X
    fn compute_second_permutation_summand(
        &self,
        a: DensePolynomial<E::ScalarField>,
        b: DensePolynomial<E::ScalarField>,
        c: DensePolynomial<E::ScalarField>,
        sigma_polynomials: &(
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
        ),
        beta: E::ScalarField,
        gamma: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let (sigma_1, sigma_2, sigma_3) = sigma_polynomials;

        let gamma_constant_poly = vec_to_poly(vec![gamma]);
        let a_multiply = a + sigma_1.mul(beta) + gamma_constant_poly.clone();
        let b_multiply = b + sigma_2.mul(beta) + gamma_constant_poly.clone();
        let c_multiply = c + sigma_3.mul(beta) + gamma_constant_poly;

        let mut shifted_z = z.clone();
        let omega = self.domain[1];
        for i in 1..shifted_z.len() {
            shifted_z[i] *= omega.pow(&[i as u64]);
        }

        a_multiply
            .naive_mul(&b_multiply)
            .naive_mul(&c_multiply)
            .naive_mul(&shifted_z)
    }

    /// Enforces Z(1) = 1 by constraining the Lagrange base polynomial L₁
    ///
    /// Computes α² · (Z(X) - 1) · L₁(X) / Z_H(X)
    ///
    /// Ensures Z starts correctly in the grand product construction
    fn compute_init_z_summand(
        &self,
        alpha: E::ScalarField,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
        z: &DensePolynomial<E::ScalarField>,
        lagrange_basis_1: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let last_summand = z
            .add(&constant_polynomial(-E::ScalarField::one()))
            .naive_mul(lagrange_basis_1)
            .mul(alpha.square());

        let (last_summand, r) = DenseOrSparsePolynomial::from(last_summand)
            .divide_with_q_and_r(&vanishing_polynomial.into())
            .unwrap();
        assert!(r.is_zero(), "last summand 1 remainder is nonzero: {:?}", r);

        last_summand
    }

    /// Split t(X) into t_lo + X^n · t_mid + X^{2n} · t_hi
    /// to allow commitment using a trusted setup with degree < 3n.
    fn split_quotient_polynomial(
        t: &DensePolynomial<E::ScalarField>,
        b10: E::ScalarField,
        b11: E::ScalarField,
        n: usize,
    ) -> (
        DensePolynomial<E::ScalarField>,
        DensePolynomial<E::ScalarField>,
        DensePolynomial<E::ScalarField>,
    ) {
        let mut t_lo_prime = Vec::new();
        let mut t_mid_prime = Vec::new();
        let mut t_hi_prime = Vec::new();

        for (i, c) in t.coeffs().iter().enumerate() {
            if i < n {
                t_lo_prime.push(c.clone());
            } else if i >= n && i < 2 * n {
                t_mid_prime.push(c.clone());
            } else if i >= 2 * n {
                t_hi_prime.push(c.clone());
            }
        }

        let mut t_lo = vec![E::ScalarField::zero(); n + 1];
        t_lo[n] = b10;
        let t_lo = vec_to_poly(t_lo).add(vec_to_poly(t_lo_prime));

        let mut t_mid = vec![E::ScalarField::zero(); n + 1];
        t_mid[0] = -b10;
        t_mid[n] = b11;
        let t_mid = vec_to_poly(t_mid).add(vec_to_poly(t_mid_prime));

        let t_hi = vec![-b11];
        let t_hi = vec_to_poly(t_hi).add(vec_to_poly(t_hi_prime));

        (t_lo, t_mid, t_hi)
    }

    /// Constructs the permutation polynomial Z(X) = blinding(X) + rolling_product(X)
    ///
    /// Blinding hides the structure and ensures zero-knowledge.  
    /// Rolling product encodes the grand product argument.
    fn compute_permutation_polynomial(
        b7: E::ScalarField,
        b8: E::ScalarField,
        b9: E::ScalarField,
        vanishing_poly: &DensePolynomial<E::ScalarField>,
        rolling_product: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let blinding_poly = vec_to_poly(vec![b9, b8, b7]);
        let blinding_poly = blinding_poly.naive_mul(vanishing_poly);
        blinding_poly + rolling_product.clone()
    }

    /// Blinds a witness polynomial and returns its coefficient form
    ///
    /// Form: a(X) = witness(X) + b₁·X + b₂·X² multiplied by the vanishing polynomial
    ///
    /// This is to ensure hiding and independence of wire polynomials in KZG
    fn compute_wire_coefficients_form(
        b1: E::ScalarField,
        b2: E::ScalarField,
        witness: &DensePolynomial<E::ScalarField>,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let blinding_poly = vec_to_poly(vec![b2, b1]).naive_mul(&vanishing_polynomial);
        blinding_poly + witness.clone()
    }

    /// Commit to a polynomial using CRS powers of τ in G1.
    /// This performs a KZG commitment.
    pub fn commit_polynomial(
        polynomial: &DensePolynomial<E::ScalarField>,
        crs: &[E::G1Affine],
        g1: E::G1Affine, // this is [1]_1
    ) -> E::G1Affine {
        if polynomial.coeffs.len() < 1 {
            return E::G1Affine::zero();
        }
        let mut acc = g1.into_group() * polynomial.coeffs[0]; // constant term

        for (i, coeff) in polynomial.coeffs.iter().skip(1).enumerate() {
            acc += crs[i + 1].into_group() * coeff;
        }

        acc.into_affine()
    }
}

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::Witness;
    use crate::gate::Gate;
    use crate::permutation::Permutation;
    use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Projective};
    use ark_ec::Group;
    use ark_ff::{FftField, Field, One};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_std::UniformRand;
    use ark_std::test_rng;

    fn fr(n: u64) -> Fr {
        Fr::from(n)
    }

    fn dummy_crs<E: Pairing>(degree: usize) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
        let mut rng = test_rng();
        let tau = E::ScalarField::rand(&mut rng);
        let g1_gen = E::G1::generator();
        let g2_gen = E::G2::generator();

        let crs_g1 = (0..=degree)
            .map(|i| (g1_gen * tau.pow(&[i as u64])).into_affine())
            .collect();

        let crs_g2 = (0..=degree)
            .map(|i| (g2_gen * tau.pow(&[i as u64])).into_affine())
            .collect();

        (crs_g1, crs_g2)
    }

    fn dummy_crs_tau<E: Pairing>(
        degree: usize,
    ) -> (Vec<E::G1Affine>, Vec<E::G2Affine>, E::ScalarField) {
        let mut rng = test_rng();
        let tau = E::ScalarField::rand(&mut rng);
        let g1_gen = E::G1::generator();
        let g2_gen = E::G2::generator();

        let crs_g1 = (0..=degree)
            .map(|i| (g1_gen * tau.pow(&[i as u64])).into_affine())
            .collect();

        let crs_g2 = (0..=degree)
            .map(|i| (g2_gen * tau.pow(&[i as u64])).into_affine())
            .collect();

        (crs_g1, crs_g2, tau)
    }

    #[test]
    fn test_compute_wire_poly_adds_blinding() {
        let domain = vec![fr(1), fr(2), fr(3)];
        let zh = Circuit::vanishing_poly(&domain);

        let witness_poly = DensePolynomial::from_coefficients_vec(vec![fr(5), fr(0), fr(0)]);
        let blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            fr(2),
            fr(3),
            &witness_poly,
            &zh,
        );

        // blinded(X) = witness + (3 + 2X)*Zh(X)
        let blind_term = DensePolynomial::from_coefficients_vec(vec![fr(3), fr(2)]).naive_mul(&zh);
        let expected = blind_term + witness_poly;

        assert_eq!(blinded, expected);
    }

    #[test]
    fn test_commit_polynomial_linear_combination() {
        let mut rng = test_rng();
        let (crs1, crs2) = dummy_crs::<Bls12_381>(3);
        let g1 = G1Projective::rand(&mut rng).into_affine();

        let poly1 = DensePolynomial::from_coefficients_vec(vec![fr(1), fr(2)]);
        let poly2 = DensePolynomial::from_coefficients_vec(vec![fr(3), fr(4)]);
        let sum = &poly1 + &poly2;

        let c1 = KZGProver::<Bls12_381>::commit_polynomial(&poly1, &crs1, g1);
        let c2 = KZGProver::<Bls12_381>::commit_polynomial(&poly2, &crs1, g1);
        let c_sum = KZGProver::<Bls12_381>::commit_polynomial(&sum, &crs1, g1);

        let c1_proj = c1.into_group();
        let c2_proj = c2.into_group();
        let c_sum_proj = c_sum.into_group();

        assert_eq!(c1_proj + c2_proj, c_sum_proj);
    }

    #[test]
    fn test_permutation_argument_rolling_product_validity() {
        let n = 4;
        let omega = Fr::get_root_of_unity(n as u64).unwrap();
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow([i as u64])).collect();
        let k1 = fr(5);
        let k2 = fr(7);
        let mut rng = test_rng();
        let beta = Fr::rand(&mut rng);
        let gamma = Fr::rand(&mut rng);

        // Simple wire permutation — swap a few values
        let wiring = vec![
            vec![0, 3], // a_0 == a_1 => σ(0) = 4
            vec![4, 7],
            vec![8, 11],
        ];

        let witness = Witness {
            a: vec![Fr::from(3); n],
            b: vec![Fr::from(4); n],
            c: vec![Fr::from(5); n],
        };

        let permutation = Permutation::new(witness.clone(), wiring.clone(), k1, k2);

        let sigma_maps = permutation.get_sigma_maps();
        let (sigma_a, sigma_b, sigma_c) =
            permutation.generate_sigma_polynomials(sigma_maps, &domain);

        let rolling_product_poly = permutation.get_rolling_product(gamma, beta, &domain);
        let z_evals = crate::fft::fft(&rolling_product_poly.coeffs, domain[1]);

        assert_eq!(z_evals[0], Fr::one());

        // 2. Recurrence holds for 1..n-1
        let mut expected = vec![Fr::one()];
        for i in 0..n - 1 {
            let x = domain[i];
            let a = witness.a[i];
            let b = witness.b[i];
            let c = witness.c[i];

            let num =
                (a + beta * x + gamma) * (b + beta * k1 * x + gamma) * (c + beta * k2 * x + gamma);

            let a_s = sigma_a.evaluate(&x);
            let b_s = sigma_b.evaluate(&x);
            let c_s = sigma_c.evaluate(&x);

            let den =
                (a + beta * a_s + gamma) * (b + beta * b_s + gamma) * (c + beta * c_s + gamma);

            expected.push(expected[i] * num / den);
            let actual = z_evals[i + 1];

            assert_eq!(expected[i + 1], actual, "Mismatch at i = {}", i);
        }
    }

    #[test]
    fn test_split_and_reconstruct_quotient_poly() {
        use ark_std::UniformRand;

        let mut rng = test_rng();
        let n = 8;

        // Construct t(X) with deg < 3n
        let t_coeffs: Vec<Fr> = (0..(3 * n - 1)).map(|_| Fr::rand(&mut rng)).collect();
        let t = DensePolynomial::from_coefficients_vec(t_coeffs.clone());

        let b10 = Fr::rand(&mut rng);
        let b11 = Fr::rand(&mut rng);

        let (t_lo, t_mid, t_hi) =
            KZGProver::<Bls12_381>::split_quotient_polynomial(&t, b10, b11, n);

        // Reconstruct t(X) = t_lo(X) + X^n * t_mid(X) + X^{2n} * t_hi(X)
        let mut xn = vec![Fr::zero(); n];
        xn.push(Fr::one());
        let xn_poly = DensePolynomial::from_coefficients_vec(xn.clone());

        let mut x2n = vec![Fr::zero(); 2 * n];
        x2n.push(Fr::one());
        let x2n_poly = DensePolynomial::from_coefficients_vec(x2n);

        let t_reconstructed = &t_lo + &(xn_poly.mul(&t_mid)) + (x2n_poly.mul(&t_hi));

        assert_eq!(
            t, t_reconstructed,
            "Reconstructed t(X) does not match original"
        );
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

    fn dummy_circuit2() -> Circuit<Fr> {
        let n = 4;
        let omega = Fr::get_root_of_unity(n as u64).unwrap();
        let domain: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        let gate0 = Gate::new(fr(1), fr(0), fr(0), fr(0), fr(0)); // a = x₀
        let gate1 = Gate::new(fr(1), fr(0), fr(0), fr(0), fr(0)); // a = x₁
        let gate2 = Gate::new(fr(1), fr(1), fr(1), -fr(1), fr(0)); // a + b + ab = c
        let gate3 = Gate::new(fr(1), fr(1), fr(0), -fr(1), fr(0)); // padding (optional)

        let gates = vec![gate0, gate1, gate2, gate3];

        // Witness
        let witness = Witness {
            a: vec![fr(3), fr(4), fr(3), fr(3)], // a[0]=x₀, a[1]=x₁, a[2]=x₀ for gate2
            b: vec![fr(0), fr(0), fr(4), fr(4)], // b[2]=x₁ for gate2
            c: vec![fr(0), fr(0), fr(19), fr(7)], // c[2] = 3 + 4 + 3 * 4 = 19
        };

        let public_inputs = vec![fr(3), fr(4)];

        let circuit = Circuit::new(
            gates,
            witness,
            public_inputs,
            domain.clone(),
            vec![vec![1, 6, 7], vec![0, 2, 3]],
            fr(3),
            fr(5),
        );
        assert!(circuit.is_gate_constraint_polynomial_zero_over_h(
            &circuit.get_selector_polynomials(),
            &circuit.get_witness_polynomials(),
            &circuit.compute_public_input_polynomial()
        ));
        circuit
    }
    #[test]
    fn test_gate_constraint_summand_divides_zh() {
        let circuit = dummy_circuit(); // Build a simple known-valid circuit
        let zh = Circuit::vanishing_poly(&circuit.domain);

        // get_gate_constraint_polynomial should be Q(X)
        let q = circuit.get_gate_constraint_polynomial();

        let zh_poly = DensePolynomial::from_coefficients_vec(zh.coeffs().to_vec());
        let gcp = DenseOrSparsePolynomial::from(q.clone());

        let (q_divided, r) = gcp.divide_with_q_and_r(&zh_poly.clone().into()).unwrap();
        assert!(
            r.is_zero(),
            "Gate constraint polynomial does not divide Z_H(X)"
        );

        // Optional: recompute the remainder manually
        let recomposed = q_divided.mul(&zh_poly);
        assert_eq!(
            recomposed, q,
            "Gate constraint summand incorrect: Q(X) != q(X) * Z_H(X)"
        );
    }

    #[test]
    fn test_permutation_summand_correctness_against_raw_constraint() {
        let circuit = dummy_circuit();
        let zh = Circuit::vanishing_poly(&circuit.domain);
        let witness_polys = circuit.get_witness_polynomials();

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, false);

        let mut rng = test_rng();
        let beta = Fr::rand(&mut rng);
        let gamma = Fr::rand(&mut rng);
        let alpha = Fr::rand(&mut rng);

        let vanishing_poly = Circuit::vanishing_poly(&circuit.domain);

        let a = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            fr(2),
            fr(3),
            &witness_polys.a,
            &vanishing_poly,
        );
        let b = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            fr(2),
            fr(3),
            &witness_polys.b,
            &vanishing_poly,
        );
        let c = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            fr(2),
            fr(3),
            &witness_polys.c,
            &vanishing_poly,
        );
        let z = circuit
            .permutation
            .get_rolling_product(gamma, beta, &circuit.domain);

        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        // Compute raw LHS - RHS (not divided)
        let lhs = prover.compute_first_permutation_summand(
            a.clone(),
            b.clone(),
            c.clone(),
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
        );

        let rhs = prover.compute_second_permutation_summand(
            a.clone(),
            b.clone(),
            c.clone(),
            &sigma_polys.clone(),
            beta,
            gamma,
            &z,
        );

        let raw_constraint = lhs.sub(&rhs).mul(alpha);

        // Now compute your summand (already divided by Z_H)
        let quotient_summand = prover.compute_permutation_summand(
            &a,
            &b,
            &c,
            &sigma_polys,
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            &zh,
            alpha,
        );

        let recomposed = quotient_summand.mul(&zh);

        assert_eq!(
            raw_constraint, recomposed,
            "Computed quotient summand is incorrect: q_perm(X) * Z_H(X) != raw_constraint"
        );
    }

    #[test]
    fn test_init_z_summand_correctness_against_raw_constraint() {
        let circuit = dummy_circuit();
        let zh = Circuit::vanishing_poly(&circuit.domain);

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, false);

        let mut rng = test_rng();
        let alpha = Fr::rand(&mut rng);
        let beta = Fr::rand(&mut rng);
        let gamma = Fr::rand(&mut rng);

        let z = circuit
            .permutation
            .get_rolling_product(gamma, beta, &circuit.domain);

        let lagrange_basis_1 = compute_lagrange_base(1, &circuit.domain);

        // raw form: α² · L₁(X) · (Z(X) - 1)
        let one = vec_to_poly(vec![Fr::one()]);
        let raw_constraint = lagrange_basis_1
            .clone()
            .naive_mul(&z.clone().sub(&one))
            .mul(&vec_to_poly(vec![alpha.pow([2])]));

        // already-divided form
        let quotient_summand = prover.compute_init_z_summand(alpha, &zh, &z, &lagrange_basis_1);

        let recomposed = quotient_summand.mul(&zh);

        assert_eq!(
            raw_constraint, recomposed,
            "Last quotient summand incorrect: q_last(X) * Z_H(X) != α² · L₁(X) · (Z(X) - 1)"
        );
    }

    #[test]
    fn test_constraint_linearisation_summand_correctness() {
        let circuit = dummy_circuit();
        let selector = circuit.get_selector_polynomials();
        let pi = circuit.compute_public_input_polynomial();
        let zeta = circuit.domain[2];

        let witness_polys = circuit.get_witness_polynomials();
        let a_bar = witness_polys.a.evaluate(&zeta);
        let b_bar = witness_polys.b.evaluate(&zeta);
        let c_bar = witness_polys.c.evaluate(&zeta);

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, false);
        let poly = prover
            .compute_constraint_linearisation_summand(a_bar, b_bar, c_bar, &selector, &pi, zeta);
        let eval = poly.evaluate(&zeta);

        // Manual check
        let expected = a_bar * selector.q_l.evaluate(&zeta)
            + b_bar * selector.q_r.evaluate(&zeta)
            + a_bar * b_bar * selector.q_m.evaluate(&zeta)
            + c_bar * selector.q_o.evaluate(&zeta)
            + selector.q_c.evaluate(&zeta)
            + pi.evaluate(&zeta);

        assert_eq!(eval, expected);
    }

    #[test]
    fn test_quotient_linearisation_summand_evaluation_matches_expected() {
        let circuit = dummy_circuit();
        let zh = Circuit::vanishing_poly(&circuit.domain);
        let zeta = circuit.domain[2];

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, false);

        // Use known t_lo, t_mid, t_hi
        let t_lo = DensePolynomial::from_coefficients_vec(vec![fr(1); 4]);
        let t_mid = DensePolynomial::from_coefficients_vec(vec![fr(2); 4]);
        let t_hi = DensePolynomial::from_coefficients_vec(vec![fr(3); 4]);

        let eval_t = t_lo.evaluate(&zeta)
            + zeta.pow([4]) * t_mid.evaluate(&zeta)
            + zeta.pow([8]) * t_hi.evaluate(&zeta);
        let zh_eval = zh.evaluate(&zeta);
        let expected = eval_t * zh_eval;

        let poly =
            prover.compute_quotient_linearization_summand(4, zeta, &zh, &t_lo, &t_mid, &t_hi);
        let actual = poly.evaluate(&zeta);

        assert_eq!(actual, expected, "quotient lin summand mismatch");
    }

    #[test]
    fn test_init_z_linearisation_summand_matches_expected() {
        let circuit = dummy_circuit();
        let zeta = circuit.domain[1];
        let lagrange_base = compute_lagrange_base(1, &circuit.domain);

        let mut rng = test_rng();
        let alpha = Fr::rand(&mut rng);
        let z = circuit.permutation.get_rolling_product(
            Fr::rand(&mut rng),
            Fr::rand(&mut rng),
            &circuit.domain,
        );

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, false);

        let lin = prover.compute_init_z_linearization_summand(alpha, zeta, &z, &lagrange_base);
        let eval = lin.evaluate(&zeta);

        let expected = (z.evaluate(&zeta) - fr(1)) * lagrange_base.evaluate(&zeta) * alpha.pow([2]);
        assert_eq!(eval, expected);
    }
    #[test]
    fn test_permutation_linearisation_summand_correctness() {
        let circuit = dummy_circuit();

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, false);
        let witness = circuit.get_witness_polynomials();
        let zeta = circuit.domain[1];

        let mut rng = test_rng();
        let alpha = Fr::rand(&mut rng);
        let beta = Fr::rand(&mut rng);
        let gamma = Fr::rand(&mut rng);

        let z = circuit
            .permutation
            .get_rolling_product(gamma, beta, &circuit.domain);
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        let a_bar = witness.a.evaluate(&zeta);
        let b_bar = witness.b.evaluate(&zeta);
        let c_bar = witness.c.evaluate(&zeta);
        let sigma_bar_1 = sigma_polys.0.evaluate(&zeta);
        let sigma_bar_2 = sigma_polys.1.evaluate(&zeta);
        let z_omega_bar = z.evaluate(&(zeta * circuit.domain[1]));

        let lin = prover.compute_permutation_linearization_summand(
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
            alpha,
            beta,
            gamma,
            zeta,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            &sigma_polys.2,
        );

        let lhs = z.evaluate(&zeta)
            * (a_bar + beta * zeta + gamma)
            * (b_bar + beta * circuit.permutation.k1 * zeta + gamma)
            * (c_bar + beta * circuit.permutation.k2 * zeta + gamma);

        let rhs = z_omega_bar
            * (a_bar + beta * sigma_bar_1 + gamma)
            * (b_bar + beta * sigma_bar_2 + gamma)
            * (c_bar + beta * sigma_polys.2.evaluate(&zeta) + gamma);

        let expected = (lhs - rhs) * alpha;
        let actual = lin.evaluate(&zeta);
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_quotient_summands_over_linearisations() {
        let circuit = dummy_circuit2();

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();
        let mut rng = test_rng();

        let mut prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, true);
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();
        let witness = circuit.get_witness_polynomials();

        prover.generate_proof(circuit.clone(), &blinding_scalars);
        let debug_info = prover.prover_debug_info.clone().unwrap();
        let alpha = debug_info.alpha;
        let beta = debug_info.beta;
        let gamma = debug_info.gamma;
        let zeta = debug_info.zeta;

        let vanishing_poly = &Circuit::vanishing_poly(&circuit.domain);
        let rolling_product = circuit
            .permutation
            .get_rolling_product(gamma, beta, &circuit.domain);
        let z = KZGProver::<Bls12_381>::compute_permutation_polynomial(
            blinding_scalars[6],
            blinding_scalars[7],
            blinding_scalars[8],
            vanishing_poly,
            &rolling_product,
        );
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        let a_blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            blinding_scalars[0],
            blinding_scalars[1],
            &witness.a,
            vanishing_poly,
        );
        let b_blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            blinding_scalars[2],
            blinding_scalars[3],
            &witness.b,
            vanishing_poly,
        );
        let c_blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            blinding_scalars[4],
            blinding_scalars[5],
            &witness.c,
            vanishing_poly,
        );
        let a_bar = a_blinded.evaluate(&zeta);
        let b_bar = b_blinded.evaluate(&zeta);
        let c_bar = c_blinded.evaluate(&zeta);
        let sigma_bar_1 = sigma_polys.0.evaluate(&zeta);
        let sigma_bar_2 = sigma_polys.1.evaluate(&zeta);
        let z_omega_bar = z.evaluate(&(zeta * circuit.domain[1]));
        let lagrange_basis = &compute_lagrange_base(1, &circuit.domain);
        let pi = circuit.compute_public_input_polynomial();
        let selector_polys = circuit.get_selector_polynomials();

        let constraint_summand = prover.compute_constraint_summand(
            &a_blinded,
            &b_blinded,
            &c_blinded,
            pi.clone(),
            vanishing_poly,
            &selector_polys,
        );

        let lin_constraint_summand = prover.compute_constraint_linearisation_summand(
            a_bar,
            b_bar,
            c_bar,
            &selector_polys,
            &pi,
            zeta,
        );
        assert_eq!(
            constraint_summand
                .evaluate(&zeta)
                .mul(vanishing_poly.evaluate(&zeta)),
            lin_constraint_summand.evaluate(&zeta),
            "Constraint portion differs at zeta"
        );

        let permutation_summand = prover.compute_permutation_summand(
            &a_blinded,
            &b_blinded,
            &c_blinded,
            &sigma_polys,
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            vanishing_poly,
            alpha,
        );

        let lin_permut_summand = prover.compute_permutation_linearization_summand(
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
            alpha,
            beta,
            gamma,
            zeta,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            &sigma_polys.2,
        );
        assert_eq!(
            permutation_summand
                .evaluate(&zeta)
                .mul(vanishing_poly.evaluate(&zeta)),
            lin_permut_summand.evaluate(&zeta)
        );

        let init_zsummand =
            prover.compute_init_z_summand(alpha, &vanishing_poly, &z, &lagrange_basis);
        let init_z_summand_lin =
            prover.compute_init_z_linearization_summand(alpha, zeta, &z, lagrange_basis);
        assert_eq!(
            init_zsummand
                .evaluate(&zeta)
                .mul(vanishing_poly.evaluate(&zeta)),
            init_z_summand_lin.evaluate(&zeta)
        );

        let t = constraint_summand + permutation_summand + init_zsummand;

        let (t_lo, t_mid, t_hi) = KZGProver::<Bls12_381>::split_quotient_polynomial(
            &t,
            blinding_scalars[9],
            blinding_scalars[10],
            circuit.domain.len(),
        );

        let quotient_summand = prover.compute_quotient_linearization_summand(
            circuit.domain.len(),
            zeta,
            vanishing_poly,
            &t_lo,
            &t_mid,
            &t_hi,
        );
        let zeta_n = zeta.pow(&[circuit.domain.len() as u64]);
        let zeta_2n = zeta_n.square();
        let t_reconstructed_z = t_lo.clone() + t_mid.mul(zeta_n) + t_hi.mul(zeta_2n);

        assert_eq!(
            quotient_summand.evaluate(&zeta),
            t_reconstructed_z
                .evaluate(&zeta)
                .mul(vanishing_poly.evaluate(&zeta))
        );

        let r = prover.compute_linearisation_polynomial(
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
            alpha,
            beta,
            gamma,
            zeta,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &selector_polys,
            &pi,
            &z,
            &sigma_polys.2,
            lagrange_basis,
            vanishing_poly,
            &t_lo,
            &t_mid,
            &t_hi,
        );

        assert!(r.evaluate(&zeta).is_zero())
    }

    #[test]
    fn test_linearisation_polynomial() {
        let circuit = dummy_circuit2();

        let (crs_g1, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let g2 = G2Projective::generator().into_affine();
        let mut rng = test_rng();

        let mut prover = KZGProver::<Bls12_381>::new(crs_g1, circuit.domain.clone(), g1, true);
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();
        let witness = circuit.get_witness_polynomials();

        prover.generate_proof(circuit.clone(), &blinding_scalars);
        let debug_info = prover.prover_debug_info.clone().unwrap();
        let alpha = debug_info.alpha;
        let beta = debug_info.beta;
        let gamma = debug_info.gamma;
        let zeta = debug_info.zeta;

        let vanishing_poly = &Circuit::vanishing_poly(&circuit.domain);
        let rolling_product = circuit
            .permutation
            .get_rolling_product(gamma, beta, &circuit.domain);
        let z = KZGProver::<Bls12_381>::compute_permutation_polynomial(
            blinding_scalars[6],
            blinding_scalars[7],
            blinding_scalars[8],
            vanishing_poly,
            &rolling_product,
        );
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        let a_blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            blinding_scalars[0],
            blinding_scalars[1],
            &witness.a,
            vanishing_poly,
        );
        let b_blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            blinding_scalars[2],
            blinding_scalars[3],
            &witness.b,
            vanishing_poly,
        );
        let c_blinded = KZGProver::<Bls12_381>::compute_wire_coefficients_form(
            blinding_scalars[4],
            blinding_scalars[5],
            &witness.c,
            vanishing_poly,
        );
        let a_bar = a_blinded.evaluate(&zeta);
        let b_bar = b_blinded.evaluate(&zeta);
        let c_bar = c_blinded.evaluate(&zeta);
        let sigma_bar_1 = sigma_polys.0.evaluate(&zeta);
        let sigma_bar_2 = sigma_polys.1.evaluate(&zeta);
        let z_omega_bar = z.evaluate(&(zeta * circuit.domain[1]));
        let lagrange_basis = &compute_lagrange_base(1, &circuit.domain);
        let pi = circuit.compute_public_input_polynomial();
        let selector_polys = circuit.get_selector_polynomials();

        let constraint_summand = prover.compute_constraint_summand(
            &a_blinded,
            &b_blinded,
            &c_blinded,
            pi.clone(),
            vanishing_poly,
            &selector_polys,
        );
        let permutation_summand = prover.compute_permutation_summand(
            &a_blinded,
            &b_blinded,
            &c_blinded,
            &sigma_polys,
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            vanishing_poly,
            alpha,
        );
        let init_zsummand =
            prover.compute_init_z_summand(alpha, &vanishing_poly, &z, &lagrange_basis);
        let t = constraint_summand + permutation_summand + init_zsummand;
        let (t_lo, t_mid, t_hi) = KZGProver::<Bls12_381>::split_quotient_polynomial(
            &t,
            blinding_scalars[9],
            blinding_scalars[10],
            circuit.domain.len(),
        );

        let r = prover.compute_linearisation_polynomial(
            a_bar,
            b_bar,
            c_bar,
            sigma_bar_1,
            sigma_bar_2,
            z_omega_bar,
            alpha,
            beta,
            gamma,
            zeta,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &selector_polys,
            &pi,
            &z,
            &sigma_polys.2,
            lagrange_basis,
            vanishing_poly,
            &t_lo,
            &t_mid,
            &t_hi,
        );
        assert_eq!(debug_info.linearisation_poly, r);
        assert!(r.evaluate(&zeta).is_zero(), "R is not zero at z");
        assert!(debug_info.linearisation_poly.evaluate(&zeta).is_zero());
    }

    #[test]
    fn test_opening_polynomial_reconstruction() {
        let circuit = dummy_circuit();
        let domain = circuit.domain.clone();
        let (crs, _) = dummy_crs::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();

        let mut prover = KZGProver::<Bls12_381>::new(crs, domain.clone(), g1, true);

        let mut rng = test_rng();
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();

        let proof = prover.generate_proof(circuit.clone(), &blinding_scalars);
        let debug_info = prover.prover_debug_info.as_ref().unwrap();

        let zeta = debug_info.zeta;
        let v = debug_info.v;

        // First, let's check if r(ζ) = 0
        let r_at_zeta = debug_info.linearisation_poly.evaluate(&zeta);
        println!("r(ζ) = {}", r_at_zeta);

        // Check individual terms at ζ
        let a_term_at_zeta = (debug_info.a_poly_blinded.evaluate(&zeta) - proof.a_bar) * v;
        let b_term_at_zeta = (debug_info.b_poly_blinded.evaluate(&zeta) - proof.b_bar) * v.pow([2]);
        let c_term_at_zeta = (debug_info.c_poly_blinded.evaluate(&zeta) - proof.c_bar) * v.pow([3]);

        println!("a_term(ζ) = {}", a_term_at_zeta);
        println!("b_term(ζ) = {}", b_term_at_zeta);
        println!("c_term(ζ) = {}", c_term_at_zeta);

        // These should all be zero!
        assert!(
            a_term_at_zeta.is_zero(),
            "a_poly_blinded(ζ) should equal a_bar"
        );
        assert!(
            b_term_at_zeta.is_zero(),
            "b_poly_blinded(ζ) should equal b_bar"
        );
        assert!(
            c_term_at_zeta.is_zero(),
            "c_poly_blinded(ζ) should equal c_bar"
        );

        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &domain);

        let sigma1_term_at_zeta = (sigma_polys.0.evaluate(&zeta) - proof.sigma_bar_1) * v.pow([4]);
        let sigma2_term_at_zeta = (sigma_polys.1.evaluate(&zeta) - proof.sigma_bar_2) * v.pow([5]);

        println!("sigma1_term(ζ) = {}", sigma1_term_at_zeta);
        println!("sigma2_term(ζ) = {}", sigma2_term_at_zeta);

        assert!(
            sigma1_term_at_zeta.is_zero(),
            "sigma1(ζ) should equal sigma_bar_1"
        );
        assert!(
            sigma2_term_at_zeta.is_zero(),
            "sigma2(ζ) should equal sigma_bar_2"
        );

        // Now check the total numerator
        let numerator = debug_info.linearisation_poly.clone()
            + debug_info
                .a_poly_blinded
                .clone()
                .sub(&constant_polynomial(proof.a_bar))
                .mul(v)
            + debug_info
                .b_poly_blinded
                .clone()
                .sub(&constant_polynomial(proof.b_bar))
                .mul(v.pow([2]))
            + debug_info
                .c_poly_blinded
                .clone()
                .sub(&constant_polynomial(proof.c_bar))
                .mul(v.pow([3]))
            + sigma_polys
                .0
                .clone()
                .sub(&constant_polynomial(proof.sigma_bar_1))
                .mul(v.pow([4]))
            + sigma_polys
                .1
                .clone()
                .sub(&constant_polynomial(proof.sigma_bar_2))
                .mul(v.pow([5]));

        let numerator_at_zeta = numerator.evaluate(&zeta);
        println!("Total numerator(ζ) = {}", numerator_at_zeta);

        assert!(numerator_at_zeta.is_zero(), "Numerator should be zero at ζ");

        // Only proceed with reconstruction if numerator is zero
        let denominator = vec_to_poly(vec![-zeta, Fr::one()]);
        let reconstructed_w_zeta = numerator.clone().div(&denominator);

        assert_eq!(
            debug_info.opening_poly, reconstructed_w_zeta,
            "Reconstructed opening polynomial doesn't match computed one"
        );
    }

    #[test]
    fn test_opening_polynomial_omega() {
        let circuit = dummy_circuit();
        let (crs_g1, _, x) = dummy_crs_tau::<Bls12_381>(circuit.domain.len() + 5);
        let g1 = G1Projective::generator().into_affine();
        let mut rng = test_rng();

        let mut prover =
            KZGProver::<Bls12_381>::new(crs_g1.clone(), circuit.domain.clone(), g1, true);
        let blinding_scalars: Vec<Fr> = (0..11).map(|_| Fr::rand(&mut rng)).collect();

        let proof = prover.generate_proof(circuit.clone(), &blinding_scalars);
        let debug_info = prover.prover_debug_info.as_ref().unwrap();

        let zeta = debug_info.zeta;
        let omega = circuit.domain[1];

        // For W_ζω, test similar properties
        let zeta_omega = zeta * omega;

        // Reconstruct Z polynomial
        let vanishing_poly = Circuit::vanishing_poly(&circuit.domain);
        let rolling_product = circuit.permutation.get_rolling_product(
            debug_info.gamma,
            debug_info.beta,
            &circuit.domain,
        );
        let z = KZGProver::<Bls12_381>::compute_permutation_polynomial(
            blinding_scalars[6],
            blinding_scalars[7],
            blinding_scalars[8],
            &vanishing_poly,
            &rolling_product,
        );

        let numerator_omega = z.clone().sub(&constant_polynomial(proof.z_omega_bar));
        assert!(
            numerator_omega.evaluate(&zeta_omega).is_zero(),
            "Omega numerator should be zero at ζω"
        );

        let denominator_omega = vec_to_poly(vec![-zeta_omega, Fr::one()]);
        let reconstructed_opening_omega = numerator_omega.clone().div(&denominator_omega);
        assert_eq!(
            debug_info.opening_omega_poly, reconstructed_opening_omega,
            "Opening omega polynomial division incorrect"
        );

        let opening_commitment =
            KZGProver::<Bls12_381>::commit_polynomial(&reconstructed_opening_omega, &crs_g1, g1);
        assert_eq!(opening_commitment, proof.w_zeta_omega);
        assert_eq!(
            g1 * reconstructed_opening_omega.evaluate(&x),
            proof.w_zeta_omega
        );
    }
}
