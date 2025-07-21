use crate::circuit::{Circuit, SelectorPolynomials, WitnessPolynomials};
use crate::fft::{compute_lagrange_base, vec_to_poly};
use crate::transccript::hash_to_field;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, One, Zero};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::{DenseUVPolynomial, Polynomial};
use std::ops::{Add, Mul, Sub};

pub struct KZGProver<E: Pairing> {
    crs: Vec<E::G1Affine>,
    domain: Vec<E::ScalarField>,
    g1: E::G1Affine,
    g2: E::G2Affine,
}

impl<E: Pairing> KZGProver<E> {
    pub fn generate_proof(
        &self,
        circuit: Circuit<E::ScalarField>,
        blinding_scalars: &[E::ScalarField],
    ) {
        let n = circuit.gates.len();
        assert_eq!(blinding_scalars.len(), 11);
        assert_eq!(self.crs.len(), n + 5);

        let mut commitment_buffer = Vec::new();
        let vanishing_poly = Circuit::vanishing_poly(&self.domain);
        let witness_polynomial = circuit.get_witness_polynomials();
        let a = witness_polynomial.a.clone();
        let b = witness_polynomial.b.clone();
        let c = witness_polynomial.c.clone();

        // round one, computing wire polynomials
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

        // round two, computing permutation poly
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

        // round three, computing quotient polynomial
        // first summand
        let gcp = DenseOrSparsePolynomial::from(circuit.get_gate_constraint_polynomial());
        let zh = DenseOrSparsePolynomial::from(vanishing_poly.clone());
        let (first_summand, r) = gcp.divide_with_q_and_r(&zh).unwrap();
        assert!(
            r.is_zero(),
            "gate constraint has non zero remainder: {:?}",
            r
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
            witness_polynomial,
            &sigma_polynomials,
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            &vanishing_poly,
            alpha,
        );

        let quotient_polynomial = first_summand + permutation_summand + init_z_summand;

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

        // round 4
        let zeta = hash_to_field("zeta", &commitment_buffer);

        let a_bar = a.evaluate(&zeta);
        let b_bar = b.evaluate(&zeta);
        let c_bar = c.evaluate(&zeta);

        let sigma_bar_1 = sigma_polynomials.0.evaluate(&zeta);
        let sigma_bar_2 = sigma_polynomials.1.evaluate(&zeta);

        let z_omega_bar = z.evaluate(&zeta.mul(self.domain[1]));

        // round 5
        // previous round's output could be added to transcript
        // however, I will just add v for now.
        // let v = hash_to_field("v", &commitment_buffer);
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
            &circuit.get_selector_polynomials(),
            &circuit.compute_public_input_polynomial(),
            &z,
            &sigma_polynomials.2,
            &lagrange_base_1,
            &vanishing_poly,
            &t_lo,
            &t_mid,
            &t_hi,
        );
    }

    fn compute_linearisation_polynomial(
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
        let constraint_lin_summand = self.calculate_constraint_linearisation_summand(
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

        constraint_lin_summand + permutation_lin_summand + init_z_lin_summand + quotient_summand
    }

    fn compute_quotient_linearization_summand(
        &self,
        n: usize,
        zeta: E::ScalarField,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
        t_lo: &DensePolynomial<E::ScalarField>,
        t_mid: &DensePolynomial<E::ScalarField>,
        t_hi: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let vanishing_evaluations = vanishing_polynomial.evaluate(&zeta);

        let zeta_n = zeta.pow(&[n as u64]);
        let zeta_2n = zeta_n.square();
        let mut t_mid_shifted = vec![E::ScalarField::zero(); n];
        let mut t_hi_shifted = vec![E::ScalarField::zero(); n * 2];
        t_mid.iter().for_each(|&x| t_mid_shifted.push(x * zeta_n));
        t_hi.iter().for_each(|&x| t_hi_shifted.push(x * zeta_2n));

        (t_lo.clone() + vec_to_poly(t_mid_shifted) + vec_to_poly(t_hi_shifted))
            .mul(vanishing_evaluations)
    }
    fn compute_init_z_linearization_summand(
        &self,
        alpha: E::ScalarField,
        zeta: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
        lagrange_base_1: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let lagrange_base_evaluation = lagrange_base_1.evaluate(&zeta);
        // Z(x) - 1
        z.sub(&vec_to_poly(vec![-E::ScalarField::one()]))
            // (Z(x) - 1)L1(z)
            .mul(lagrange_base_evaluation)
            // a^2 *[(Z(x) - 1)L1(z)]
            .mul(alpha.pow(&[2u64]))
    }

    fn compute_permutation_linearization_summand(
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
        let permut_2_c = sigma_3.mul(beta) + vec_to_poly(vec![c_bar + gamma]);
        let permut_summand_2 = permut_2_c.mul(
            (a_bar + beta * sigma_bar_1 + gamma)
                * (b_bar + beta * sigma_bar_2 + gamma)
                * z_omega_bar,
        );

        permut_summand_1.sub(&permut_summand_2).mul(alpha)
    }

    fn calculate_constraint_linearisation_summand(
        &self,
        a_bar: E::ScalarField,
        b_bar: E::ScalarField,
        c_bar: E::ScalarField,
        selector_polynomials: &SelectorPolynomials<E::ScalarField>,
        public_polynomial: &DensePolynomial<E::ScalarField>,
        zeta: E::ScalarField,
    ) -> DensePolynomial<E::ScalarField> {
        let pi_eval = vec_to_poly(vec![public_polynomial.evaluate(&zeta)]);
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

    fn compute_permutation_summand(
        &self,
        witness_polynomials: WitnessPolynomials<E::ScalarField>,
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
            witness_polynomials.clone(),
            beta,
            gamma,
            k1,
            k2,
            &z,
        );

        let permutation_summand_2 = self.compute_second_permutation_summand(
            witness_polynomials.clone(),
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

    fn compute_first_permutation_summand(
        &self,
        witness_polynomials: WitnessPolynomials<E::ScalarField>,
        beta: E::ScalarField,
        gamma: E::ScalarField,
        k1: E::ScalarField,
        k2: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let WitnessPolynomials { a, b, c } = witness_polynomials;
        let a_multiply = vec_to_poly(vec![gamma, beta]) + a;
        let b_multiply = vec_to_poly(vec![gamma, beta * k1]) + b;
        let c_multiply = vec_to_poly(vec![gamma, beta * k2]) + c;

        a_multiply
            .naive_mul(&b_multiply)
            .naive_mul(&c_multiply)
            .naive_mul(&z)
    }

    fn compute_second_permutation_summand(
        &self,
        witness_polynomials: WitnessPolynomials<E::ScalarField>,
        sigma_polynomials: &(
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
        ),
        beta: E::ScalarField,
        gamma: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let WitnessPolynomials { a, b, c } = witness_polynomials;
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

    fn compute_init_z_summand(
        &self,
        alpha: E::ScalarField,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
        z: &DensePolynomial<E::ScalarField>,
        lagrange_basis_1: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let last_summand = z
            .add(&vec_to_poly(vec![-E::ScalarField::one()]))
            .naive_mul(lagrange_basis_1)
            .mul(&vec_to_poly(vec![alpha.pow(&[2u64])]));

        let (last_summand, r) = DenseOrSparsePolynomial::from(last_summand)
            .divide_with_q_and_r(&vanishing_polynomial.into())
            .unwrap();
        assert!(r.is_zero(), "last summand 1 remainder is nonzero: {:?}", r);

        last_summand
    }

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
            }
            if i >= n && i < 2 * n {
                t_mid_prime.push(c.clone());
            }
            if i >= 2 * n {
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

    fn compute_wire_coefficients_form(
        b1: E::ScalarField,
        b2: E::ScalarField,
        witness: &DensePolynomial<E::ScalarField>,
        vanishing_polynomial: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let blinding_poly = vec_to_poly(vec![b2, b1]).naive_mul(&vanishing_polynomial);
        blinding_poly + witness.clone()
    }

    fn commit_polynomial(
        polynomial: &DensePolynomial<E::ScalarField>,
        crs: &[E::G1Affine],
        g1: E::G1Affine, // this is [1]_1
    ) -> E::G1Affine {
        let mut acc = g1.into_group() * polynomial.coeffs[0]; // constant term

        for (i, coeff) in polynomial.coeffs.iter().skip(1).enumerate() {
            acc += crs[i].into_group() * coeff; // crs[i] = [x^{i+1}]_1
        }

        acc.into_affine()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gate::Gate;
    use crate::permutation::Permutation;
    use crate::witness::Witness;
    use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
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
        use ark_ec::Group;
        use ark_ff::UniformRand;
        use ark_std::test_rng;

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

        // 1. Last value is 1
        println!("z evals: {:?}", z_evals);
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

        let prover = KZGProver::<Bls12_381> {
            crs: crs_g1,
            domain: circuit.domain.clone(),
            g1,
            g2,
        };

        let mut rng = test_rng();
        let beta = Fr::rand(&mut rng);
        let gamma = Fr::rand(&mut rng);
        let alpha = Fr::rand(&mut rng);

        let z = circuit
            .permutation
            .get_rolling_product(gamma, beta, &circuit.domain);

        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polys = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &circuit.domain);

        // Compute raw LHS - RHS (not divided)
        let lhs = prover.compute_first_permutation_summand(
            witness_polys.clone(),
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
        );

        let rhs = prover.compute_second_permutation_summand(
            witness_polys.clone(),
            &sigma_polys.clone(),
            beta,
            gamma,
            &z,
        );

        let raw_constraint = lhs.sub(&rhs).mul(alpha);

        // Now compute your summand (already divided by Z_H)
        let quotient_summand = prover.compute_permutation_summand(
            witness_polys,
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
        let g2 = G2Projective::generator().into_affine();

        let prover = KZGProver::<Bls12_381> {
            crs: crs_g1,
            domain: circuit.domain.clone(),
            g1,
            g2,
        };

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
}
