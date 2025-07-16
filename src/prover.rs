use crate::circuit::{Circuit, WitnessPolynomials};
use crate::fft::{compute_lagrange_base, vec_to_poly};
use crate::witness;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{Field, One, Zero};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use std::ops::{Add, Mul, Sub};

pub struct KZGProver<E: Pairing> {
    crs: Vec<E::G1Affine>,
    domain: Vec<E::ScalarField>,
    g1: E::G1Affine,
    g2: E::G1Affine,
}

impl<E: Pairing> KZGProver<E> {
    pub fn generate_proof(
        &self,
        circuit: Circuit<E::ScalarField>,
        blinding_scalars: &[E::ScalarField],
        beta: E::ScalarField,
        gamma: E::ScalarField,
        alpha: E::ScalarField,
    ) {
        let n = circuit.gates.len();
        assert_eq!(blinding_scalars.len(), 11);
        assert_eq!(self.crs.len(), n + 5);

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

        let b_poly = Self::compute_wire_coefficients_form(
            blinding_scalars[2],
            blinding_scalars[3],
            &b,
            &vanishing_poly,
        );
        let b_commitment = Self::commit_polynomial(&b_poly, &self.crs, self.g1);

        let c_poly = Self::compute_wire_coefficients_form(
            blinding_scalars[4],
            blinding_scalars[5],
            &c,
            &vanishing_poly,
        );
        let c_commitment = Self::commit_polynomial(&c_poly, &self.crs, self.g1);

        // round two, computing permutation poly
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

        let permutation_summand_1 = self.compute_first_permutation_summand(
            witness_polynomial.clone(),
            beta,
            gamma,
            circuit.permutation.k1,
            circuit.permutation.k2,
            &z,
            &vanishing_poly,
            alpha,
        );

        // second summand, permutation 2
        let sigma_maps = circuit.permutation.get_sigma_maps();
        let sigma_polynomials = circuit
            .permutation
            .generate_sigma_polynomials(sigma_maps, &self.domain);
        let permutation_summand_2 = self.compute_second_permutation_summand(
            witness_polynomial.clone(),
            sigma_polynomials,
            beta,
            gamma,
            &z,
            &vanishing_poly,
            alpha,
        );

        let lagrange_base_1 = compute_lagrange_base(1, &self.domain);
        let last_summand = self.compute_last_summand(alpha, &vanishing_poly, &z, &lagrange_base_1);

        let quotient_polynomial = first_summand
            .add(permutation_summand_1)
            .sub(&permutation_summand_2)
            + last_summand;

        let (t_lo, t_mid, t_hi) = Self::split_quotient_polynomial(
            &quotient_polynomial,
            blinding_scalars[9],
            blinding_scalars[10],
            self.domain.len(),
        );

        let t_lo_commitment = Self::commit_polynomial(&t_lo, &self.crs, self.g1);
        let t_mid_commitment = Self::commit_polynomial(&t_mid, &self.crs, self.g1);
        let t_hi_commitment = Self::commit_polynomial(&t_hi, &self.crs, self.g1);
    }

    fn compute_first_permutation_summand(
        &self,
        witness_polynomials: WitnessPolynomials<E::ScalarField>,
        beta: E::ScalarField,
        gamma: E::ScalarField,
        k1: E::ScalarField,
        k2: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
        vanishing_poly: &DensePolynomial<E::ScalarField>,
        alpha: E::ScalarField,
    ) -> DensePolynomial<E::ScalarField> {
        let WitnessPolynomials { a, b, c } = witness_polynomials;
        let a_multiply = vec_to_poly(vec![gamma, beta]) + a;
        let b_multiply = vec_to_poly(vec![gamma, beta * k1]) + b;
        let c_multiply = vec_to_poly(vec![gamma, beta * k2]) + c;

        let permutation_portion = a_multiply
            .naive_mul(&b_multiply)
            .naive_mul(&c_multiply)
            .naive_mul(&z)
            .mul(alpha);

        let (permutation_summand, r) = DenseOrSparsePolynomial::from(permutation_portion)
            .divide_with_q_and_r(&vanishing_poly.into())
            .unwrap();
        assert!(
            r.is_zero(),
            "permutation summand 1 remainder is nonzero: {:?}",
            r
        );

        permutation_summand
    }

    fn compute_second_permutation_summand(
        &self,
        witness_polynomials: WitnessPolynomials<E::ScalarField>,
        sigma_polynomials: (
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
            DensePolynomial<E::ScalarField>,
        ),
        beta: E::ScalarField,
        gamma: E::ScalarField,
        z: &DensePolynomial<E::ScalarField>,
        vanishing_poly: &DensePolynomial<E::ScalarField>,
        alpha: E::ScalarField,
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

        let permutation_portion = a_multiply
            .naive_mul(&b_multiply)
            .naive_mul(&c_multiply)
            .naive_mul(&shifted_z)
            .mul(alpha);

        let (permutation_summand, r) = DenseOrSparsePolynomial::from(permutation_portion)
            .divide_with_q_and_r(&vanishing_poly.into())
            .unwrap();
        assert!(
            r.is_zero(),
            "permutation summand 1 remainder is nonzero: {:?}",
            r
        );

        permutation_summand
    }

    fn compute_last_summand(
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
    use crate::permutation::Permutation;
    use crate::witness::Witness;
    use ark_bls12_381::{Bls12_381, Fr, G1Projective};
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
}
