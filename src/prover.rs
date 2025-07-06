use crate::circuit::{Circuit, WitnessPolynomials};
use ark_ec::pairing::{Pairing};
use ark_ec::{AffineRepr, CurveGroup};
use ark_poly::univariate::DensePolynomial;
use crate::fft::vec_to_poly;

pub struct KZGProver<E: Pairing> {
    crs: Vec<E::G1Affine>,
    domain: Vec<E::ScalarField>,
    g1: E::G1Affine,
    g2: E::G1Affine,
}


impl<E: Pairing> KZGProver<E> {
    pub fn generate_proof(&self, circuit: Circuit<E::ScalarField>, blinding_scalars: &[E::ScalarField]) {
        let n = circuit.gates.len();
        assert_eq!(blinding_scalars.len(), 9);
        assert_eq!(self.crs.len(), n+5);

        let selector_polynomials = circuit.get_selector_polynomials();
        let vanishing_poly = Circuit::vanishing_poly(&self.domain);
        let WitnessPolynomials{ a, b, c }= circuit.get_witness_polynomials();

        // round one, computing wire polynomials
        let a_poly = Self::compute_wire_coefficients_form(blinding_scalars[0], blinding_scalars[1], &a, &vanishing_poly);
        let a_commitment = Self::commit_polynomial(&a_poly, &self.crs, self.g1);

        let b_poly = Self::compute_wire_coefficients_form(blinding_scalars[2], blinding_scalars[3], &b, &vanishing_poly);
        let b_commitment = Self::commit_polynomial(&b_poly, &self.crs, self.g1);

        let c_poly = Self::compute_wire_coefficients_form(blinding_scalars[4], blinding_scalars[5], &c, &vanishing_poly);
        let c_commitment = Self::commit_polynomial(&c_poly, &self.crs, self.g1);
    }

    fn compute_wire_coefficients_form(
        b1: E::ScalarField,
        b2: E::ScalarField,
        witness: &DensePolynomial<E::ScalarField>,
        zh_eval: &DensePolynomial<E::ScalarField>,
    ) -> DensePolynomial<E::ScalarField> {
        let blinding_poly = vec_to_poly(vec![b2, b1]).naive_mul(&zh_eval);
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
