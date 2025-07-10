use ark_ff::Field;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use itertools::Itertools;

/// Recursive Cooley–Tukey FFT over a finite field
///
/// Input:
/// - `coefficients`: vector of polynomial coefficients (in increasing order)
/// - `omega`: a primitive n-th root of unity in the field
///
/// Output:
/// - Vector of `f(ω^0), f(ω^1), ..., f(ω^{n-1})`
pub fn fft<F: Field>(coefficients: &[F], omega: F) -> Vec<F> {
    let n = coefficients.len();
    assert!(n.is_power_of_two(), "Input size must be power of two");
    if n == 1 {
        return coefficients.to_vec(); // base case
    }

    // Split coefficients into even and odd powers:
    // f(x) = f_even(x²) + x·f_odd(x²)
    let even = coefficients.iter().step_by(2).cloned().collect::<Vec<_>>();
    let odd = coefficients
        .iter()
        .skip(1)
        .step_by(2)
        .cloned()
        .collect::<Vec<_>>();

    // Recursively evaluate even and odd parts on ω² domain
    let even_eval = fft(&even, omega.square());
    let odd_eval = fft(&odd, omega.square());

    let mut r = vec![F::zero(); n];

    // Combine the even and odd evaluations:
    // f(ω^k)     = even(ω²^k) + ω^k · odd(ω²^k)
    // f(ω^{k+n/2}) = even(ω²^k) - ω^k · odd(ω²^k)
    let mut w = F::one();
    for i in 0..n / 2 {
        let t = w * odd_eval[i];
        r[i] = even_eval[i] + t;
        r[i + n / 2] = even_eval[i] - t;
        w *= omega;
    }

    r
}

pub fn inverse_fft<F: Field>(evaluations: &[F], omega: F) -> Vec<F> {
    let n = evaluations.len();
    let omega_inv = omega.inverse().unwrap();
    let mut result = fft(evaluations, omega_inv);
    let n_inv = F::from(n as u64).inverse().unwrap();
    for x in &mut result {
        *x *= n_inv;
    }
    result
}

/// Convert a coefficient vector into a DensePolynomial, trimming trailing zeros.
pub fn vec_to_poly<F: Field>(mut coeffs: Vec<F>) -> DensePolynomial<F> {
    for i in (0..coeffs.len()).rev() {
        if !coeffs[i].is_zero() {
            break;
        }
        coeffs.pop();
    }

    // Avoid creating zero-degree zero polynomial (arkworks panics)
    if coeffs.is_empty() {
        coeffs.push(F::zero());
    }

    DensePolynomial::from_coefficients_vec(coeffs)
}

pub fn pad_with_zeroes<F: Field>(coeffs: &[F], domain_size: usize) -> Vec<F> {
    coeffs
        .to_vec()
        .into_iter()
        .pad_using(domain_size, |_| F::zero())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_poly::EvaluationDomain;
    use ark_poly::domain::general::GeneralEvaluationDomain;

    /// Evaluate a polynomial at a given x value
    fn naive_eval<F: Field>(coeffs: &[F], x: F) -> F {
        let mut acc = F::zero();
        for &c in coeffs.iter().rev() {
            acc = acc * x + c;
        }
        acc
    }

    #[test]
    fn test_fft_roundtrip() {
        let domain = GeneralEvaluationDomain::<Fr>::new(8).unwrap();
        let omega = domain.group_gen();

        let original_coeffs = vec![
            Fr::from(3u64),
            Fr::from(2u64),
            Fr::from(0u64),
            Fr::from(1u64),
            Fr::from(4u64),
            Fr::from(0u64),
            Fr::from(0u64),
            Fr::from(0u64),
        ];

        let evals = fft(&original_coeffs, omega);
        let recovered = inverse_fft(&evals, omega);

        assert_eq!(original_coeffs, recovered);
    }

    #[test]
    fn test_fft_matches_naive_eval() {
        use ark_bls12_381::Fr;
        use ark_poly::EvaluationDomain;
        use ark_poly::domain::general::GeneralEvaluationDomain;

        let domain = GeneralEvaluationDomain::<Fr>::new(8).unwrap();
        let omega = domain.group_gen();

        let coeffs = vec![
            Fr::from(3u64),
            Fr::from(2u64),
            Fr::from(0u64),
            Fr::from(1u64),
            Fr::from(4u64),
            Fr::from(0u64),
            Fr::from(0u64),
            Fr::from(0u64),
        ];

        let fft_result = fft(&coeffs, omega);

        for i in 0..8 {
            let x = omega.pow([i as u64]);
            let expected = naive_eval(&coeffs, x);
            assert_eq!(fft_result[i], expected, "Mismatch at ω^{}", i);
        }
    }

    #[test]
    fn test_ifft_roundtrip() {
        use ark_bls12_381::Fr;
        use ark_poly::EvaluationDomain;
        use ark_poly::domain::general::GeneralEvaluationDomain;

        let domain = GeneralEvaluationDomain::<Fr>::new(8).unwrap();
        let omega = domain.group_gen();

        let coeffs = vec![
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(0u64),
            Fr::from(4u64),
            Fr::from(1u64),
            Fr::from(0u64),
            Fr::from(0u64),
            Fr::from(0u64),
        ];

        let evals = fft(&coeffs, omega);
        let result = inverse_fft(&evals, omega);

        assert_eq!(result, coeffs, "IFFT(FFT(f)) != f");
    }
}
