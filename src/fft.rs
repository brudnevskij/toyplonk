use ark_ff::Field;


/// Recursive Cooley–Tukey FFT over a finite field
///
/// Input:
/// - `coefficients`: vector of polynomial coefficients (in increasing order)
/// - `omega`: a primitive n-th root of unity in the field
///
/// Output:
/// - Vector of `f(ω^0), f(ω^1), ..., f(ω^{n-1})`
fn fft<F: Field>(coefficients: &[F], omega: F) -> Vec<F> {
    let n = coefficients.len();
    if n == 1 {
        return coefficients.to_vec(); // base case
    }

    // Split coefficients into even and odd powers:
    // f(x) = f_even(x²) + x·f_odd(x²)
    let even = coefficients.iter().step_by(2).cloned().collect::<Vec<_>>();
    let odd  = coefficients.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>();

    // Recursively evaluate even and odd parts on ω² domain
    let even_eval = fft(&even, omega.square());
    let odd_eval  = fft(&odd,  omega.square());

    let mut r = vec![F::zero(); n];

    // Combine the even and odd evaluations:
    // f(ω^k)     = even(ω²^k) + ω^k · odd(ω²^k)
    // f(ω^{k+n/2}) = even(ω²^k) - ω^k · odd(ω²^k)
    let mut w = F::one();
    for i in 0..n / 2 {
        let t = w * odd_eval[i];
        r[i]        = even_eval[i] + t;
        r[i + n/2]  = even_eval[i] - t;
        w *= omega;
    }

    r
}

