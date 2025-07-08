use crate::fft::{fft, inverse_fft, vec_to_poly};
use crate::witness::Witness;
use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_std::iterable::Iterable;
use itertools::izip;

pub struct Permutation<F: Field> {
    pub witness: Witness<F>,
    pub wiring: Vec<Vec<usize>>,
    k1: F,
    k2: F,
}

impl<F: Field> Permutation<F> {
    pub fn new(witness: Witness<F>, wiring: Vec<Vec<usize>>, k1: F, k2: F) -> Permutation<F> {
        Self {
            witness,
            wiring,
            k1,
            k2,
        }
    }

    pub fn generate_sigma_mapping(&self) -> Vec<usize> {
        let n = self.witness.a.len();
        let mut sigma_map: Vec<_> = (0..3 * n).collect();

        for equal_wires in &self.wiring {
            if equal_wires.len() >= 2 {
                let mut rotated = equal_wires.clone();
                rotated.rotate_left(1);
                for (&from, &to) in equal_wires.iter().zip(rotated.iter()) {
                    sigma_map[from] = to;
                }
            }
        }

        sigma_map
    }

    pub fn get_sigma_maps(&self) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let sigma = self.generate_sigma_mapping();
        let n = self.witness.a.len();
        let a_sigma = sigma[..n].to_vec();
        let b_sigma = sigma[n..2 * n].to_vec();
        let c_sigma = sigma[2 * n..].to_vec();

        (a_sigma, b_sigma, c_sigma)
    }

    pub fn get_sigma_polynomials(
        &self,
        mappings: (Vec<usize>, Vec<usize>, Vec<usize>),
        domain: &[F],
    ) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        let omega = domain[1];

        let b_coset = domain.iter().map(|&x| self.k1 * x).collect::<Vec<_>>();
        let c_coset = domain.iter().map(|&x| self.k2 * x).collect::<Vec<_>>();

        let mut h_prime = Vec::from(domain.clone()); // h
        h_prime.extend(b_coset); // h ∪ k1·h
        h_prime.extend(c_coset); // ∪ k2·h

        let a_evaluations = mappings.0.iter().map(|&i| h_prime[i]).collect::<Vec<_>>();
        let b_evaluations = mappings.1.iter().map(|&i| h_prime[i]).collect::<Vec<_>>();
        let c_evaluations = mappings.2.iter().map(|&i| h_prime[i]).collect::<Vec<_>>();

        let a_sigma = vec_to_poly(inverse_fft(&a_evaluations, omega));
        let b_sigma = vec_to_poly(inverse_fft(&b_evaluations, omega));
        let c_sigma = vec_to_poly(inverse_fft(&c_evaluations, omega));

        (a_sigma, b_sigma, c_sigma)
    }

    // Calculates rolling product and outputs interpolation with
    // [Z0,..., Zn-1] values.
    pub fn calculate_rolling_product(
        &self,
        sigma_polys: (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>),
        domain: &[F],
        gamma: F,
        beta: F,
    ) -> DensePolynomial<F> {
        let omega = domain[1];
        let (a_sigma, b_sigma, c_sigma) = sigma_polys;

        let mut z_evaluations = vec![F::one()];
        let a = &self.witness.a;
        let b = &self.witness.b;
        let c = &self.witness.c;

        // calculating all Z(x) evals except last one.
        for i in 0..self.witness.a.len() {
            let a_num = a[i] + beta * domain[i] + gamma;
            let b_num = b[i] + beta * domain[i] * self.k1 + gamma;
            let c_num = c[i] + beta * domain[i] * self.k2 + gamma;
            let numerator = a_num * b_num * c_num;

            let a_denom = a[i] + beta * a_sigma.evaluate(&domain[i]) + gamma;
            let b_denom = b[i] + beta * b_sigma.evaluate(&domain[i]) + gamma;
            let c_denom = c[i] + beta * c_sigma.evaluate(&domain[i]) + gamma;
            let denominator = a_denom * b_denom * c_denom;

            z_evaluations.push(z_evaluations.last().unwrap().clone() * (numerator / denominator));
        }

        z_evaluations = z_evaluations[..z_evaluations.len() - 1].to_vec();

        vec_to_poly(inverse_fft(&z_evaluations, omega))
    }

    // convenience function producing rolling product
    pub fn get_rolling_product(&self, gamma: F, beta: F, domain: &[F]) -> DensePolynomial<F> {
        let sigma_maps = self.get_sigma_maps();
        let sigma_polys = self.get_sigma_polynomials(sigma_maps.clone(), domain);
        self.calculate_rolling_product(sigma_polys, &domain, gamma, beta)
    }
}

pub fn verify_permutation_argument<F: Field>(
    domain: &[F],
    z_evals: &[F],
    witness: &Witness<F>,
    sigma_maps: &(Vec<usize>, Vec<usize>, Vec<usize>),
    k1: F,
    k2: F,
    beta: F,
    gamma: F,
) -> Result<(), String> {
    let n = domain.len();
    let a = &witness.a;
    let b = &witness.b;
    let c = &witness.c;

    // H' calculation
    let b_coset = domain.iter().map(|&x| k1 * x).collect::<Vec<_>>();
    let c_coset = domain.iter().map(|&x| k2 * x).collect::<Vec<_>>();
    let mut h_prime = domain.to_vec();
    h_prime.extend(b_coset);
    h_prime.extend(c_coset);

    let mut expected = vec![F::one()];
    for i in 0..n {
        let x = domain[i];
        let a_i = a[i];
        let b_i = b[i];
        let c_i = c[i];

        let numerator = (a_i + beta * x + gamma)
            * (b_i + beta * k1 * x + gamma)
            * (c_i + beta * k2 * x + gamma);

        let a_sigma = sigma_maps.0[i];
        let b_sigma = sigma_maps.1[i];
        let c_sigma = sigma_maps.2[i];

        let a_s = h_prime[a_sigma];
        let b_s = h_prime[b_sigma];
        let c_s = h_prime[c_sigma];

        let denominator =
            (a_i + beta * a_s + gamma) * (b_i + beta * b_s + gamma) * (c_i + beta * c_s + gamma);

        let next = expected[i] * (numerator / denominator);
        expected.push(next);

        if z_evals[i] != expected[i] {
            return Err(format!(
                "Permutation recurrence mismatch at i = {}:\nexpected {}\nbut got {}",
                i, next, z_evals[i]
            ));
        }
    }

    if expected.last() != Some(&F::one()) {
        return Err(format!(
            "Last Z(w^n) != 1, Z(w^n) = {}",
            expected.last().unwrap()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::{FftField, One};

    fn fr(n: u64) -> Fr {
        Fr::from(n)
    }

    fn get_domain(n: usize) -> Vec<Fr> {
        let omega = Fr::get_root_of_unity(n as u64).unwrap();
        (0..n).map(|i| omega.pow(&[i as u64])).collect()
    }

    #[test]
    fn test_generate_sigma_mapping_and_get_sigma_maps() {
        // Witness with 3 gates → 3 * 3 = 9 wires
        let witness = Witness {
            a: vec![fr(1), fr(2), fr(3)],
            b: vec![fr(2), fr(1), fr(6)],
            c: vec![fr(7), fr(8), fr(1)],
        };

        // Wiring groups:
        // a0 = b1 = c2 → [0, 4, 8]
        // b0 = a1      → [1, 3]
        let wiring = vec![vec![0, 4, 8], vec![1, 3]];

        let permutation = Permutation::new(witness, wiring, fr(2), fr(3));

        // Test full sigma mapping
        let sigma = permutation.generate_sigma_mapping();

        // Expected:
        // 0 → 4, 4 → 8, 8 → 0
        // 1 → 3, 3 → 1
        let expected_sigma = vec![4, 3, 2, 1, 8, 5, 6, 7, 0];
        assert_eq!(sigma, expected_sigma);

        // Test chunked sigma maps
        let (a_sigma, b_sigma, c_sigma) = permutation.get_sigma_maps();
        assert_eq!(a_sigma, vec![4, 3, 2]);
        assert_eq!(b_sigma, vec![1, 8, 5]);
        assert_eq!(c_sigma, vec![6, 7, 0]);
    }

    #[test]
    fn test_no_wiring() {
        let witness = Witness {
            a: vec![fr(1); 2],
            b: vec![fr(2); 2],
            c: vec![fr(3); 2],
        };
        let wiring = vec![];

        let permutation = Permutation::new(witness, wiring, fr(2), fr(3));

        let sigma = permutation.generate_sigma_mapping();
        assert_eq!(sigma, (0..6).collect::<Vec<_>>());

        let (a_sigma, b_sigma, c_sigma) = permutation.get_sigma_maps();
        assert_eq!(a_sigma, vec![0, 1]);
        assert_eq!(b_sigma, vec![2, 3]);
        assert_eq!(c_sigma, vec![4, 5]);
    }

    #[test]
    fn test_get_sigma_polynomials_correctness() {
        let domain = get_domain(4); // size 4
        let padding_index = domain.len() - 1;
        let mappings = (
            vec![2, 0, 1, 0], // a_sigma map
            vec![1, 2, 0, 0], // b_sigma map
            vec![0, 2, 1, 0], // c_sigma map
        );

        let perm = Permutation {
            witness: Witness {
                a: vec![],
                b: vec![],
                c: vec![],
            },
            wiring: vec![],
            k1: fr(2),
            k2: fr(3),
        };

        let (a_sigma, b_sigma, c_sigma) = perm.get_sigma_polynomials(mappings.clone(), &domain);

        let eval_a = fft(&a_sigma, domain[1]);
        let eval_b = fft(&b_sigma, domain[1]);
        let eval_c = fft(&c_sigma, domain[1]);

        assert_eq!(
            eval_a,
            mappings.0.iter().map(|&i| domain[i]).collect::<Vec<_>>()
        );
        assert_eq!(
            eval_b,
            mappings.1.iter().map(|&i| domain[i]).collect::<Vec<_>>()
        );
        assert_eq!(
            eval_c,
            mappings.2.iter().map(|&i| domain[i]).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_sigma_maps_match_expected_permutation_cycles() {
        let n = 4;
        let domain = get_domain(n);

        // Dummy witness with distinct values
        let witness = Witness {
            a: vec![fr(1), fr(2), fr(10), fr(4)],
            b: vec![fr(5), fr(1), fr(7), fr(10)],
            c: vec![fr(9), fr(10), fr(11), fr(12)],
        };

        let k1 = Fr::from(2u64);
        let k2 = Fr::from(3u64);
        let wiring = vec![vec![0, 5], vec![2, 7, 9]]; // i.e., sigma[i] = i
        let perm = Permutation {
            witness: witness.clone(),
            wiring,
            k1,
            k2,
        };

        // Identity sigma maps
        let sigma_maps = perm.get_sigma_maps();

        let sigma_polys = perm.get_sigma_polynomials(sigma_maps.clone(), &domain);

        let gamma = fr(13);
        let beta = fr(17);
        let z_poly = perm.calculate_rolling_product(sigma_polys, &domain, gamma, beta);

        // Forward FFT of Z(X) to get back evaluations
        let z_eval = fft(&z_poly, domain[1]);

        // Check recurrence:
        assert_eq!(z_eval.len(), domain.len());
        assert_eq!(z_eval[0], Fr::one());

        let result = verify_permutation_argument(
            &domain,
            &z_eval,
            &witness,
            &sigma_maps,
            k1,
            k2,
            beta,
            gamma,
        );
        assert_eq!(Ok(()), result);
    }

    #[test]
    fn test_zx_satisfies_recurrence_on_padded_domain() {
        let n = 8;
        let domain = get_domain(n);

        // witness padded with zeroes
        let witness = Witness {
            a: vec![fr(1), fr(2), fr(3), fr(4), fr(5), fr(6), fr(0), fr(0)],
            b: vec![fr(7), fr(2), fr(1), fr(10), fr(11), fr(12), fr(0), fr(0)],
            c: vec![fr(9), fr(8), fr(13), fr(14), fr(1), fr(9), fr(0), fr(0)],
        };

        // Copy constraints:
        // a0 = b2 = c4 → [0, 10, 20]
        // a1 = b1      → [1, 9]
        // c0 = c5      → [16, 21]
        // wire indexing: a_i = i, b_i = n + i, c_i = 2n + i
        let wiring = vec![vec![0, 10, 20], vec![1, 9], vec![16, 21]];

        let k1 = fr(2);
        let k2 = fr(3);
        let perm = Permutation::new(witness.clone(), wiring, k1, k2);

        let sigma_maps = perm.get_sigma_maps();

        let gamma = fr(9);
        let beta = fr(6);

        let sigma_polys = perm.get_sigma_polynomials(sigma_maps.clone(), &domain);
        let z_poly = perm.calculate_rolling_product(sigma_polys, &domain, gamma, beta);

        // Z evals
        let z_eval = fft(&z_poly, domain[1]);
        assert_eq!(z_eval.len(), domain.len());
        assert_eq!(z_eval[0], fr(1));

        let result = verify_permutation_argument(
            &domain,
            &z_eval,
            &witness,
            &sigma_maps,
            k1,
            k2,
            beta,
            gamma,
        );
        assert_eq!(Ok(()), result);
    }

    #[test]
    fn test_rolling_product_failure() {
        let n = 4;
        let domain = get_domain(n);

        let witness = Witness {
            a: vec![fr(1), fr(2), fr(3), fr(4)],
            b: vec![fr(1), fr(6), fr(7), fr(8)],
            c: vec![fr(9), fr(10), fr(11), fr(12)],
        };

        let wiring = vec![vec![0, 5]];

        let k1 = fr(2);
        let k2 = fr(3);
        let perm = Permutation::new(witness.clone(), wiring, k1, k2);
        let sigma_maps = perm.get_sigma_maps();

        let beta = fr(6);
        let gamma = fr(9);

        let sigma_polys = perm.get_sigma_polynomials(sigma_maps.clone(), &domain);
        let mut z_poly = perm.calculate_rolling_product(sigma_polys, &domain, gamma, beta);

        // Tamper Z polynomial
        z_poly[1] += fr(1);

        let z_eval = fft(&z_poly, domain[1]);

        let result = verify_permutation_argument(
            &domain,
            &z_eval,
            &witness,
            &sigma_maps,
            k1,
            k2,
            beta,
            gamma,
        );
        assert!(result.is_err(), "Expected verification failure, but got Ok");
    }
}
