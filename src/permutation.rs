use crate::witness::Witness;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_std::iterable::Iterable;
use itertools::izip;
use crate::fft::{fft, inverse_fft};

pub struct Permutation<F: Field> {
    pub witness: Witness<F>,
    pub wiring: Vec<Vec<usize>>,
}

impl<F: Field> Permutation<F> {
    fn new(witness: Witness<F>, wiring: Vec<Vec<usize>>) -> Permutation<F> {
        Self { witness, wiring }
    }

    fn generate_sigma_mapping(&self) -> Vec<usize> {
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

    fn get_sigma_maps(&self) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let sigma = self.generate_sigma_mapping();
        let mut a_sigma = vec![];
        let mut b_sigma = vec![];
        let mut c_sigma = vec![];

        for chunk in sigma.chunks_exact(3) {
            a_sigma.push(chunk[0]);
            b_sigma.push(chunk[1]);
            c_sigma.push(chunk[2]);
        }

        (a_sigma, b_sigma, c_sigma)
    }

    fn interpolate_sigma_from_mapping(&self, mapping: Vec<usize>, domain: &Vec<F>) -> Vec<F> {
        if mapping.len() != domain.len() {
            panic!("Evaluation domain and mapping must be the same length, domain.len() = {}, mapping.len() = {}", domain.len(), mapping.len());
        }

        let evaluations = mapping.iter().map(|x| domain[*x]).collect::<Vec<_>>();
        inverse_fft(&evaluations, domain[1])
    }

    fn get_sigma_polynomials(&self, mappings: (Vec<usize>, Vec<usize>, Vec<usize>), domain: &Vec<F>) -> (Vec<F>, Vec<F>, Vec<F>) {
        let a_sigma = self.interpolate_sigma_from_mapping(mappings.0, domain);
        let b_sigma = self.interpolate_sigma_from_mapping(mappings.1, domain);
        let c_sigma = self.interpolate_sigma_from_mapping(mappings.2, domain);
        (a_sigma, b_sigma, c_sigma)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;

    fn fr(n: u64) -> Fr {
        Fr::from(n)
    }

    #[test]
    fn test_generate_sigma_mapping_and_get_sigma_maps() {
        // Witness with 3 gates → 3 * 3 = 9 wires
        let witness = Witness {
            a: vec![fr(1), fr(2), fr(3)],
            b: vec![fr(4), fr(5), fr(6)],
            c: vec![fr(7), fr(8), fr(9)],
        };

        // Wiring groups:
        // a0 = b1 = c2 → [0, 4, 8]
        // b0 = a1      → [1, 3]
        let wiring = vec![vec![0, 4, 8], vec![1, 3]];

        let permutation = Permutation::new(witness, wiring);

        // Test full sigma mapping
        let sigma = permutation.generate_sigma_mapping();

        // Expected:
        // 0 → 4, 4 → 8, 8 → 0
        // 1 → 3, 3 → 1
        let expected_sigma = vec![4, 3, 2, 1, 8, 5, 6, 7, 0];
        assert_eq!(sigma, expected_sigma);

        // Test chunked sigma maps
        let (a_sigma, b_sigma, c_sigma) = permutation.get_sigma_maps();
        assert_eq!(a_sigma, vec![4, 1, 6]);
        assert_eq!(b_sigma, vec![3, 8, 7]);
        assert_eq!(c_sigma, vec![2, 5, 0]);
    }

    #[test]
    fn test_no_wiring() {
        let witness = Witness {
            a: vec![fr(1); 2],
            b: vec![fr(2); 2],
            c: vec![fr(3); 2],
        };
        let wiring = vec![];

        let permutation = Permutation::new(witness, wiring);

        let sigma = permutation.generate_sigma_mapping();
        assert_eq!(sigma, (0..6).collect::<Vec<_>>());

        let (a_sigma, b_sigma, c_sigma) = permutation.get_sigma_maps();
        assert_eq!(a_sigma, vec![0, 3]);
        assert_eq!(b_sigma, vec![1, 4]);
        assert_eq!(c_sigma, vec![2, 5]);
    }
}

