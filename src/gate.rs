use ark_ff::Field;

#[derive(Clone, Debug)]
pub struct Gate<F: Field> {
    pub
    q_l: F,
    pub q_r: F,
    pub q_m: F,
    pub q_o: F,
    pub q_c: F,
}

impl<F: Field> Gate<F> {
    pub fn new(q_l: F, q_r: F, q_m: F, q_o: F, q_c: F) -> Self {
        Self {
            q_l,
            q_r,
            q_m,
            q_o,
            q_c,
        }
    }
    pub fn addition_gate(
        left_coefficient: F,
        right_coefficient: F,
        output_coefficient: F,
    ) -> Gate<F> {
        Self {
            q_l: left_coefficient,
            q_r: right_coefficient,
            q_m: F::zero(),
            q_o: -output_coefficient,
            q_c: F::zero(),
        }
    }
    pub fn simple_addition_gate() -> Gate<F> {
        Self::addition_gate(F::one(), F::one(), F::one())
    }

    pub fn mul_gate(mul_coefficient: F, output_coefficient: F) -> Self {
        Self {
            q_l: F::zero(),
            q_r: F::zero(),
            q_m: mul_coefficient,
            q_o: -output_coefficient,
            q_c: F::zero(),
        }
    }

    pub fn simple_mul_gate() -> Gate<F> {
        Self::mul_gate(F::one(), F::one())
    }

    pub fn is_satisfied(&self, a: F, b: F, c: F) -> bool {
        let res = self.q_l * a + self.q_r * b + self.q_m * a * b + self.q_o * c + self.q_c;
        res.is_zero()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Fr;

    #[test]
    fn test_simple_addition_gate_is_satisfied() {
        use ark_bls12_381::Fr;
        let gate = Gate::<Fr>::simple_addition_gate();

        let a = Fr::from(3);
        let b = Fr::from(4);
        let c = Fr::from(7);

        assert!(gate.is_satisfied(a, b, c));
    }

    #[test]
    fn test_simple_mul_gate_is_satisfied() {
        use ark_bls12_381::Fr;
        let gate = Gate::<Fr>::simple_mul_gate();

        let a = Fr::from(2);
        let b = Fr::from(5);
        let c = Fr::from(10);

        assert!(gate.is_satisfied(a, b, c));
    }

    #[test]
    fn test_add_gate_fails_on_wrong_values() {
        use ark_bls12_381::Fr;
        let gate = Gate::<Fr>::simple_addition_gate();

        let a = Fr::from(1);
        let b = Fr::from(2);
        let c = Fr::from(5); // wrong!

        assert!(!gate.is_satisfied(a, b, c));
    }

    #[test]
    fn test_weighted_mul_gate() {
        use ark_bls12_381::Fr;
        let gate = Gate::<Fr>::mul_gate(Fr::from(3), Fr::from(2)); // 3ab = 2c

        let a = Fr::from(2);
        let b = Fr::from(5);
        let c = Fr::from(15); // 3·2·5 = 2·15 = 30, OK

        assert!(gate.is_satisfied(a, b, c));
    }
}
