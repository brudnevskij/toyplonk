use ark_ff::Field;

#[derive(Clone, Debug)]
pub struct Gate<F: Field> {
    q_l: F,
    q_r: F,
    q_m: F,
    q_o: F,
    q_c: F,
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
