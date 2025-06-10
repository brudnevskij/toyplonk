use ark_ff::Field;

#[derive(Clone, Debug)]
pub struct Witness<F: Field> {
    pub a: Vec<F>,
    pub b: Vec<F>,
    pub c: Vec<F>,
}
