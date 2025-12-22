# 🔒 toyplonk – A Minimal Rust Implementation of PLONK zk-SNARK 

This project implements the full proving and verification pipeline of the [PLONK](https://eprint.iacr.org/2019/953) zero-knowledge proof system using the Rust `arkworks` ecosystem.

Built from scratch for clarity and education, it aims to demystify PLONK by walking through its polynomial-based components: constraint encoding, permutation argument, quotient polynomial, KZG commitments, and pairing-based verification.

---

## 🚀 Features

- ✅ Gate-based constraint system with selector polynomials
- ✅ Support for public input encoding via identity gates
- ✅ Full PLONK prover with permutation & quotient polynomials
- ✅ Linearization and batched opening proofs
- ✅ Verifier with KZG pairing equation
- ✅ Generic over any pairing-friendly curve via arkworks traits
- ✅ Example circuit proving `c = ab + a + b + (a + b)^2`

---

## 🧠 Motivation & Scope

This was built to **understand the internals** of PLONK from scratch, without macros or abstractions hiding the protocol steps.
I didn’t focus on performance or general-purpose circuits. The goal is clarity, debuggability, and correctness.

---

## 🧪 Example: Proving c = ab + a + b + (a + b)^2

The demo in `main.rs` builds an 8-gate circuit with 3 public inputs and verifies a PLONK proof for `a=2`, `b=3`, `c=36`.

```bash
cargo run
```

Expected output:
```
📄 PLONK proof: ACCEPTED ✅
```

---

## 📚 Learning Resources

Here are the references that helped me build and understand this implementation:

- [The PLONK paper](https://eprint.iacr.org/2019/953) - Original specification.
- [How PLONK Works: Part 2](https://scryptplatform.medium.com/how-plonk-works-part-2-1072dcd7634a) - Helped me in understanding permutation check.
- [PLONK by Hand](https://research.metastate.dev/plonk-by-hand-part-1/) - Small PLONK example, fun and easy to read.
- [All you wanted to know about Plonk](https://blog.lambdaclass.com/all-you-wanted-to-know-about-plonk/) - Article by LambdaClass, helped me understand how public input is integrated to plonk, also has some cool matrix based examples.
- [Notes on Plonk Prover's and Verifier's Algorithm](https://hackmd.io/@aztec-network/ByiUK_Plt) - Helped me identify core algebraic properties during proof construction, for testing purposes.
---

## 🧱 Structure

| Module             | Role                                      |
|--------------------|-------------------------------------------|
| `circuit.rs`       | Circuit structure, selector & witness polynomials |
| `gate.rs`          | Basic arithmetic gates (add, mul, identity)       |
| `permutation.rs`   | Permutation argument & σ-polynomials       |
| `prover.rs`        | Prover implementation (quotient, linearization, KZG) |
| `verifier.rs`      | Verifier using KZG opening and batching    |
| `fft.rs`           | FFT & IFFT routines                        |
| `main.rs`          | Example proving circuit                    |

---

Feel free to explore, debug, break, or extend. The goal is understanding.

This project implements the full proving and verification pipeline of the [PLONK](https://eprint.iacr.org/2019/953) zero-knowledge proof system using the Rust `arkworks` ecosystem.

Built from scratch for clarity and education, it aims to demystify PLONK by walking through its polynomial-based components: constraint encoding, permutation argument, quotient polynomial, KZG commitments, and pairing-based verification.
