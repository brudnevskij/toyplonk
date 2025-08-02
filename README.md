# 🔒 toyplonk

A minimal and educational implementation of the [PLONK](https://eprint.iacr.org/2019/953) zero-knowledge proving system written from scratch in Rust.

This project is designed for learning and exploring the internals of zk-SNARKs like PLONK, focusing on clarity, correctness, and a hands-on understanding of polynomial-based proof systems.

---

## 🧱 Components

| Module          | Description                      |
|------------------|----------------------------------|
| `circuit.rs`     | Circuit definition, selector polynomials, public inputs, witness wires |
| `gate.rs`        | Gate types (addition, multiplication, public input) and interface logic |
| `permutation.rs` | Copy constraints via permutation argument and σ-polynomials |
| `fft.rs`         | Cooley–Tukey FFT and IFFT for domain interpolation and evaluation |
| `prover.rs`      | Implements the prover side of PLONK: quotient poly, linearization, etc. |
| `verifier.rs`    | Verifies PLONK proofs using KZG commitments and pairing checks |
| `witness.rs`     | Stores witness wire values for each gate |
| `transccript.rs` | Fiat–Shamir heuristic            |
| `main.rs`        | A full example: builds and proves a non-trivial arithmetic circuit |

---

## ✨ Example Circuit

The main function proves knowledge of private inputs satisfying:

```
c = ab + a + b + (a + b)^2
```

for public inputs `a = 2`, `b = 3`, `c = 36`, over 8 gates, including 3 public input gates.  
The full proof is generated and successfully verified using KZG commitments over BLS12-381.

---

## 🎯 Goals

This project is:

✅ **Educational** - built for learning PLONK internals, not for speed  
✅ **Minimal** - fewer abstractions, clearer math  
✅ **Hands-on** - full pipeline: build circuit → generate proof → verify it

Not intended for:

❌ Production use  
❌ Generic circuit support or universality  
❌ High performance or prover speed

---

## ✅ Status

- ✅ Fully working PLONK prover & verifier
- ✅ KZG polynomial commitment scheme via `arkworks`
- ✅ Public input encoding and gate arithmetization
- ✅ FFT and inverse FFT implemented and tested
- ✅ Full example in `main.rs`