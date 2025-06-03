
# 🔒 toy-plonk

A minimal, educational implementation of the [PLONK](https://eprint.iacr.org/2019/953) zero-knowledge proving system built from scratch in Rust.

This project is designed to help understand how modern zk-SNARKs like PLONK work by implementing them step by step, starting from constraint systems all the way to proof generation and verification.

---

## 🧱 Components

| Module         | Description |
|----------------|-------------|
| `constraint_system.rs` | Gate definitions, wire selectors, and arithmetic constraints |
| `circuit.rs`           | Circuit construction and witness assignment |
| `permutation.rs`       | Implements the copy-constraint permutation argument |
| `fft.rs`               | Cooley–Tukey FFT and inverse FFT for evaluation domains |
| `commitment.rs`        | Polynomial commitment scheme (planned: KZG via `arkworks`) |
| `prover.rs`            | Proof generation logic (PLONK-style arithmetization) |
| `verifier.rs`          | Proof verification logic |
| `tests/`               | Unit and integration tests for all components |

---

## 🔍 Scope & Focus

This is a **toy project**, focused on:
- Learning and understanding ZK proving systems
- Learning `arkworks` ecosystem
- Keeping code readable and educational

**Not** focused on:
- Performance
- Production cryptographic security
- General-purpose circuits or universality

---

## 🧪 Current Status

- ✅ FFT & IFFT tested and working
- ⏳ Permutation argument in progress
- ⏳ Gate constraint evaluation
- ⏳ Proof construction and verification
