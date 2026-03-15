# slicot-rs

> *"There is nothing more frightful than ignorance in action."*  
> — Goethe  
>  
> *"The fundamental cause of trouble in the world today is that the stupid are cocksure while the intelligent are full of doubt."*  
> — Bertrand Russell  

This crate exists because someone thought it was a good idea to rewrite **SLICOT**—the Subroutine Library in Control Theory, a decades-old, battle-tested Fortran pile used in aerospace, process control, and systems engineering—in **Rust**, and to do it with **AI**. The hubris required to aim for “bulletproof” parity with SLICOT in a new language and with generated code is considerable. We’re doing it anyway. Use at your own risk; validate against the [official SLICOT library](https://www.slicot.org/) and documentation.

---

## What is this?

**slicot-rs** is a pure Rust, 1:1 mapping of SLICOT routines. It uses [nalgebra](https://crates.io/crates/nalgebra) for linear algebra and does **not** call LAPACK/BLAS via FFI. The goal is to provide the same numerical control and systems primitives (real dense LTI systems, transformations, model reduction, identification, etc.) in a safe, dependency-light Rust API.

Routine status and remaining work are tracked in **[plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md)**. The full SLICOT → Rust mapping is in [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md).

---

## Examples and real-world use cases

### Control design and analysis

- **Controllability / observability**  
  Use **AB01ND** / **AB01OD** (and related) to reduce \((A, B)\) or \((A, C)\) to staircase or block Hessenberg form and detect controllable/unobservable subspaces. Essential for minimal realizations and before placing poles or designing observers.

- **Transfer-function rank and structure**  
  **AB08MD** computes the normal rank of the transfer matrix \(G(s) = C(sI - A)^{-1}B + D\). Used in decoupling control, feedforward design, and checking invertibility of systems.

- **Model reduction (balanced truncation, singular perturbation)**  
  The **AB09\*** family (AB09AD, AB09BD, …) implements balanced truncation and related methods to reduce state dimension while approximating input–output behavior. Used in controller simplification and large-scale simulation.

### Signal and polynomial operations

- **Convolution / deconvolution**  
  **DE01OD** performs convolution or deconvolution of two real sequences (e.g. impulse responses or polynomial coefficients). **DE01PD** does weighted convolution/deconvolution.

- **Discrete-time transforms**  
  **DG01MD**, **DG01ND**, **DG01NY**, **DG01OD** provide direct/inverse transforms and weighting for discrete-time signals (e.g. pre-processing for identification).

- **Polynomial / rational operations**  
  **MC01\***-style operations (when available) and **DG01** help with polynomial arithmetic and frequency-domain manipulations used in transfer-function and filter design.

### Matrix and linear algebra utilities

- **Symmetric / skew-symmetric helpers**  
  **MA02ED** completes the opposite triangle of a symmetric matrix from the given triangle. **MB01MD** computes \(y = \alpha A x + \beta y\) for skew-symmetric \(A\). These underpin many Lyapunov and Riccati solvers and structure-preserving factorizations.

- **Matrix copy and layout**  
  **DLACPY_SLC** copies full or triangular parts of a matrix (LAPACK-style); **DLATZM** applies a Householder-like transformation. Used inside higher-level SLICOT routines.

### Identification and filtering

- **Subspace identification**  
  The **IB01\*** suite (IB01AD–IB01RD) and **IB03\*** (IB03AD, IB03BD) support subspace-based identification of state-space models from input–output data—e.g. identifying a MIMO LTI model from measured sequences.

- **Kalman filtering**  
  **FB01\*** (FB01QD, FB01RD, FB01SD, FB01TD, FB01VD) provide one-step Kalman filter updates (covariance and state). Used in state estimation and recursive identification.

### Benchmarks and system norms

- **Norms and distances**  
  **AB13DD**, **AB13ED**, **AB13FD**, **AB13HD**, **AB13MD** compute \(\mathcal{H}_2\)/\(\mathcal{H}_\infty\) norms and related measures for LTI systems. Used in robustness analysis and model-reduction error bounds.

- **Benchmark drivers**  
  **BB01AD–BB04AD** and **BD01AD**, **BD02AD** are benchmark-style drivers for performance and regression testing against SLICOT.

---

## Quick example

```rust
use nalgebra::DMatrix;
use slicot_rs::ma02::{ma02ed, Ma02EdUplo};

// Build a 2×2 symmetric matrix from upper triangle only
let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]);
let info = ma02ed(Ma02EdUplo::Upper, &mut a);
assert_eq!(info, 0);
// Lower triangle filled by symmetry: a(1,0) == 2.0
```

More examples and a function-by-function reference (with MathJax and SLICOT doc pointers) are in **[docs/FUNCTIONS.md](docs/FUNCTIONS.md)**.

---

## Documentation and references

| Resource | Description |
|----------|-------------|
| [docs/FUNCTIONS.md](docs/FUNCTIONS.md) | Function index with examples, math, and SLICOT references |
| [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md) | SLICOT routine → Rust module/function mapping |
| [plans/remaining-slicot-functions-todo.md](plans/remaining-slicot-functions-todo.md) | TODO list for remaining SLICOT functions |
| [SLICOT Library](https://www.slicot.org/) | Official SLICOT site and documentation |
| [SLICOT Routine Index](https://www.slicot.org/objects/software/shared/libindex.html) | Alphabetical list of routines |

---

## License

BSD-3-Clause. See [LICENSE](LICENSE) or the crate root.
