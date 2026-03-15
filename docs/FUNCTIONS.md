# slicot-rs function reference

This document lists every implemented routine in **slicot-rs**, with short descriptions, mathematical notation (MathJax), and usage examples. For full specifications and algorithms, see the [SLICOT Library](https://www.slicot.org/) and the [SLICOT Routine Index](https://www.slicot.org/objects/software/shared/libindex.html).

**Rendering MathJax:** Use a viewer that supports MathJax (e.g. GitHub with a MathJax extension, or a local Markdown renderer with MathJax enabled) to see matrices and equations.

---

## AB01 — State-space analysis (canonical / staircase forms)

### AB01ND — Orthogonal controllability staircase form (multi-input)

- **Rust:** `ab01::ab01nd`
- **SLICOT:** [AB01ND](https://www.slicot.org/objects/software/shared/libindex.html) — Orthogonal controllability staircase form for multi-input system.

Reduces the pair \((A, B)\) to an orthogonal canonical form (block Hessenberg). In state-space form:

\[
\dot{x} = A x + B u, \quad A \in \mathbb{R}^{n \times n}, \quad B \in \mathbb{R}^{n \times m}.
\]

The routine computes an orthogonal \(Z\) such that \((Z^\top A Z, Z^\top B)\) is in staircase form, revealing the controllable subspace dimension `ncont`.

**Example:**

```rust
use nalgebra::DMatrix;
use slicot_rs::ab01::{ab01nd, JobZ};

let n = 4;
let m = 2;
let mut a = DMatrix::zeros(n, n);
let mut b = DMatrix::zeros(n, m);
let mut ncont = 0;
let mut indcon = 0;
let mut nblk = vec![0i32; n];
let info = ab01nd(JobZ::No, n, m, &mut a, &mut b, &mut ncont, &mut indcon, &mut nblk, None, 1e-10);
// info: 0 = success, 1 = not implemented, < 0 = invalid argument
```

---

### AB01OD — Staircase form for (A,B) with forward/backward stages

- **Rust:** `ab01::ab01od`
- **SLICOT:** [AB01OD](https://www.slicot.org/objects/software/shared/libindex.html) — Staircase form for multi-input system using orthogonal transformations.

Reduces \((A, B)\) to upper staircase form via orthogonal \(U, V\). Same state-space pair as above; options control forward/backward stages and whether to accumulate \(U\) and \(V\).

**Example:**

```rust
use slicot_rs::ab01::{ab01od, Stages, JobUV};
use nalgebra::DMatrix;

let (n, m) = (3, 1);
let mut a = DMatrix::zeros(n, n);
let mut b = DMatrix::zeros(n, m);
let mut ncont = 0;
let mut indcon = 0;
let mut kstair = vec![0i32; n + 1];
let info = ab01od(
    Stages::Forward,
    JobUV::No,
    JobUV::No,
    n, m,
    &mut a, &mut b,
    None, None,
    &mut ncont, &mut indcon,
    &mut kstair,
    0.0,
);
```

---

## AB04 — Continuous/discrete time conversion

### AB04MD — Discrete ↔ continuous time by bilinear transformation

- **Rust:** `ab04::ab04md`
- **SLICOT:** [AB04MD](https://www.slicot.org/objects/software/shared/libindex.html) — Discrete-time ↔ continuous-time conversion by bilinear transformation.

Converts between continuous-time \(\dot{x}=Ax+Bu\) and discrete-time \(x_{k+1}=A_d x_k + B_d u_k\) using the bilinear (Tustin) transformation. The mapping involves matrices \(A,B,C,D\) and sampling time.

**Example:** (see crate tests or `docs/SLICOT_MAPPING.md` for full usage.)

---

## AB05 — Interconnections of subsystems

### AB05ND — Feedback interconnection

- **Rust:** `ab05::ab05nd`
- **SLICOT:** [AB05ND](https://www.slicot.org/objects/software/shared/libindex.html) — Feedback inter-connection of two systems in state-space form.

Forms the closed-loop system when two LTI systems are connected in feedback. Given \(G_1\) and \(G_2\), computes the state-space representation of the feedback loop.

### AB05OD — Rowwise concatenation

- **Rust:** `ab05::ab05od`
- **SLICOT:** [AB05OD](https://www.slicot.org/objects/software/shared/libindex.html) — Rowwise concatenation of two systems.

### AB05PD — Parallel interconnection

- **Rust:** `ab05::ab05pd`
- **SLICOT:** [AB05PD](https://www.slicot.org/objects/software/shared/libindex.html) — Parallel inter-connection of two systems.

### AB05QD — Appending two systems

- **Rust:** `ab05::ab05qd`
- **SLICOT:** [AB05QD](https://www.slicot.org/objects/software/shared/libindex.html) — Appending two systems in state-space form.

### AB05RD — Closed-loop (output + state feedback)

- **Rust:** `ab05::ab05rd`
- **SLICOT:** [AB05RD](https://www.slicot.org/objects/software/shared/libindex.html) — Closed-loop system for mixed output and state feedback.

### AB05SD — Closed-loop (output feedback)

- **Rust:** `ab05::ab05sd`
- **SLICOT:** [AB05SD](https://www.slicot.org/objects/software/shared/libindex.html) — Closed-loop system for output feedback control law.

---

## AB07 — Inverse and dual systems

### AB07MD — Dual system

- **Rust:** `ab07::ab07md`
- **SLICOT:** [AB07MD](https://www.slicot.org/objects/software/shared/libindex.html) — Dual of a given state-space representation.

Given \((A, B, C, D)\), forms the dual realization \((A^\top, C^\top, B^\top, D^\top)\).

\[
G(s) = C(sI - A)^{-1} B + D \quad \Rightarrow \quad G^T(s) = B^\top (sI - A^\top)^{-1} C^\top + D^\top.
\]

### AB07ND — Inverse system

- **Rust:** `ab07::ab07nd`
- **SLICOT:** [AB07ND](https://www.slicot.org/objects/software/shared/libindex.html) — Inverse of a given state-space representation.

Computes a state-space realization of \(G^{-1}(s)\) when \(D\) is invertible.

---

## AB08 — Poles, zeros, normal rank

### AB08MD — Normal rank of transfer-function matrix

- **Rust:** `ab08::ab08md`
- **SLICOT:** [AB08MD](https://www.slicot.org/objects/software/shared/libindex.html) — Normal rank of the transfer-function matrix of a state space model.

Computes the normal rank of
\[
G(s) = C(sI - A)^{-1} B + D \in \mathbb{R}^{p \times m}(s).
\]
When \(n=0\), rank is the rank of \(D\); otherwise the routine uses the system pencil. Used for decoupling and feedforward design.

**Example:**

```rust
use slicot_rs::ab08::{ab08md, Ab08MdEquil};

let (n, m, p) = (0, 2, 2);
let d = [1.0, 0.0, 0.0, 1.0]; // 2×2 identity
let mut rank = -1i32;
let mut iwork = [0i32; 4];
let mut dwork = vec![0.0; 64];
let info = ab08md(
    Ab08MdEquil::No, n, m, p,
    &[], 1, &[], 1, &[], 1, &d, 2,
    &mut rank, 0.0, &mut iwork, &mut dwork, 64,
);
assert_eq!(info, 0);
assert_eq!(rank, 2);
```

### AB08MZ, AB08ND, AB08NW, AB08NX, AB08NY, AB08NZ

- **Rust:** `ab08::ab08mz`, `ab08nd`, `ab08nw`, `ab08nx`, `ab08ny`, `ab08nz`
- **SLICOT:** [AB08*](https://www.slicot.org/objects/software/shared/libindex.html) — System zeros and Kronecker structure (real/complex, various pencils).

Routines for zeros and Kronecker structure of the system pencil; used in pole-zero analysis and minimality.

---

## AB09 — Model reduction

### AB09AD through AB09ND

- **Rust:** `ab09::ab09ad`, `ab09ax`, `ab09bd`, `ab09bx`, … `ab09nd`
- **SLICOT:** [AB09*](https://www.slicot.org/objects/software/shared/libindex.html) — Balanced truncation and related model reduction.

Family for balanced truncation, singular perturbation approximation, and Hankel-norm approximation. Given a full-order model \((A,B,C,D)\), produce a reduced model \((\hat{A},\hat{B},\hat{C},\hat{D})\) of order \(r\) with error bounds. For example, balanced truncation uses the controllability and observability Gramians \(P,Q\) satisfying
\[
A P + P A^\top + B B^\top = 0, \qquad A^\top Q + Q A + C^\top C = 0.
\]

---

## AB13 — Norms and distances

### AB13DD, AB13ED, AB13FD, AB13HD, AB13MD

- **Rust:** `ab13::ab13dd`, `ab13ed`, `ab13fd`, `ab13hd`, `ab13md`
- **SLICOT:** [AB13*](https://www.slicot.org/objects/software/shared/libindex.html) — \(\mathcal{H}_2\), \(\mathcal{H}_\infty\) norms and related measures.

Compute system norms (e.g. \(\|G\|_2\), \(\|G\|_\infty\)) and distances for LTI systems; used in robustness and model-reduction validation.

---

## AB8N — Extended AB08 (complex)

### AB8NXZ

- **Rust:** `ab8n::ab8nxz`
- **SLICOT:** [AB8NXZ](https://www.slicot.org/objects/software/shared/libindex.html) — Wrapper/extension for AB08-style computations (complex).

---

## AG07, AG08, AG8B — Generalized state-space and descriptor systems

### AG07BD, AG08BD, AG08BY, AG08BZ, AG8BYZ

- **Rust:** `ag07::ag07bd`; `ag08::ag08bd`, `ag08by`, `ag08bz`; `ag8b::ag8byz`
- **SLICOT:** [AG07/AG08/AG8B](https://www.slicot.org/objects/software/shared/libindex.html) — Descriptor systems and generalized pencils.

Routines for systems \(E\dot{x}=Ax+Bu\), \(y=Cx+Du\), and related pencil manipulations.

---

## BB01–BB04, BD01–BD02 — Benchmarks

### BB01AD, BB02AD, BB03AD, BB04AD, BD01AD, BD02AD

- **Rust:** `bb01::bb01ad`, `bb02::bb02ad`, `bb03::bb03ad`, `bb04::bb04ad`; `bd01::bd01ad`, `bd02::bd02ad`
- **SLICOT:** [BB01–BD02](https://www.slicot.org/objects/software/shared/libindex.html) — Benchmark drivers.

Used for performance and correctness regression against reference SLICOT outputs.

---

## DE01 — Convolution and deconvolution

### DE01OD — Convolution or deconvolution of two sequences

- **Rust:** `de01::de01od`
- **SLICOT:** [DE01OD](https://www.slicot.org/objects/software/shared/libindex.html) — Convolution/deconvolution of two real sequences.

For sequences \((a_k)\), \((b_k)\) of length \(n\), computes either convolution \(c = a * b\) or deconvolution. In polynomial form, multiplication or division of polynomials with these coefficients.

**Example:**

```rust
use slicot_rs::de01::{de01od, De01OdConv};

let n = 4;
let mut a = vec![1.0, 2.0, 3.0, 4.0];
let mut b = vec![0.0; 4];
let info = de01od(De01OdConv::Convolution, n, &mut a, &mut b);
```

### DE01PD — Weighted convolution/deconvolution

- **Rust:** `de01::de01pd`
- **SLICOT:** [DE01PD](https://www.slicot.org/objects/software/shared/libindex.html) — Weighted convolution/deconvolution.

---

## DF01, DG01 — Discrete-time and transform utilities

### DF01MD — Discrete-time frequency response (sine/cosine)

- **Rust:** `df01::df01md`
- **SLICOT:** [DF01MD](https://www.slicot.org/objects/software/shared/libindex.html) — Discrete-time conversion (sine/cosine form).

### DG01MD — Direct transform (real sequences)

- **Rust:** `dg01::dg01md`
- **SLICOT:** [DG01MD](https://www.slicot.org/objects/software/shared/libindex.html) — Direct transform for real sequences.

### DG01ND — Inverse transform (real sequences)

- **Rust:** `dg01::dg01nd`
- **SLICOT:** [DG01ND](https://www.slicot.org/objects/software/shared/libindex.html) — Inverse transform.

### DG01NY — Another direct/inverse transform variant

- **Rust:** `dg01::dg01ny`
- **SLICOT:** [DG01NY](https://www.slicot.org/objects/software/shared/libindex.html) — Transform variant for real sequences.

### DG01OD — Orthogonal transformation (weighting)

- **Rust:** `dg01::dg01od`
- **SLICOT:** [DG01OD](https://www.slicot.org/objects/software/shared/libindex.html) — Orthogonal transformation with optional weighting.

---

## DGEG — Generalized eigenvalue problems (LAPACK-style)

### DGEGS, DGEGV

- **Rust:** `dgeg::dgegs`, `dgeg::dgegv`
- **SLICOT:** [DGEGS/DGEGV](https://www.slicot.org/objects/software/shared/libindex.html) — Generalized eigenvalue problem for matrix pairs \((A,B)\).

Solve \(A x = \lambda B x\). Used inside many SLICOT routines (e.g. pole/zero computation).

---

## DK01 — Discrete-time filtering

### DK01MD — Digital filter (certain filter types)

- **Rust:** `dk01::dk01md`
- **SLICOT:** [DK01MD](https://www.slicot.org/objects/software/shared/libindex.html) — Digital filter implementation.

Applies a specified filter type to a real sequence stored in `a`.

**Example:**

```rust
use slicot_rs::dk01::{dk01md, Dk01MdType};

let n = 8;
let mut a = vec![1.0_f64; n];
let info = dk01md(Dk01MdType::Average, n, &mut a);
```

---

## DLAC, DLAT — LAPACK-style auxiliaries

### DLACPY_SLC — Copy matrix (full or triangle)

- **Rust:** `dlac::dlacpy_slc`
- **SLICOT/LAPACK:** Copy all or part of a matrix (upper, lower, or full).

\[
B \gets \text{copy of } A \text{ (optionally only upper/lower triangle)}.
\]

**Example:**

```rust
use nalgebra::DMatrix;
use slicot_rs::dlac::{dlacpy_slc, DlacpyUplo};

let a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
let mut b = DMatrix::zeros(2, 2);
let info = dlacpy_slc(DlacpyUplo::All, &a, &mut b);
assert_eq!(info, 0);
assert_eq!(b, a);
```

### DLATZM — Apply Householder-like transformation

- **Rust:** `dlat::dlatzm`
- **SLICOT/LAPACK:** Apply a block reflector (Householder-style) to a matrix.

Used in orthogonal factorizations and staircase algorithms.

---

## FB01 — Kalman filter updates

### FB01QD, FB01RD, FB01SD, FB01TD, FB01VD

- **Rust:** `fb01::fb01qd`, `fb01rd`, `fb01sd`, `fb01td`, `fb01vd`
- **SLICOT:** [FB01*](https://www.slicot.org/objects/software/shared/libindex.html) — One-step Kalman filter (covariance and state update).

Used in recursive estimation and subspace identification. Typical update:
\[
\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k (z_k - H_k \hat{x}_{k|k-1}), \qquad P_{k|k} = (I - K_k H_k) P_{k|k-1}.
\]

---

## FD01 — Filter design

### FD01AD

- **Rust:** `fd01::fd01ad`
- **SLICOT:** [FD01AD](https://www.slicot.org/objects/software/shared/libindex.html) — Filter design utility.

---

## IB01 — Subspace identification

### IB01AD through IB01RD

- **Rust:** `ib01::ib01ad`, `ib01bd`, `ib01cd`, `ib01md`, `ib01my`, `ib01nd`, `ib01od`, `ib01oy`, `ib01pd`, `ib01px`, `ib01py`, `ib01qd`, `ib01rd`
- **SLICOT:** [IB01*](https://www.slicot.org/objects/software/shared/libindex.html) — Subspace identification.

Estimate state-space models \((A,B,C,D)\) from input-output data using subspace methods (e.g. N4SID-style). Input: batches of \(u_k\), \(y_k\); output: order and system matrices.

---

## IB03 — Subspace identification (continued)

### IB03AD, IB03BD

- **Rust:** `ib03::ib03ad`, `ib03::ib03bd`
- **SLICOT:** [IB03*](https://www.slicot.org/objects/software/shared/libindex.html) — Additional subspace identification steps.

---

## MA01 — Matrix utilities (scaling, norms)

### MA01AD, MA01BD, MA01BZ, MA01DD, MA01DZ

- **Rust:** `ma01::ma01ad`, `ma01bd`, `ma01bz`, `ma01dd`, `ma01dz`
- **SLICOT:** [MA01*](https://www.slicot.org/objects/software/shared/libindex.html) — Matrix scaling and norm-like operations.

---

## MA02 — Symmetry and storage

### MA02ED — Store by symmetry (complete symmetric matrix)

- **Rust:** `ma02::ma02ed`
- **SLICOT:** [MA02ED](https://www.slicot.org/objects/software/shared/libindex.html) — Store by symmetry; fill the opposite triangle of a symmetric matrix.

Given \(A\) with either upper or lower triangle stored, fills the other triangle so that \(A = A^\top\):
\[
A_{ji} \gets A_{ij} \quad \text{(if upper given)} \qquad \text{or} \qquad A_{ij} \gets A_{ji} \quad \text{(if lower given)}.
\]

**Example:**

```rust
use nalgebra::DMatrix;
use slicot_rs::ma02::{ma02ed, Ma02EdUplo};

let mut a = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 0.0, 3.0]); // upper triangle
let info = ma02ed(Ma02EdUplo::Upper, &mut a);
assert_eq!(info, 0);
assert_eq!(a[(1, 0)], 2.0);
```

### MA02AZ, MA02BZ, MA02CD, MA02CZ, MA02DD, MA02ES, MA02EZ, MA02FD, MA02GD, MA02GZ, MA02NZ, MA02PD, MA02PZ, MA02RD

- **Rust:** `ma02::ma02az`, `ma02bz`, `ma02cd`, `ma02cz`, `ma02dd`, `ma02es`, `ma02ez`, `ma02fd`, `ma02gd`, `ma02gz`, `ma02nz`, `ma02pd`, `ma02pz`, `ma02rd`
- **SLICOT:** [MA02*](https://www.slicot.org/objects/software/shared/libindex.html) — Transpose, conjugate, copy, norm, symmetry and packing utilities.

---

## MB01 — Matrix operations (products, symmetry, packing)

### MB01MD — Skew-symmetric matrix–vector product

- **Rust:** `mb01::mb01md`
- **SLICOT:** [MB01MD](https://www.slicot.org/objects/software/shared/libindex.html) — Skew-symmetric matrix-vector product.

Computes \(y := \alpha A x + \beta y\) with \(A = -A^\top\) (only one triangle of \(A\) stored):
\[
A^\top = -A \quad \Rightarrow \quad (Ax)_i = \sum_j A_{ij} x_j.
\]

**Example:**

```rust
use nalgebra::{DMatrix, DVector};
use slicot_rs::mb01::{mb01md, Mb01MdUplo};

let n = 3;
// Skew-symmetric: store strictly upper triangle (e.g. 0, a12, a13; 0, 0, a23)
let a = DMatrix::from_row_slice(3, 3, &[
    0.0,  1.0, -2.0,
    0.0,  0.0,  1.0,
    0.0,  0.0,  0.0,
]);
let x = DVector::from_vec(vec![1.0, 0.0, 0.0]);
let mut y = DVector::zeros(n);
let info = mb01md(Mb01MdUplo::Upper, 1.0, &a, &x, 0.0, &mut y);
assert_eq!(info, 0);
// y = A*x with A skew-symmetric
```

### MB01LD, MB01ND, MB01OC, MB01OD, MB01OE, MB01OH, MB01OO, MB01OS, MB01OT, MB01RB, MB01RD, MB01RH, MB01RT, MB01RU, MB01RW, MB01RX, MB01RY, MB01SD, MB01SS, MB01TD, MB01UD, MB01UW, MB01UX, MB01UY, MB01UZ, MB01VD, MB01WD, MB01XD, MB01XY, MB01KD

- **Rust:** `mb01::mb01ld`, `mb01nd`, `mb01oc`, `mb01od`, … `mb01xy`, `mb01kd`
- **SLICOT:** [MB01*](https://www.slicot.org/objects/software/shared/libindex.html) — Symmetric/skew-symmetric products, rank-k updates, packing, scaling.

Used throughout SLICOT for Lyapunov/Riccati steps and structure-preserving linear algebra (e.g. \(A X B^\top + B X A^\top\), rank-k updates).

---

## TB01 — Input/output transformations

### TB01MD

- **Rust:** `tb01::tb01md`
- **SLICOT:** [TB01MD](https://www.slicot.org/objects/software/shared/libindex.html) — Input/output scaling or transformation.

---

## Summary table (implemented routines)

| SLICOT   | Rust module | Rust function   | Area                    |
|----------|-------------|-----------------|-------------------------|
| AB01ND   | ab01        | ab01nd          | Staircase (A,B)         |
| AB01OD   | ab01        | ab01od          | Staircase (A,B)         |
| AB04MD   | ab04        | ab04md          | Discrete/continuous     |
| AB05*    | ab05        | ab05nd,…        | Interconnections        |
| AB07MD/ND| ab07        | ab07md, ab07nd  | Dual, inverse           |
| AB08*    | ab08        | ab08md,…        | Rank, zeros             |
| AB09*    | ab09        | ab09ad,…        | Model reduction         |
| AB13*    | ab13        | ab13dd,…        | Norms                   |
| AB8NXZ   | ab8n        | ab8nxz          | AB08 extension          |
| AG07BD,… | ag07, ag08, ag8b | ag07bd,…   | Descriptor systems      |
| BB01AD–BD02AD | bb01–bd02 | bb01ad,…   | Benchmarks              |
| DE01OD/PD| de01        | de01od, de01pd  | Convolution             |
| DF01MD   | df01        | df01md          | Discrete-time           |
| DG01*    | dg01        | dg01md,…        | Transforms              |
| DGEGS/V  | dgeg        | dgegs, dgegv    | Generalized eigen       |
| DK01MD   | dk01        | dk01md          | Digital filter          |
| DLACPY_SLC | dlac      | dlacpy_slc      | Matrix copy             |
| DLATZM   | dlat        | dlatzm          | Block reflector         |
| FB01*    | fb01        | fb01qd,…        | Kalman filter           |
| FD01AD   | fd01        | fd01ad          | Filter design           |
| IB01*    | ib01        | ib01ad,…        | Subspace ID             |
| IB03AD/BD| ib03        | ib03ad, ib03bd  | Subspace ID             |
| MA01*    | ma01        | ma01ad,…        | Matrix utilities        |
| MA02*    | ma02        | ma02ed,…        | Symmetry, transpose     |
| MB01*    | mb01        | mb01md,…        | Products, updates       |
| TB01MD   | tb01        | tb01md          | I/O transformation      |

For the exact list of “done” routines and remaining TODOs, see [plans/remaining-slicot-functions-todo.md](../plans/remaining-slicot-functions-todo.md) and [docs/SLICOT_MAPPING.md](SLICOT_MAPPING.md).
