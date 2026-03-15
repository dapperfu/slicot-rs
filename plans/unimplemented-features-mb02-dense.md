# Unimplemented features: MB02 dense fallbacks

## Goal
Implement unimplemented MB02 stub routines with dense fallbacks so each routine has a working implementation (correctness over performance).

## Scope
- **MB02FD**: Cholesky factorization of symmetric positive definite matrix (stub API: n, a, b → b = upper Cholesky R, R'*R = A).
- **MB02GD**: Same dense Cholesky fallback for banded/s.p.d. (stub API identical: n, a, b).
- Other MB02 stubs (CU, CV, CX, CY, DD, HD, …): leave as stubs or add same-pattern dense fallbacks in follow-up.

## Decisions
- Stub API kept: `(n, a: &DMatrix, b: &mut DMatrix)`; a = input matrix, b = output.
- MB02FD/MB02GD: interpret a as n×n s.p.d. matrix; on success b gets upper Cholesky R (b = R where A = R'*R). Return 0 on success, 1 if not positive definite.
- No new dependencies; use nalgebra's `cholesky()`.

## TODOs
- [x] Save plan to plans/
- [x] Implement MB02FD
- [x] Implement MB02GD
- [x] Update SLICOT_MAPPING and FEATURES
- [ ] Commit and push per workflow
