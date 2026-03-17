# Full SLICOT 1:1 Implementation (No Stubs)

## Requirement
DO NOT use stubs. Implement FULL functionality 1:1 from SLICOT Fortran.

## Dependency chain for AB09AD (balanced truncation)

1. **SB03OV** – complex plane rotation (small; uses DLAPY3)
2. **SB03OY** – 2×2 Lyapunov solver (uses DLANV2, SB03OV, DLAPY2, DLAMCH, DLABAD)
3. **SB03OR** – Sylvester solver (called by SB03OT)
4. **SB03OT** – Lyapunov for Schur form (calls SB03OY, SB03OR, MB04ND, MB04OD)
5. **SB03OU** – Lyapunov wrapper (QR/RQ of B, SB03OT, sign fix)
6. **MB03UD** – SVD of upper triangular matrix (for Hankel singular values)
7. **AB09AX** – balance & truncate core (SB03OU×2, MA02DD, MA02AD, MB03UD, DGEQRF, …)
8. **AB09AD** – top-level (TB01ID, TB01WD, AB09AX)

## Auxiliary needed
- DLAPY2, DLAPY3, DLAMCH, DLABAD (machine constants / norms)
- DLANV2 (2×2 eigenvalue; LAPACK) or equivalent
- DGEQRF, DGERQF, DORMQR (QR/RQ) – nalgebra or manual
- MB04ND, MB04OD (used in SB03OT)

## Implementation order

**See [plans/slicot-resolution-order.md](slicot-resolution-order.md)** for the full dependency graph and resolution order (all callees before callers).

Abbreviated sequence:
1. Tier 0: Primitives (DLAMCH, DLABAD, DLAPY2, DLAPY3, DLANV2, DLASY2, DLARFG, BLAS)
2. Tier 1: SB03OV ✅, SB04PX ✅, MB04NY ✅, MB04OY ✅
3. Tier 2: SB03OY ✅, SB03OR ✅, MB04ND ✅, MB04OD ✅
4. Tier 3: SB03OT ✅, SB03OU
5. Tier 4: MB03UD ✅
6. Tier 5: AB09AX
7. Tier 6: AB09AD (TB01ID, TB01WD already in crate)

## Status
- In progress: SB03OU, AB09AX, AB09AD.
- SB03OY (discrete alpha fix), SB03OT (full port), MB03UD (SVD via nalgebra), SB03OV + DLAPY3: done.
