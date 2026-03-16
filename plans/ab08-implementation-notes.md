# AB08 implementation notes

AB08NW, AB08NY, AB08NZ are required for FORTRAN validation (AB08ND.dat, AB08NW.dat, AB08NZ.dat exist).

## Dependency order

- **AB08NY** is called by AB08NW (and possibly AB08NZ). Implement AB08NY first.
- **AB08NY** (Fortran ~590 lines) calls: MB04ID, DORMQR, DORMRQ, MB03OY, DLAPMT, MB03PY, DLASET.
- **AB08NW** (Fortran ~586 lines) calls: AB08NY, TB01ID, TB01XD, MA02BD, DTZRZF, DORMRZ, DLACPY, DLASET.
- **AB08NZ** uses COMPLEX*16 (complex matrices); different from NW/NY.

Rust crate already has: TB01ID, TB01XD, MA02BD, MB03OY, MB03PY, MB04ID (simplified APIs).

## Work items

1. Port AB08NY: compound matrix (B A; D C) reduction, rank decisions (RCOND, SVLMAX), QR/RQ with pivoting. Map LAPACK DORMQR/DORMRQ to nalgebra Q application.
2. Port AB08NW: build system pencil, call AB08NY, then DTZRZF/DORMRZ (RQ on block), extract Af, Ef.
3. AB08NZ: complex version; either port with num_complex or defer.
4. Add ab08/io.rs: parse AB08NW.dat (N, M, P, TOL, EQUIL, A, B, C, D), parse Fortran output (KRONL, NFZ, A reduced, etc.).
5. Add validate_ab08nw (and ny, nz if needed) in tests/fortran_validation.rs.

## AB08NW .dat format

Line 1: title. Line 2: N, M, P, TOL, EQUIL. Then A (N×N), B (N×M), C (P×N), D (P×M).

## Status

- AB07: Done (ab07md, ab07nd full implementation + validation).
- AB08: Pending full port of AB08NY then AB08NW (and AB08NZ if real .dat).
