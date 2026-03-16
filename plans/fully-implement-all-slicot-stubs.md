# Fully implement every SLICOT stub and validate against FORTRAN

## Decisions (from questions)

- **Scope**: All stubs including NF01 (16), MD03 (7); plus complete partial implementations (TB01MD lower path, MB01QD triangular/Hessenberg options). Search codebase for any other "not yet implemented" branches and include in plan.
- **Order**: By module (e.g. finish all AB07, then AB08, then AB09, …).
- **Validation**: Dedicated `validate_*` in tests/fortran_validation.rs for each routine that has `.dat` in SLICOT-Reference/examples/data; full compare INFO and outputs (like AB01MD/AB01ND).
- **I/O**: Split by module: new parsers in module-specific files (e.g. src/ab09/io.rs) rather than a single slicot_io.rs.
- **Tolerance**: Allow per-routine tolerance where documented (default REL_TOL=1e-10; relax for iterative routines if needed).
- **Dependencies**: May introduce new crates (e.g. rustlapack) if needed for correctness; prefer existing nalgebra + crate helpers where sufficient.
- **Status**: Add status **validated** in docs/SLICOT_MAPPING.md for routines that pass FORTRAN validation; keep **done** for implemented but not yet validated.
- **Tracking**: Reuse plans/remaining-slicot-functions-todo.md with a "validated" column or sub-items; no new tracking file.
- **Fortran build**: Assume contributors run ./tools/slicot-fortran/build_fortran.sh when SLICOT-Reference is present; document in plan.

## Stub inventory (routines to fully implement)

| Module | Routines | Reference .dat |
|--------|----------|----------------|
| **ab07** | ab07md, ab07nd | AB07MD.dat, AB07ND.dat |
| **ab08** | ab08nw, ab08ny, ab08nz | AB08ND.dat, AB08NW.dat, AB08NZ.dat |
| **ab09** | ab09ad … ab09nd (24 routines) | Many |
| **ab13** | ab13ed, ab13fd, ab13hd, ab13md | AB13*.dat |
| **ab8n** | ab8nxz | — |
| **ag07** | ag07bd | — |
| **ag08** | ag08bd, ag08by, ag08bz | AG08BD.dat, AG08BZ.dat |
| **ag8b** | ag8byz | — |
| **fb01** | fb01qd, fb01rd, fb01sd, fb01td, fb01vd | FB01*.dat |
| **fd01** | fd01ad | FD01AD.dat |
| **ib01** | ib01ad … ib01rd (13 routines) | IB01AD, IB01BD, IB01CD.dat |
| **ib03** | ib03ad, ib03bd | IB03AD.dat, IB03BD.dat |
| **mc03** | mc03nd, mc03ny | MC03ND.dat |
| **md03** | md03ad … md03by (7) | MD03AD.dat, MD03BD.dat |
| **nf01** | all nf01* (16) | — |
| **Partial** | TB01MD (lower path), MB01QD (triangular/Hessenberg) | — |

## Additional partial branches (from codebase search)

- **TB01MD** (`src/tb01/tb01md.rs`): returns -2 for "Lower not yet implemented in pilot".
- **MB01QD** (`src/mb01/mb01qd.rs`): "Lower triangular", "Upper triangular", "Upper Hessenberg" options not yet implemented.
- **FB01*** (fb01qd, fb01rd, fb01sd, fb01td, fb01vd): SLICOT stubs "not yet implemented".
- **AB07MD/AB07ND**: "Validated stub: returns 0 when N=0 and M=0; 1 (not yet implemented) otherwise."

(MA02RD uses partial_cmp for sorting; MB02SD comment is descriptive, not a stub.)

## Implementation pattern (per routine)

1. Port algorithm from SLICOT-Reference/src/<ROUTINE>.f to Rust.
2. Add module I/O parsers if .dat exists.
3. Add dedicated validate_* in tests/fortran_validation.rs.
4. Unit tests; set status validated in SLICOT_MAPPING when FORTRAN passes.

## Module order

1. AB07 2. AB08 3. AB09 4. AB13 5. AB8N 6. AG07 7. AG08 8. AG8B 9. FB01 10. FD01 11. IB01 12. IB03 13. MC03 14. MD03 15. NF01 16. Partial (TB01MD, MB01QD)

## Completion criteria

- Every stub and partial implementation fully implemented; every routine with .dat has validate_* and passes; SLICOT_MAPPING uses **validated** where applicable.
