# Analysis: Unimplemented and stubbed functions (100% Rust implementation)

## Sources of truth

- **Canonical list**: [docs/SLICOT_MAPPING.md](docs/SLICOT_MAPPING.md) — column **Implementation**: `stub` | `partial` | `full`.
- **In-code**: Routines that return `1` (INFO=1) for "not yet implemented" or document "not implemented" in module comments.

No `unimplemented!()` or `todo!()` macros appear in the codebase; stubs are implemented as input validation + `return 1` for the main path.

---

## 1. Stub routines (Implementation = stub)

These have a Rust API and trivial/minimal behavior (e.g. only N=0/M=0 or argument checks), and return **INFO=1** for the real algorithm. Each must be replaced by a full algorithm port from SLICOT Fortran.

| Module   | Routines                                                                                                                                                                               | Count |
| -------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----- |
| **ab01** | ab01nd, ab01od                                                                                                                                                                         | 2     |
| **ab08** | ab08nz                                                                                                                                                                                 | 1     |
| **ab8n** | ab8nxz                                                                                                                                                                                 | 1     |
| **ag08** | ag08bd, ag08by, ag08bz                                                                                                                                                                 | 3     |
| **ag8b** | ag8byz                                                                                                                                                                                 | 1     |
| **fd01** | fd01ad                                                                                                                                                                                 | 1     |
| **ib01** | ib01ad, ib01bd, ib01cd, ib01md, ib01my, ib01nd, ib01od, ib01oy, ib01pd, ib01px, ib01py, ib01qd, ib01rd                                                                                 | 13    |
| **ib03** | ib03ad, ib03bd                                                                                                                                                                         | 2     |
| **mc03** | mc03ny                                                                                                                                                                                 | 1     |
| **md03** | md03ad, md03ba, md03bb, md03bd, md03bf, md03bx, md03by                                                                                                                                 | 7     |
| **nf01** | nf01ba, nf01bb, nf01bd, nf01be, nf01bf, nf01bp, nf01bq, nf01br, nf01bs, nf01bu, nf01bv, nf01bw, nf01bx, nf01by                                                                         | 14    |
| **sb04** | sb04md, sb04mr, sb04mu, sb04mw, sb04my, sb04nd, sb04nv, sb04nw, sb04nx, sb04ny, sb04od, sb04ow, sb04pd, sb04py, sb04qd, sb04qr, sb04qu, sb04qy, sb04rd, sb04rv, sb04rw, sb04rx, sb04ry | 23    |
| **sb06** | sb06nd                                                                                                                                                                                 | 1     |
| **sb08** | sb08cd, sb08dd, sb08ed, sb08fd, sb08gd, sb08hd, sb08md, sb08my, sb08nd, sb08ny                                                                                                         | 10    |
| **sb10** | sb10jd, sb10zp                                                                                                                                                                         | 2     |
| **tg01** | tg01dd, tg01wd                                                                                                                                                                         | 2     |

**Total stub routines: 84.**

---

## 2. Partial routines (Implementation = partial)

Main path is implemented; some options or code paths return INFO=1 or are documented as "not implemented". Those branches must be implemented so the routine is 100% complete.

| SLICOT | Module | Rust function | Missing / partial behavior                                                                   |
| ------ | ------ | ------------- | -------------------------------------------------------------------------------------------- |
| AB08MD | ab08   | ab08md        | Some paths/options partial                                                                   |
| AB08MZ | ab08   | ab08mz        | Some paths/options partial                                                                   |
| AB09HY | ab09   | ab09hy        | D-based stochastic path (Riccati, Cw) not implemented; D must be full row rank               |
| AB09KD | ab09   | ab09kd        | With weighting: AB07ND (inverses V,W) not implemented                                        |
| AB09KX | ab09   | ab09kx        | WEIGHT='N' no-op; L/R/B: TB01WD on weights then return (projection formulas not implemented) |
| AB13AD | ab13   | ab13ad        | Full additive decomposition not implemented; equilibration not implemented                   |
| AB13HD | ab13   | ab13hd        | Some options return 1 (not implemented)                                                      |
| AB13MD | ab13   | ab13md        | Some options return 1 (not implemented)                                                      |
| AG07BD | ag07   | ag07bd        | Algorithm not implemented (returns 1)                                                        |
| SB02OD | sb02   | sb02od        | Some options partial                                                                         |
| TB05AD | tb05   | tb05ad        | `evre`/`evim` (eigenvalues) not implemented                                                  |

**Total partial routines: 11.**

---

## 3. In-code branches that return 1 (not yet implemented)

- **ab05/ab05sd.rs**: two branches return 1.
- **ab08/ab08nw.rs**, **ab08/ab08ny.rs**: Full path not yet implemented.
- **ab09** (multiple): various option branches return 1.
- **ab13**: ab13ad, ab13dx, ab13hd, ab13md — unsupported options return 1.
- **ag07/ag07bd.rs**: algorithm not implemented.
- **ib03**: Full LM identification / MINPACK-like LM not yet implemented.
- **dg01/dg01od.rs**, **df01/df01md.rs**, **fb01/fb01qd.rs**: Full transform/covariance update not yet implemented.
- **tb05/tb05ad.rs**: eigenvalues output not implemented.

Routines that document **valid** INFO=1 (e.g. SB04PX "equation perturbed") are **not** stubs.

---

## 4. Summary and requirement

- **Stub (mapping)**: 84 — Replace with full algorithm from SLICOT Fortran.
- **Partial (mapping)**: 11 — Implement missing options/branches.
- **In-code return 1**: Implement missing branches or document as valid SLICOT codes.

Every function SHALL be implemented completely; there SHALL be no stubs.

---

## 5. Implementation approach

- **Order**: Follow [plans/slicot-resolution-order.md](plans/slicot-resolution-order.md).
- **Reference**: Port from `SLICOT-Reference/src/<ROUTINE>.f`.
- **Tracking**: Update SLICOT_MAPPING Implementation from `stub`/`partial` → `full` when done.

Existing plans: [plans/fully-implement-all-slicot-stubs.md](plans/fully-implement-all-slicot-stubs.md), [plans/full-implementation-requirement.md](plans/full-implementation-requirement.md).
