# SLICOT 100% feature parity roadmap

This plan lists all **475** unimplemented SLICOT routines by phase and module. When all are implemented in pure Rust and marked `done` in [docs/SLICOT_MAPPING.md](../docs/SLICOT_MAPPING.md), the crate will have 100% feature coverage of the mapping.

Feature status summary: [docs/FEATURES.md](../docs/FEATURES.md). Regenerate with `./tools/gen_features_table.sh` after mapping changes.

## Phases (by dependency and usage)

| Phase | Module(s) | Routines | Description |
|-------|-----------|----------|-------------|
| 1 | ab13, ma01, ma02, mb01 | 20 | Complete partially-done modules (norms, basic matrix, BLAS-like) |
| 2 | mb02 | 35 | Matrix operations (structure-exploiting linear algebra) |
| 3 | mb03, mb04, mb05, mb3j, mb3l, mb3o, mb3p, mb4d | 161 | Eigenvalue/Schur, Hessenberg, and related factorizations |
| 4 | mc01, mc03, md03 | 26 | Polynomial operations, matrix pencils, data fitting |
| 5 | nf01 | 16 | Nonlinear optimization (NF01) |
| 6 | sb01, sb02, sb03, sb04, sb06, sb08, sb09, sb10, sb16 | 113 | Riccati equations, LQR, LQG, H-infinity, etc. |
| 7 | sg02, sg03 | 18 | Generalized Riccati (discrete/generalized) |
| 8 | tb01, tb03, tb04, tb05, tc01, tc04, tc05, td03, td04, td05, tf01 | 46 | System transformations, frequency response, polynomial forms |
| 9 | tg01 | 30 | Descriptor system transformations |
| 10 | ud01, ue01, zgeg, zlat | 10 | Utility (UD/UE), complex generalized eigen (ZGEG), complex ZLATZM |

**Total unimplemented: 475.**

## Per-routine checklist (by module)

Check off each when the routine is implemented in pure Rust, has at least one test, and status in `docs/SLICOT_MAPPING.md` is set to `done`. Run `./tools/validate_slicot_done.sh` after each batch.

### Phase 1: Complete partial modules (20)

#### ab13 (6)
- [x] AB13AD
- [x] AB13AX
- [x] AB13BD
- [x] AB13CD
- [x] AB13DX
- [x] AB13ID

#### ma01 (1)
- [x] MA01CD

#### ma02 (11)
- [x] MA02HD
- [x] MA02HZ
- [x] MA02ID
- [x] MA02IZ
- [x] MA02JD
- [x] MA02JZ
- [x] MA02MD
- [x] MA02MZ
- [x] MA02OD
- [x] MA02OZ
- [x] MA02SD

#### mb01 (2)
- [x] MB01YD
- [x] MB01ZD

### Phase 2: mb02 (35)

- [x] MB02CD
- [ ] MB02CU
- [ ] MB02CV
- [ ] MB02CX
- [ ] MB02CY
- [ ] MB02DD
- [x] MB02ED
- [ ] MB02FD
- [ ] MB02GD
- [ ] MB02HD
- [ ] MB02ID
- [ ] MB02JD
- [ ] MB02JX
- [ ] MB02KD
- [ ] MB02MD
- [ ] MB02ND
- [ ] MB02NY
- [ ] MB02OD
- [ ] MB02PD
- [ ] MB02QD
- [ ] MB02QY
- [ ] MB02RD
- [ ] MB02RZ
- [ ] MB02SD
- [ ] MB02SZ
- [ ] MB02TD
- [ ] MB02TZ
- [ ] MB02UD
- [ ] MB02UU
- [ ] MB02UV
- [ ] MB02UW
- [ ] MB02VD
- [ ] MB02WD
- [ ] MB02XD
- [ ] MB02YD

### Phase 3: mb03, mb04, mb05, mb3j, mb3l, mb3o, mb3p, mb4d (160)

#### mb03 (79)
- [ ] MB03AB MB03AD MB03AE MB03AF MB03AG MB03AH MB03AI
- [ ] MB03BA MB03BB MB03BC MB03BD MB03BE MB03BF MB03BG MB03BZ
- [ ] MB03CD MB03CZ MB03DD MB03DZ MB03ED MB03FD MB03FZ MB03GD MB03GZ
- [ ] MB03HD MB03HZ MB03ID MB03IZ MB03JD MB03JP MB03JZ
- [ ] MB03KA MB03KB MB03KC MB03KD MB03KE
- [ ] MB03LD MB03LF MB03LP MB03LZ MB03MD MB03MY MB03ND MB03NY
- [ ] MB03OD MB03OY MB03PD MB03PY MB03QD MB03QG MB03QV MB03QW MB03QX MB03QY
- [ ] MB03RD MB03RW MB03RX MB03RY MB03RZ MB03SD MB03TD MB03TS
- [ ] MB03UD MB03VD MB03VW MB03VY MB03WA MB03WD MB03WX
- [ ] MB03XD MB03XP MB03XS MB03XU MB03XZ MB03YA MB03YD MB03YT MB03ZA MB03ZD

#### mb04 (70)
- [ ] MB04AD MB04AZ MB04BD MB04BP MB04BZ MB04CD MB04DB MB04DD MB04DI MB04DL MB04DP MB04DS MB04DY MB04DZ
- [ ] MB04ED MB04FD MB04FP MB04GD MB04HD MB04ID MB04IY MB04IZ MB04JD MB04KD MB04LD MB04MD MB04ND MB04NY
- [ ] MB04OD MB04OW MB04OX MB04OY MB04PA MB04PB MB04PU MB04PY
- [ ] MB04QB MB04QC MB04QF MB04QS MB04QU MB04RB MB04RD MB04RS MB04RT MB04RU MB04RV MB04RW MB04RZ
- [ ] MB04SU MB04TB MB04TS MB04TT MB04TU MB04TV MB04TW MB04TX MB04TY
- [ ] MB04UD MB04VD MB04VX MB04WD MB04WP MB04WR MB04WU MB04XD MB04XY MB04YD MB04YW MB04ZD

#### mb05 (5)
- [ ] MB05MD MB05MY MB05ND MB05OD MB05OY

#### mb3j, mb3l, mb3o, mb3p (4)
- [ ] MB3JZP MB3LZP MB3OYZ MB3PYZ

#### mb4d (3)
- [ ] MB4DBZ MB4DLZ MB4DPZ

### Phase 4: mc01, mc03, md03 (26)

#### mc01 (15)
- [ ] MC01MD MC01ND MC01OD MC01PD MC01PY MC01QD MC01RD MC01SD MC01SW MC01SX MC01SY MC01TD MC01VD MC01WD MC01XD

#### mc03 (4)
- [ ] MC03MD MC03ND MC03NX MC03NY

#### md03 (7)
- [ ] MD03AD MD03BA MD03BB MD03BD MD03BF MD03BX MD03BY

### Phase 5: nf01 (16)

- [ ] NF01AD NF01AY NF01BA NF01BB NF01BD NF01BE NF01BF NF01BP NF01BQ NF01BR NF01BS NF01BU NF01BV NF01BW NF01BX NF01BY

### Phase 6: sb01–sb16 (113)

#### sb01 (6)
- [ ] SB01BD SB01BX SB01BY SB01DD SB01FY SB01MD

#### sb02 (21)
- [ ] SB02CX SB02MD SB02MR SB02MS SB02MT SB02MU SB02MV SB02MW SB02MX SB02ND SB02OD SB02OU SB02OV SB02OW SB02OX SB02OY SB02PD SB02QD SB02RD SB02RU SB02SD

#### sb03 (24)
- [ ] SB03MD SB03MU SB03MV SB03MW SB03MX SB03MY SB03OD SB03OR SB03OS SB03OT SB03OU SB03OV SB03OY SB03OZ SB03PD SB03QD SB03QX SB03QY SB03RD SB03SD SB03SX SB03SY SB03TD SB03UD

#### sb04 (24)
- [ ] SB04MD SB04MR SB04MU SB04MW SB04MY SB04ND SB04NV SB04NW SB04NX SB04NY SB04OD SB04OW SB04PD SB04PX SB04PY SB04QD SB04QR SB04QU SB04QY SB04RD SB04RV SB04RW SB04RX SB04RY

#### sb06, sb08, sb09 (12)
- [ ] SB06ND
- [ ] SB08CD SB08DD SB08ED SB08FD SB08GD SB08HD SB08MD SB08MY SB08ND SB08NY
- [ ] SB09MD

#### sb10 (21)
- [ ] SB10AD SB10DD SB10ED SB10FD SB10HD SB10ID SB10JD SB10KD SB10LD SB10MD SB10PD SB10QD SB10RD SB10SD SB10TD SB10UD SB10VD SB10WD SB10YD SB10ZD SB10ZP

#### sb16 (5)
- [ ] SB16AD SB16AY SB16BD SB16CD SB16CY

### Phase 7: sg02, sg03 (18)

#### sg02 (5)
- [ ] SG02AD SG02CV SG02CW SG02CX SG02ND

#### sg03 (13)
- [ ] SG03AD SG03AX SG03AY SG03BD SG03BR SG03BS SG03BT SG03BU SG03BV SG03BW SG03BX SG03BY SG03BZ

### Phase 8: tb01, tb03, tb04, tb05, tc01, tc04, tc05, td03, td04, td05, tf01 (52)

#### tb01 (21)
- [ ] TB01ID TB01IZ TB01KD TB01KX TB01LD TB01ND TB01PD TB01PX TB01TD TB01TY TB01UD TB01UX TB01UY TB01VD TB01VY TB01WD TB01WX TB01XD TB01XZ TB01YD TB01ZD

#### tb03, tb04, tb05 (10)
- [ ] TB03AD TB03AY
- [ ] TB04AD TB04AY TB04BD TB04BV TB04BW TB04BX TB04CD
- [ ] TB05AD

#### tc01, tc04, tc05 (3)
- [ ] TC01OD TC04AD TC05AD

#### td03, td04, td05 (4)
- [ ] TD03AD TD03AY TD04AD TD05AD

#### tf01 (8)
- [ ] TF01MD TF01MX TF01MY TF01ND TF01OD TF01PD TF01QD TF01RD

### Phase 9: tg01 (30)

- [ ] TG01AD TG01AZ TG01BD TG01CD TG01DD TG01ED TG01FD TG01FZ TG01GD TG01HD TG01HU TG01HX TG01HY TG01ID TG01JD TG01JY TG01KD TG01KZ TG01LD TG01LY TG01MD TG01ND TG01NX TG01OA TG01OB TG01OD TG01OZ TG01PD TG01QD TG01WD

### Phase 10: ud01, ue01, zgeg, zlat (10)

#### ud01 (6)
- [ ] UD01BD UD01CD UD01DD UD01MD UD01MZ UD01ND

#### ue01 (1)
- [ ] UE01MD

#### zgeg (2)
- [ ] ZGEGS ZGEGV

#### zlat (1)
- [ ] ZLATZM

---

## Workflow per routine

1. Add Rust module/function under `src/<module>/<function>.rs` matching [docs/SLICOT_MAPPING.md](../docs/SLICOT_MAPPING.md).
2. Implement in pure Rust (no FFI).
3. Add at least one `#[cfg(test)]` with `#[test]` in the implementation file.
4. Set status to `done` in docs/SLICOT_MAPPING.md.
5. Run `./tools/validate_slicot_done.sh` and fix any failures.
6. Regenerate features table: `./tools/gen_features_table.sh` (optional; update [docs/FEATURES.md](../docs/FEATURES.md)).
7. Check off the routine in this roadmap.

## Completion

100% feature parity is reached when all 625 routines have status `done` in SLICOT_MAPPING.md and this checklist is fully checked. Then [docs/FEATURES.md](../docs/FEATURES.md) will show 625 implemented, 0 unimplemented.
