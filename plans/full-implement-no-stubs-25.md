# Full implementation of 25 stubbed routines (no main path left undone)

## Requirement
Fully implement every function listed; DO NOT leave main path undone.

## Routines (from terminal selection)

### FB01 (5)
- [ ] FB01QD
- [ ] FB01RD
- [ ] FB01SD
- [ ] FB01TD
- [ ] FB01VD

### AB09 (20)
- [ ] AB09IX, AB09HD, AB09JV, AB09KX, AB09IY, AB09GD
- [ ] AB09CX, AB09CD
- [ ] AB09JD, AB09ED, AB09HX, AB09HY, AB09JW
- [x] AB09DD
- [ ] AB09KD, AB09FD
- [ ] AB09JX, AB09ND, AB09ID, AB09MD
- [x] AB09BX, AB09BD

## Dependency order (AB09)
1. **AB09DD** — SPA formulas. **DONE.**
2. **AB09BX** — SPA (calls AB09DD, SB03OU, MB03UD). **DONE.**
3. **AB09BD** — SPA with D (calls AB09BX, TB01ID, TB01WD). **DONE.**
4. AB09CD — (calls AB09DD or similar).
5. AB09MD — (calls AB09AX, TB01KD).
6. AB09ND — (calls AB09BX, TB01KD).
7. AB09FD, AB09GD — (call AB09AX/AB09BX, SB08*).
8. AB09ID, AB09JD, AB09KD, AB09ED, AB09HD — (various cores).
9. AB09CX, AB09HX, AB09IX, AB09HY, AB09IY, AB09JV, AB09JW, AB09JX, AB09KX.

## Status
- AB09DD: full port.
- AB09BX: full port (ab09bx_core, ab09bx_full, ab09bx(n,m)); SB03OU×2, MB03UD, NMINR, AB09DD.
- AB09BD: full port (TB01ID, TB01WD, AB09BX with D).
