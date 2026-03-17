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
- [x] AB09CX, AB09CD
- [ ] AB09JD, AB09HX, AB09HY, AB09JW
- [x] AB09ED
- [x] AB09DD
- [ ] AB09KD, AB09FD
- [ ] AB09JX, AB09ID
- [x] AB09ND, AB09MD, AB09BX, AB09BD

## Dependency order (AB09)
1. **AB09DD** — SPA formulas. **DONE.**
2. **AB09BX** — SPA (calls AB09DD, SB03OU, MB03UD). **DONE.**
3. **AB09BD** — SPA with D (calls AB09BX, TB01ID, TB01WD). **DONE.**
4. **AB09CX** — Hankel-norm (AB09AX, TB01WD, TB01KD, AB04MD). **DONE.**
5. **AB09CD** — TB01ID, TB01WD, AB09CX. **DONE.**
6. **AB09MD** — TB01ID, TB01KD, AB09AX. **DONE.**
7. **AB09ND** — TB01ID, TB01KD, AB09BX. **DONE.**
8. **AB09ED** — TB01ID, TB01KD, AB09CX. **DONE.**
9. AB09FD, AB09GD — (call AB09AX/AB09BX, SB08*).
10. AB09ID, AB09JD, AB09KD, AB09HD — (various cores).
11. AB09HX, AB09IX, AB09HY, AB09IY, AB09JV, AB09JW, AB09JX, AB09KX.

## Status
- AB09DD, AB09BX, AB09BD: full port.
- AB09CX: full port (AB09AX, order selection, TB01WD or Hankel step with pinv, MB01SD, TB01KD, AB04MD).
- AB09CD: full port (TB01ID, TB01WD, AB09CX).
- AB09MD, AB09ND: full port (TB01ID, TB01KD, AB09AX/AB09BX on stable part).
- AB09ED: full port (TB01ID, TB01KD, AB09CX on stable part).
- Remaining: AB09HY, IX, HX, HD, JX, JV, JW, JD, KX, KD, IY, ID, FD, GD (14 routines).
