# MB02 remaining routines

- **All 35 MB02 routines are now implemented** (dense fallbacks where applicable).
- **Done:** MB02CD, ED (block Toeplitz); MB02SD, RD, TD (Hessenberg LU); MB02FD, GD (Cholesky); MB02CU (in-place Cholesky); MB02CV (copy); MB02CX, CY, DD, HD, JD, JX, KD, MD, ND, OD, PD, QD, UD, UU, UV, UW, VD, WD, XD, YD (solve A*X=B or op(A)*X=B); MB02NY (1-norm), MB02QY (Frobenius norm); MB02ID, RZ, SZ, TZ (in-place LU).
- Dense fallbacks are correct; fast Schur algorithm can be added later for performance.
