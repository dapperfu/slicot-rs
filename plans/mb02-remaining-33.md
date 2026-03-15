# MB02 remaining 33 routines

After AB13 and MB02CD/ED full implementations:

- **Done:** MB02CD (block Toeplitz Cholesky, dense formation), MB02ED (solve T*X=B or X*T=B, dense).
- **Done (Hessenberg LU):** MB02SD (LU factor), MB02RD (solve), MB02TD (RCOND).
- **Done (dense Cholesky):** MB02FD, MB02GD (s.p.d. Cholesky R, R'*R = A).
- **Remaining 28:** MB02CU, CV, CX, CY, DD, HD, ID, JD, JX, KD, MD, ND, NY, OD, PD, QD, QY, RZ, SZ, TZ, UD, UU, UV, UW, VD, WD, XD, YD.

Each requires either the fast Schur algorithm (Householder + hyperbolic rotations) or a dense fallback (form full matrix, then standard factor/solve). Dense fallbacks are correct but O(N³K³) vs intended O(K³N²).

Next: implement remaining MB02 using dense formation + Cholesky/solve where applicable, or full Schur algorithm for performance.
