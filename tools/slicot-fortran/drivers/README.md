# SLICOT example drivers with file I/O

These drivers are copies of the SLICOT example programs adapted to accept **optional file arguments** so the fuzzer can pass input/output paths.

**Convention:** `T* [input.dat [output.res]]`

- No args: read from stdin (5), write to stdout (6), same as original.
- One arg: read from `input.dat`, write to stdout.
- Two args: read from `input.dat`, write to `output.res`.

## TAB01ND

Copy into your SLICOT-Reference tree to enable file-based I/O for AB01ND:

```bash
cp tools/slicot-fortran/drivers/TAB01ND.f SLICOT-Reference/examples/TAB01ND.f
```

Then rebuild the examples. After that you can run e.g.:

```bash
cd SLICOT-Reference/examples
./TAB01ND ../data/AB01ND.dat out.res
```

**.dat layout (AB01ND):** Optional title line; then one line with N, M, TOL, JOBZ; then N lines of A (N values per line); then N lines of B (M values per line, Fortran order A(I,J) I=1,N for each J).

**.res layout (AB01ND):** Title line; INFO line if non-zero; else NCONT, block of A(1:NCONT,1:NCONT), NBLK(1:INDCON), block of B(1:NCONT,1:M), INDCON, optionally Z matrix.
