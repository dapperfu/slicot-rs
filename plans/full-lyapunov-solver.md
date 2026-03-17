# Full Lyapunov Solver

Implement the SB03 Lyapunov stack end-to-end, starting from the core solver and then filling in the remaining adapter routines that currently short-circuit with `INFO=1`. The goal is to make the solver family usable through the existing SB03 entrypoints without changing the broader AB09 flow.

## Approach

- Keep the work centered on `src/sb03/`, using the already-implemented kernels as building blocks.
- Preserve the current AB09 callers that already rely on `sb03ou`, `sb03oy`, and related helpers.
- Add or update validation so the new solver paths are exercised through the SB03 docs and tests.

## Execution Plan

1. **Persist and track the plan**
   - Ensure the plan is saved under `plans/` and referenced in the implementation work so the SB03 work stays traceable.

2. **Finish the core Lyapunov solver**
   - Review and complete `src/sb03/sb03md.rs` so it fully handles continuous/discrete Lyapunov solves, including the separation/error estimate path expected by the existing API.
   - Verify the low-level kernels it depends on are wired correctly: `src/sb03/sb03mx.rs`, `src/sb03/sb03my.rs`, `src/sb03/sb03mv.rs`, `src/sb03/sb03mw.rs`, and `src/sb03/sb03or.rs`.

3. **Complete the factorized/Gramian path**
   - Finish or tighten the Cholesky-factor Lyapunov routines in `src/sb03/sb03ot.rs`, `src/sb03/sb03ou.rs`, and `src/sb03/sb03oy.rs`.
   - Keep the AB09 consumers working, especially the paths in `src/ab09/ab09ax.rs`, `src/ab09/ab09bx.rs`, `src/ab09/ab09hy.rs`, and `src/ab09/ab09iy.rs`.

4. **Replace the remaining `INFO=1` SB03 adapters**
   - Implement the remaining SB03 adapter files that are currently placeholders: `src/sb03/sb03od.rs`, `src/sb03/sb03os.rs`, `src/sb03/sb03oz.rs`, `src/sb03/sb03pd.rs`, `src/sb03/sb03qd.rs`, `src/sb03/sb03qx.rs`, `src/sb03/sb03qy.rs`, `src/sb03/sb03rd.rs`, `src/sb03/sb03sd.rs`, `src/sb03/sb03sx.rs`, `src/sb03/sb03sy.rs`, `src/sb03/sb03td.rs`, `src/sb03/sb03ud.rs`, and `src/sb03/sb03mu.rs`.
   - Prefer thin forwarding wrappers where possible, but replace placeholder returns with real solver logic when the adapter is the only missing layer.

5. **Validate and document**
   - Update `validation/sb03.md` to move the SB03 entries from "adapter not implemented" to validated once the Rust paths are working.
   - Update `docs/SLICOT_MAPPING.md` if any status changes are needed.
   - Run the SB03 test slice and the AB09 tests that rely on the solver stack.

## Success Criteria

- `sb03md` and the factorized SB03 routines return real results for supported inputs instead of placeholder failures.
- The remaining SB03 adapter routines no longer return `1` for ordinary inputs.
- Existing AB09 callers continue to pass their tests without regressions.
- SB03 validation docs reflect the implemented status.
