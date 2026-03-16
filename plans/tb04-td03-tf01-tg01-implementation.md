# TB04, TD03AY, TF01, TG01 full implementation plan

## Scope (45 routines)

| Module | Routines | Action |
|--------|----------|--------|
| **TB04** | TB04AD, TB04AY, TB04BD, TB04BV, TB04BW, TB04BX, TB04CD | New module `src/tb04/` |
| **TD03** | TD03AY | Add to existing `src/td03/` |
| **TF01** | TF01MX, TF01MY, TF01ND, TF01OD, TF01PD, TF01QD, TF01RD | Add to existing `src/tf01/` |
| **TG01** | 30 routines (TG01AD through TG01WD) | New module `src/tg01/` |

**Constraints:** Full implementations only (no stubs). Every routine must have at least one test that checks meaningful numerical behavior.

## Conventions

- One file per routine: `src/<module>/<name>.rs` (lowercase).
- Module: `mod.rs` with `pub mod <name>;`.
- lib.rs: add `pub mod tb04;` and `pub mod tg01;`.
- Tests: `#[cfg(test)] mod tests` with at least one `#[test]` asserting on outputs.
- After each routine: set status to `done` in docs/SLICOT_MAPPING.md; check off in plans/remaining-slicot-functions-todo.md.

## Execution order

1. Save plan to `plans/`.
2. TB04: Create `src/tb04/`, implement TB04AD then TB04AY, TB04BD, TB04BV, TB04BW, TB04BX, TB04CD.
3. TD03: Implement TD03AY.
4. TF01: Implement TF01MX, TF01MY, TF01ND, TF01OD, TF01PD, TF01QD, TF01RD.
5. TG01: Create `src/tg01/`, implement all 30 routines in dependency order.
6. Update SLICOT_MAPPING.md, remaining-slicot-functions-todo.md; run validation.
