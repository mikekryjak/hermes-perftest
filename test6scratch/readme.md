# test6: ST40 low res test

Low-resolution ST40 at neutral flux limit 1.0, with separate limiters, run
from scratch over 0-100 ms (100 outputs of 1 ms). It is the low-resolution
counterpart of test4, which runs full-resolution ST40 at flux limit 0.2 with
combined limiters.

Flux limit 1.0 is used because it is the most robust setting on this grid:
runs at 0.2 have stalled or run slowly when regularisation, variable scaling
or limiter lagging were missing, so a ladder of settings at 0.2 would not
finish. 1.0 is also closer to production.

## Source

`BOUT.inp` is the input of seplim series 24 rung d,
`/home/mike/work/cases/fluxlim/seplim/splim24d-st40lr_seplim-unifloor_nocombined_floor10_ratio40_flim1.0_scratch`:
separate limiters (`d:combined_limiters = false`), perpendicular floor 10,
ratio 40, ceiling 100, Van Albada. The header of `BOUT.inp` lists every change
from it:

- the window is 0-100 ms (`nout = 100`, `timestep = 95788`);
- the neutral density and pressure start non-zero (1e16 m^-3 at 3 eV, as in
  test4), so `solver:scale_vars` does not floor them;
- `input:error_on_unused_options = true`;
- the grid is `grid_test6.nc` in the repo root, a copy of
  `g3e4f1-lores_widev2_nonortho_xpoint_allf.nc`.
- the stale `[mesh]` profile-geometry lines (`nx = 40`, `MXG`, `xsep_inner`,
  `x_index`) are removed: the grid is 20 cells wide and nothing read them.

The source ran on a separate-limiter worktree build. Some limiter option
names differ between builds; `error_on_unused_options` names any option a
build does not read.

## Start state

There are no restart files. Every run starts from the analytic profiles in
`BOUT.inp` (`initial_ne`, `initial_pi`, `initial_nd`, `initial_pd`), so run it
without `restart`.

## Ladder

Each row adds one change to the row above it, as on test4. Settings on top of
`BOUT.inp` as it stands:

| row | settings |
|---|---|
| SNES-MUMPS-1 | `BOUT.inp` as it stands (SNES-MUMPS-1 solver block, `target_its = 100`); to be checked against the recipe file |
| + limiter lagging | `d:lag_limiter = gradient` |
| + SNES-MUMPS-3 | recipe `recipes/SNES-MUMPS-3.txt` |
| + SNES-MUMPS-4 | recipe `recipes/SNES-MUMPS-4.txt` |
| + Van Albada tolerance | `fv:limiter_epsilon = 1e-3`, `fv:limiter_epsilon_clip = true` |
| + Jacobian colouring | `solver:force_symmetric_coloring = true`, `petsc:mat_coloring_type = id`, `solver:stencil:auto_pairs = true` |
| - limiter lagging | `d:lag_limiter = off` |

Rows before SNES-MUMPS-1 (CVODE-1, and runs without regularisation or with
the MC limiter) are decided after the rows above have run.
