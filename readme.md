# Hermes-3 performance test cases

## Test 1
Baseline is a simplified, steady state low power solution based on the ST40 tokamak (see Kryjak @ APS 2023). 
This test restarts from the baseline with x1.5 power and runs 10 output timesteps of 0.005ms each.

- Runs for 10 output timesteps of 0.2ms each.
- Needs 10 cores to run.

Original case ID: p2d3ab-power_x1.5
Original grid ID: g3e4-lores_widev2_nonortho_xpoint.nc

Known performance on M. Kryjak's machine:
 - PID timestepper + Hypre ILU: 6m 19s

## Test 2
Baseline is a full, unsimplified version of Test 1. It's a lot more computationally intensive than Test 1 
and is on the critical path for the SOLPS comparison project. The test restarts from a steady state solution
with a lower neutral pump albedo, which leads to the reduction of plasma density.

- Runs for 10 output timesteps of 0.01ms each. (20x less than Test 1)
- Needs 10 cores to run.

Original case ID: upst1ad-malamas_settings_tune
Original grid ID: g3e4f1-lores_widev2_nonortho_xpoint_allf.nc

Known performance on M. Kryjak's machine:
 - PID timestepper + Hypre ILU: 4m 51s

## Test 3
Based on DIII-D - relevant to M. Tsagkiridis' project. Needs 10 cores to run.

# Running tests
Hermes-3 needs the input (BOUT.inp) and restart files (BOUT.restart.*.nc) to run each test. The restart files are provided in the `base` subdirectory of each test. To reset the test, just copy the restart files into the case directory. Grid files need to be in the same directory that you run Hermes-3 in, so you can run the tests from the root dir.

Run command:

mpirun -np 10 hermes_dir/hermes-3 -d test1 restart


## PETSc configuration
To enable STRUMPACK, use the following configure flags for PETSc:

```
./configure --with-mpi=yes --download-hypre --download-make --with-fortran-bindings=0 --with-debugging=0 --download-strumpack --download-metis --download-parmetis --download-ptscotch --download-zfp --download-scalapack
```

# Solver settings

Standard CVODE settings:

```
[solver]
mxstep = 1e9
cvode_max_order = 3
maxl = 5
atol = 1e-12 * 1
rtol = 1e-6 * 1
use_precon = True
diagnose = false
```

Settings using the latest Hypre ILU (needs PETSc >=3.23.3) with Malamas Tsagkiridis' PID timestepper:

```
[solver]
type = snes                      # Backward Euler steady-state solver
snes_type = newtonls             # Nonlinear solver
ksp_type = gmres                 # Linear solver: gmres, cg
max_nonlinear_iterations = 16    # default: 50
pc_type = hypre                  # Preconditioner type
pc_hypre_type = ilu         
lag_jacobian = 7                 # Iterations between jacobian recalculations. default: 50
atol = 1e-12                      # Absolute tolerance
rtol = 1e-6                      # Relative tolerance
stol = 1e-12
maxf = 20000
maxl = 260
pidController = true
target_its = 5
kP = 0.65
kI = 0.30
kD = 0.15
matrix_free_operator = true
timestep = 0.001                 # Initial timestep

[petsc]

#log_view = true
                                    
pc_hypre_ilu_level = 1                            # k = 2  (default is 0, try 1 and 2)
pc_hypre_ilu_local_reordering = true              # reduces fill / improves robustness
pc_hypre_ilu_tri_solve = true                     # use triangular solve instead of smoothing
pc_hypre_ilu_print_level = true
snes_fd_color_use_mat = true
```

Example STRUMPACK settings with STRUMPACK as direct solver:

```
[solver]
type = snes                     # Backward Euler steady-state solver
snes_type = newtonls            # Nonlinear solver
ksp_type = preonly              # Linear solver
use_precon = false
max_nonlinear_iterations = 15   # default: 50
pc_type = lu                    # Preconditioner type
lag_jacobian = 500              # Iterations between jacobian recalculations. default: 50
atol = 1e-12                    # Absolute tolerance
rtol = 1e-6                     # Relative tolerance
stol = 1e-12
maxl = 20                       # default: 20
use_coloring = true
matrix_free_operator = true

[petsc]
pc_factor_mat_solver_type = strumpack
mat_strumpack_verbose = true
```

