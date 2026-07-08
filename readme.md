# Intro
This repo has standardised performance test cases for Hermes-3.

# Running tests
Hermes-3 needs the input (BOUT.inp) and restart files (BOUT.restart.*.nc) to run each test. The restart files are provided in the `base` subdirectory of each test. To reset the test, just copy the restart files into the case directory. Grid files need to be in the same directory that you run Hermes-3 in, so you can run the tests from the root dir.

Run command:

```
mpirun -np 10 hermes_dir/hermes-3 -d test1 restart
```

***Note:***
If you are running on a workstation, you may need to pin the simulation to certain cores, for example:
> taskset -c 0-9 mpirun --bind-to=none -np 10 hermes_dir/hermes-3 -d test1 restart

A wrapper script is provided in `/sdtools/cli/sdrun.py`.

## Post-processing
This repo includes M. Kryjak's personal post-processing script repo `sdtools`. The cases can be post-processed using the `cmonitor.py` tool:

```
cmonitor.py -s -solverdiags test1
```

Outputs from running this on M. Kryjak's machine are included in the repo. In the plots, the top row shows the evolution of physical quantities. The second row shows simulation speed in ms of simulation time per 24hrs of wall time. The third row shows the RMS of the LHS of the equations, i.e. ddt(Ne) etc. These plots show you which quantities are varying the most. The final row is parsed from the log file and contains the SNES diagnostics. These are on a solver timestep basis while the rest of the plot is on an output timestep basis - note that these are not necessarily the same, as the solver timesteps are taken from the log file which is overwritten per run, while the output time comes from the dataset.


# Useful tools
sdtools now has a few different command line tools to help with test running.

## Printing test results
`test_results.py` will print out test speed and some provenance information:

```
> test_results.py test5steady-example
Run statistics: /home/mike/work/cases/perftests/test5dev/test5steady-example  [BOUT.log.0]
  Originator            : mikekryjak@gmail.com   (git email; assumes script runs as the run user)
  Run started           : Mon Jun 15 15:24:13 2026
  Run finished          : Mon Jun 15 15:28:28 2026
  Run time              : 4 m 15 s  (255 s total)
  Hermes-3 commit       : 44ce2e6c5e140820601533937426c141a294e65c
  Hermes-3 branch       : auto-update-pydeps, master   (from local git, not run dir)
  Hermes-3 commit date  : 2026-06-05 15:53:59 -0700   (from local git, not run dir)
  BOUT++ version        : 5.2.1
  BOUT++ commit         : 8218f71a2fe495bc138c7ddf7c33bb1be9837525
  BOUT++ branch         : unknown   (needs --git / local repo)
  BOUT++ commit date    : 2026-05-11 18:31:52 +0100   (from local git, not run dir)
```

## Recording test results in a csv
`record_test.py` will put it in a csv database. By default the csv is made in the current working directory, but you can change it with an argument. You need to pass a branch name for the record as well as the recipe text file. It will record any differences to the recipe in the diffs column.

```
> record_test.py test5steady-example -b master -r $perftest/recipes/CVODE-2.txt
  case               : test5steady-example
  originator         : mikekryjak@gmail.com
  run_started        : Mon Jun 15 15:34:23 2026
  run_time_str       : 5 m 4 s
  run_time_s         : 304
  recipe             : CVODE-2
  diffs              : 
  note               : 
  hermes_branch      : master
  hermes_commit      : 44ce2e6c5e140820601533937426c141a294e65c
  hermes_commit_date : 2026-06-05 15:53:59 -0700
  bout_version       : 5.2.1
  bout_commit        : 8218f71a2fe495bc138c7ddf7c33bb1be9837525
  recorded_at        : 2026-06-24 11:05:25

Record this run as hermes branch 'master'? [y/N] y

Appended to /home/mike/work/cases/perftests/test5dev/run_records.csv
```

## Applying recipe
`apply_recipe.py` will copy a recipe from a text file into an input file:
```
> apply_recipe.py test5steady-example $perftest/recipes/SNES-MUMPS-1.txt
-> Applied recipe 'SNES-MUMPS-1.txt' to test5steady-example/BOUT.inp
   [solver] overwritten (22 lines)
   [petsc] inserted (9 lines)
```


## Resetting test
`reset_test.py` will copy the restart files from the `base` directory in the test directory and print a confirmation of the simulation time before and after:

```
> reset_test.py test5steady-example
WARNING: this will delete 10 dump/log/.pid file(s) from '/home/mike/work/cases/perftests/test5dev/test5steady-example':
  BOUT.dmp.0.nc
  BOUT.dmp.1.nc
  BOUT.dmp.2.nc
  BOUT.dmp.3.nc
  BOUT.dmp.4.nc
  BOUT.dmp.5.nc
  BOUT.dmp.6.nc
  BOUT.dmp.7.nc
  BOUT.dmp.8.nc
  BOUT.dmp.9.nc
Delete these files and reset the case? [yes/no] y
Deleted 10 dump/log/.pid file(s).
Reset /home/mike/work/cases/perftests/test5dev/test5steady-example: copied 10 restart file(s) from /home/mike/work/cases/perftests/test5dev/test5steady-example/base
  Previous simulation time: 59.9948 ms
  New simulation time:      9.99913 ms
```

## Copying over restart files from one test to another
`transplant.py` will copy restart files from one case to another one:

```
> transplant.py test5ref-m2ul3a-from_lin_260519_master_scratch_cvode test5steady-example/
-> Case test5steady-example/ cleaned, files removed: ['BOUT.restart.14.nc', 'BOUT.restart.13.nc', 'BOUT.restart.11.nc', 'BOUT.restart.18.nc', 'BOUT.restart.9.nc', 'BOUT.restart.17.nc', 'BOUT.restart.0.nc', 'BOUT.restart.3.nc', 'BOUT.restart.7.nc', 'BOUT.restart.1.nc', 'BOUT.restart.2.nc', 'BOUT.restart.16.nc', 'BOUT.restart.4.nc', 'BOUT.restart.6.nc', 'BOUT.restart.19.nc', 'BOUT.restart.15.nc', 'BOUT.restart.12.nc', 'BOUT.restart.8.nc', 'BOUT.restart.5.nc', 'BOUT.restart.10.nc']
-> Copied BOUT.restart.14.nc
-> Copied BOUT.restart.13.nc
-> Copied BOUT.restart.11.nc
-> Copied BOUT.dmp.1.nc
-> Copied BOUT.restart.18.nc
-> Copied BOUT.restart.9.nc
-> Copied BOUT.restart.17.nc
-> Copied BOUT.restart.0.nc
-> Copied BOUT.restart.3.nc
-> Copied BOUT.dmp.0.nc
-> Copied BOUT.restart.7.nc
-> Copied BOUT.restart.1.nc
-> Copied BOUT.restart.2.nc
-> Copied BOUT.dmp.6.nc
-> Copied BOUT.dmp.2.nc
-> Copied BOUT.restart.16.nc
-> Copied BOUT.dmp.7.nc
-> Copied BOUT.dmp.4.nc
-> Copied BOUT.dmp.9.nc
-> Copied BOUT.restart.4.nc
-> Copied BOUT.dmp.3.nc
-> Copied BOUT.restart.6.nc
-> Copied BOUT.dmp.8.nc
-> Copied BOUT.restart.19.nc
-> Copied BOUT.restart.15.nc
-> Copied BOUT.restart.12.nc
-> Copied BOUT.restart.8.nc
-> Copied BOUT.restart.5.nc
-> Copied BOUT.dmp.5.nc
-> Copied BOUT.restart.10.nc
Transplant completed from test5ref-m2ul3a-from_lin_260519_master_scratch_cvode to test5steady-example/
```


# Test cases
## Test 1
Based on a simplified, steady state low power solution based on the ST40 tokamak (see Kryjak @ APS 2023). 
This test restarts from the baseline with x1.5 power and runs 10 output timesteps of 0.005ms each.

- Runs for 10 output timesteps of 0.2ms each.
- Needs 10 cores to run.
- Original case ID: p2d3ab-power_x1.5
- Original grid ID: g3e4-lores_widev2_nonortho_xpoint.nc

Known performance on M. Kryjak's machine (see solver settings section):
 - SNES-1 settings: 5m 3s (~550 ms/24hrs) (Aug 2025)
 - SNES-MUMPS-1 settings: 4m 45s (Dec 2025)

## Test 2
Based on the full, unsimplified version of Test 1. It's a lot more computationally intensive than Test 1 
and is on the critical path for the SOLPS comparison project. The test restarts from a steady state solution
with a lower neutral pump albedo, which leads to the reduction of plasma density. A strange and potentially worrying fact is that the convergence reason is nearly always 4, which suggests it's not reaching atol/rtol.

- Runs for 10 output timesteps of 0.01ms each. (20x less than Test 1)
- Needs 10 cores to run.
- Original case ID: upst1ad-malamas_settings_tune
- Original grid ID: g3e4f1-lores_widev2_nonortho_xpoint_allf.nc

Known performance on M. Kryjak's machine (short):
 - SNES-1 settings: 4m 44s (~30 ms/24hrs) (Sept 2025)
 - SNES-STRUMPACK-2 settings: 45s  (Sept 2025)
 - SNES-MUMPS-1 settings: 1m 42s (Dec 2025)

### Test 2long
Same as test 3 but ran for 3ms, which is short enough to be quick to run using good recipes but long enough to capture both the initial transient and a bit of the steady state.

Known performance on M. Kryjak's machine (short):
 - SNES-MUMPS-1 settings: 7m 1s (Dec 2025)

### Test 2scratch
Ran from scratch for 45ms.
Default recipe: SNES-MUMPS-2 

### Test 2stiff
Restarted from test2scratch at 2ms, when the simulation is stiffer.
Default recipe: SNES-MUMPS-2
Runs for 50 timesteps to get 5ms. 

### Test 2steady
Restarted from test2scratch at 10ms, when the simulation is more steady.
Default recipe: SNES-MUMPS-2

## Test 3
Based on DIII-D. Relevant to M. Tsagkiridis' project. It's a very challenging
test because it's nearly from scratch, so everything is changing. Needs 10 cores to run.

Known performance on M. Kryjak's machine:
 - SNES-1 settings: 3m 11s  (~8 ms/24hrs) (Sept 2025)
 - SNES-MUMPS-1 settings: 26s (Dec 2025)

## Test 4scratch
Full resolution ST40 simulation (similar to Test 2, but 4x higher res). 
It's a low density, high temperature case, almost sheath limited.
Original case name: st40fl2-master_scratch
Original grid name: g4bf1-fine_nonorth_xpoint_guards_allfields.nc

Ran for 45ms from scratch.
Default recipe: SNES-MUMPS-2

### Test 4stiff
Restarted from test4scratch at 2ms, when the simulation is very stiff.
Runs for 30 output steps to get 1ms.
Default recipe: SNES-MUMPS-2

### Test 4steady
Restarted from test4scratch at 10ms, when the simulation is more steady.
Runs for 90 output steps to get 0.9ms.
Default recipe: SNES-MUMPS-2

## Test 5scratch
Low resolution MAST-U simulation. Runs OK on CVODE but very slow on MUMPS.
Full resolution is extremely slow on CVODE (~2ms/24hrs, not included as test).
Original case name: test5ref-m2ul3a-from_lin_260519_master_scratch_cvode
Original grid: CDN_46895_lowresol_260519_nowallpump_1e21.nc
Runs for 45ms from scratch. 
Default recipe: CVODE-1

### Test 5stiff
Restarted from test5scratch at 10ms, when the simulation is more steady.
Runs for 20 output steps to get 0.01ms.
Default recipe: CVODE-1

### Test 5steady
Restarted from test5scratch at 10ms, when the simulation is more steady.
Runs for 50 output steps to get 0.5ms.
Default recipe: CVODE-1

# PETSc configuration
To enable STRUMPACK, use the following configure flags for PETSc:

```
PETSC_VERSION="3.23.3"

rm -rf "petsc-$PETSC_VERSION.tar.gz"
wget "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$PETSC_VERSION.tar.gz"
rm -rf "petsc-$PETSC_VERSION"
tar -xf "petsc-$PETSC_VERSION.tar.gz"
cd "petsc-$PETSC_VERSION"
./configure \
    COPTFLAGS="-O3" \
    CXXOPTFLAGS="-O3" \
    FOPTFLAGS="-O3"\
    --with-fortran-bindings=0 \
    --with-debugging=0 \
    --with-mpi=yes \
    --download-hypre \
    --download-make \
    --download-openblas=1 \
    --download-metis \
    --download-parmetis \
    --download-zfp \
    --download-strumpack \
    --download-scalapack \
    --download-ptscotch \
    --download-mumps \
    --download-superlu \
    --download-suitesparse \
    --download-superlu_dist \
    --download-slepc \
    --download-hpddm \
    --with-make-np=32    # This makes sure it's 32 cores for all dependencies too

make PETSC_DIR=$PWD PETSC_ARCH=arch-linux-c-opt all
```

Note that you will need to have `flex` and `bison` installed for `ptscotch` which is a STRUMPACK dependency. These are additional to the usual dependency list in the Hermes-3 documentation

