# libEMMI

`libEMMI` is a C/MPI codebase for 3D controlled-source electromagnetic (CSEM) forward modeling and full-waveform inversion using a fictitious-wave-domain FDTD formulation.

The repository contains:

- an isotropic / land-oriented solver in `src/` with headers in `include/`
- a VTI solver in `src_vti/` with headers in `include_vti/`
- small standalone utilities for model generation and adjoint-source reconstruction
- reproducible example folders used to run and visualize inversion studies

Author: Pengliang Yang
  
Email: ypl.2100@gmail.com

## Related publications

1. Pengliang Yang, "3D fictitious wave domain CSEM inversion by adjoint source estimation", *Computers & Geosciences*, 180, 105441, 2023.  
   DOI: <https://doi.org/10.1016/j.cageo.2023.105441>
2. Pengliang Yang, "libEMM: A fictitious wave domain 3D CSEM modelling library bridging sequential and parallel GPU implementation", *Computer Physics Communications*, 288, 108745, 2023.  
   DOI: <https://doi.org/10.1016/j.cpc.2023.108745>
3. Pengliang Yang and Rune Mittet, "Controlled-source electromagnetics modelling using high order finite-difference time-domain method on a nonuniform grid", *Geophysics*, 88(2), E53-E67, 2023.  
   DOI: <https://doi.org/10.1190/geo2022-0134.1>

## Repository layout

- `src/`: main isotropic solver (`bin/fdtd`)
- `include/`: headers for `src/`
- `src_vti/`: VTI solver (`bin/fdtd`, built from a different source tree)
- `include_vti/`: headers for `src_vti/`
- `bin/`: output location for the compiled solver executable
- `run_adjsrc_freq2time/`: standalone adjoint-source estimation example
- `run_create_models/`: standalone model-generation utility
- `run_fwi_land_iter1/`: land CSEM inversion example, 1 iteration
- `run_fwi_land_iter30/`: land CSEM inversion example, 30 iterations
- `run_fwi_mcsem_vti/`: marine VTI CSEM inversion example
- `doc/document.pdf`: accompanying documentation

## Requirements

- Linux
- `gcc`
- `mpicc`
- `make`
- FFTW3 development library (`libfftw3`)

Optional:

- CUDA toolkit if building with `GPU=1`
- `gfortran` for the acquisition-generation helper in some example folders
- Python 3, `gnuplot`, and Madagascar / `scons` for visualization scripts

## Build

### Isotropic solver

Build the executable from `src/`:

```bash
cd src
make
```

This creates `bin/fdtd`.

To build the CUDA path:

```bash
cd src
make GPU=1 CUDA_PATH=/path/to/cuda
```

Notes:

- the CUDA Makefile uses `-arch=sm_50` by default; update it for your GPU if needed
- `GPU=1` forces `rd3=2` in the code

### VTI solver

Build from `src_vti/`:

```bash
cd src_vti
make
```

This also writes `bin/fdtd`, so only one variant is active at a time.

## Runtime model

The main executable is an MPI program:

```bash
mpirun -n <nranks> ../bin/fdtd key1=value1 key2=value2 ...
```

The code parses command-line `key=value` arguments directly. The example `run.sh` scripts generate an `inputpar.txt` file and launch the program with:

```bash
mpirun -n <nranks> ../bin/fdtd $(cat inputpar.txt)
```

By default the acquisition code assigns one source to each MPI rank. You can override the selected source indices with `shots=...`, but `nproc` must not exceed the number of provided shots.

## Main execution modes

- `mode=0`: forward modeling
- `mode=1`: full-waveform inversion
- `mode=2`: gradient output path used by the inversion code
- `mode=3`: dot-product / adjoint test

In the isotropic path:

- `src/main.c` dispatches to `do_modeling()` or `do_fwi()`
- `src/do_modeling.c` runs time stepping, DTFT accumulation, Green's-function extraction, and writes `emf_XXXX.txt`
- `src/do_fwi.c` runs l-BFGS or Gauss-Newton optimization and writes iteration history and updated models

## Required inputs

The solver expects these groups of inputs.

### 1. Survey geometry

- `fsrc=<file>`: source table
- `frec=<file>`: receiver table
- `fsrcrec=<file>`: source-receiver connection table

The ASCII source and receiver files contain rows of:

```text
x y z azimuth dip index
```

The connection table contains source/receiver index pairs and determines which receivers belong to each shot gather.

### 2. Model files

For the isotropic solver, the executable still expects three binary resistivity volumes:

- `frho11=<file>`
- `frho22=<file>`
- `frho33=<file>`

These are read as raw single-precision floats with size `n1 * n2 * n3`. In isotropic examples the three files are typically identical or paired logically.

For the VTI solver in `src_vti/`, the corresponding inputs are:

- `frho_h=<file>`: horizontal resistivity volume
- `frho_v=<file>`: vertical resistivity volume

These are also raw single-precision float volumes of size `n1 * n2 * n3`. The VTI branch derives `rho11`, `rho22`, and `rho33` internally through homogenization.

### 3. Grid and frequency setup

Common parameters include:

- `n1`, `n2`, `n3`: model size
- `d1`, `d2`, `d3`: grid spacing
- `x1min`, `x1max`, `x2min`, `x2max`, `x3min`, `x3max`: survey bounds
- `freqs=f1,f2,...`: modeled frequencies
- `chsrc=...`: active source channels, e.g. `Ex`
- `chrec=...`: active receiver channels, e.g. `Ex,Ey`

### 4. Inversion setup

Common inversion parameters visible in `src/do_fwi.c` and the example `run.sh` files:

- `niter`: maximum nonlinear iterations
- `nls`: maximum line searches
- `npair`: l-BFGS memory length
- `method=0|1`: `0` for l-BFGS, `1` for Gauss-Newton
- `preco=0|1`: preconditioning switch
- `bound=0|1`: bound constraints
- `npar`: number of inverted parameters
- `idxpar=...`: parameter indices
- `minpar=...`, `maxpar=...`: bounds when `bound=1`
- `gamma1`, `gamma2`: regularization weights in workflows that use them

## Typical outputs

Depending on mode and example, the code writes:

- `emf_XXXX.txt`: modeled or observed frequency-domain data
- `syn_XXXX.txt`: synthetic data at the current inversion iterate
- `sig_misfit_XXXX.txt`: significant misfit per receiver
- `iterate.txt`: optimization history
- `rmse_misfit.txt`: RMSE versus iteration
- `param_final`: final inverted model in binary format

The example folders also contain plotting scripts and Madagascar `SConstruct` files for reproducing paper figures.

## Included examples

### `run_adjsrc_freq2time`

Standalone demo for reconstructing a time-domain adjoint source from a limited set of frequencies.

Run:

```bash
cd run_adjsrc_freq2time
make
./main
python3 plot_basis_function.py
```

The executable writes `basis_function`, which the plotting scripts visualize.

### `run_create_models`

Standalone utility for generating binary resistivity models.

Run:

```bash
cd run_create_models
make
./main
```

Edit `main.c` in that directory to change the synthetic model definition. The folder includes Matlab scripts and example outputs for visualization.

### `run_fwi_land_iter30`

Reference land CSEM inversion example using the isotropic solver.

Suggested workflow:

```bash
cd src
make
cd ../run_fwi_land_iter30
bash run.sh
```

What `run.sh` currently does:

- writes `inputpar.txt`
- sets `mode=1`
- launches `mpirun -n 16 ../bin/fdtd $(cat inputpar.txt)`

Before running, confirm that:

- `sources.txt`, `receivers.txt`, and `src_rec_table.txt` match your survey
- `rho11_init`, `rho22_init`, and `rho33_init` exist and match `n1*n2*n3`
- your MPI rank count is consistent with the number of active shots

Visualization helpers in this folder:

- `plot_emdata.py`
- `plot_histogram.py`
- `plot_scatter_sigmisfit.py`
- `plot_survey_layout.gnu`
- `SConstruct`

### `run_fwi_land_iter1`

Same workflow as `run_fwi_land_iter30`, but configured for a single inversion iteration and useful as a smoke test.

### `run_fwi_mcsem_vti`

Reference marine VTI inversion example. Build `bin/fdtd` from `src_vti/` first, then run:

```bash
cd run_fwi_mcsem_vti
bash run.sh
```

This example also uses:

- nonuniform vertical gridding (`nugrid=1`)
- bathymetry input (`ibathy`)
- separate horizontal and vertical resistivity parameters

## Practical notes

- `bin/fdtd` is overwritten when you rebuild from `src/` or `src_vti/`
- the example folders contain generated artifacts and paper-reproduction outputs; not every file is hand-maintained source
- model binaries are raw floats with no embedded metadata, so shape consistency is critical
- some visualization steps assume Madagascar; if you do not use it, you can still post-process the binary outputs with your own tools

## Citation

If you use this code in academic work, cite the papers listed above.
