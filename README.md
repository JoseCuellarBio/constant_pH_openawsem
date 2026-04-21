# constant_pH_openawsem

Experimental constant-pH workflow for OpenAWSEM. The workflow combines OpenMM/OpenAWSEM molecular dynamics with Monte Carlo protonation-change attempts and updates the Debye-Huckel electrostatic term during the run.

The repository includes the base scripts and a prepared `1ZUG` example system.

## Repository Contents

### Root Directory

- `mm_run.py`: main simulation script.
- `forces_setup.py`: defines the OpenAWSEM force set and adds the electrostatic term that reads `charge.txt`.
- `pH_debyeHuckelTerms.py`: implements the Debye-Huckel `CustomNonbondedForce` with per-particle charges.
- `Montecarlo.py`: selects titratable residues, evaluates the local environment, and accepts or rejects protonation changes.
- `mm_run_imports.py`: alternative version of the main runner.
- `helpers/generate_charge.py`: builds an initial charge file from a FASTA sequence.
- `helpers/final_pka_plots.ipynb`: analysis notebook.

### `1ZUG/` Example

The `1ZUG/` directory contains a ready-to-run system and several generated outputs:

- input structures such as `1ZUG-openmmawsem.pdb`, `crystal_structure.pdb`, and `crystal_structure.fasta`
- AWSEM auxiliary parameter files such as `single_frags.mem`, `single_frags.npy`, `gamma.dat`, `burial_gamma.dat`, `membrane_gamma.dat`, and `ssweight`
- initial charge file `charge.txt`
- local copies of `mm_run.py`, `forces_setup.py`, `pH_debyeHuckelTerms.py`, and `Montecarlo.py`
- analysis script `mm_analyze.py`
- example outputs such as `movie.dcd`, `movie.pdb`, `output.log`, `info.dat`, `Hawsem.state`, `checkpnt.chk`, and `time.dat`

## How It Works

`mm_run.py` provides two main modes:

- `-m 0`: constant-temperature dynamics with periodic interruptions to attempt protonation changes.
- `-m 1`: temperature annealing without the constant-pH protonation cycle.

When `-m 0` is used, the workflow:

1. runs molecular dynamics for `--interruptFrequency` steps
2. builds a `Protein` object from the current geometry
3. randomly selects one titratable residue
4. evaluates a charge change with a Monte Carlo criterion that combines a pH term and a local electrostatic term
5. if the change is accepted, updates the parameters of the final force term, which is expected to be Debye-Huckel

`Hawsem.state` stores the accepted protonation-state history over the trajectory.

## Requirements

The workflow assumes a Python environment with:

- `openawsem`
- `openmm`
- `numpy`
- `pandas`
- the usual OpenAWSEM stack dependencies

Example environment activation:

```bash
conda activate openawsem
```

## Charge File

The electrostatic term and constant-pH mode expect a `charge.txt` file in the working directory. The file uses two columns:

```text
resid charge
```

where `resid` is the residue index and `charge` is its initial charge.

The helper `helpers/generate_charge.py` can create this file from a FASTA sequence:

```bash
python helpers/generate_charge.py helpers/1ZUG.fasta -o charge.txt
```

Note: the helper writes `charges.txt` by default, while the simulation looks for `charge.txt`.

## Recommended Run

With the repository in its current layout, the safest way to run the example is from `1ZUG/`, where the required input files are colocated:

- AWSEM input files
- `charge.txt`
- `forces_setup.py`
- `mm_analyze.py`

Example:

```bash
cd 1ZUG
python mm_run.py 1ZUG -s 10000 -r 100 --interruptFrequency 10 --pH 7.0 --tempStart 300 -m 0 -p CPU
```

Relevant `mm_run.py` arguments:

- `protein`: base system name
- `-s, --steps`: total number of steps
- `-r, --reportFrequency`: reporter write frequency
- `-m, --simulation_mode`: `0` for constant pH, `1` for annealing
- `--interruptFrequency`: steps between protonation attempts in mode `0`
- `--pH`: pH used by the Monte Carlo criterion
- `--tempStart` and `--tempEnd`: initial and final temperatures
- `-p, --platform`: OpenMM platform, for example `CPU`, `OpenCL`, or `CUDA`
- `-t, --thread`: number of threads when using `CPU`
- `-f, --forces`: file containing the force definitions
- `--fromCheckPoint`: restart from a checkpoint

## Typical Outputs

A run can produce:

- `output.log`: potential energy and temperature
- `movie.dcd`: trajectory
- `movie.pdb`: PDB snapshots
- `native.pdb`: reported initial structure
- `checkpnt.chk`: OpenMM checkpoint
- `time.dat`: total simulation time
- `Hawsem.state`: protonation-state history
- `info.dat`: post-run energetic analysis over the trajectory

## Current Limitations

- The repository is not packaged as a Python module; the scripts expect to be run from a directory containing specific files with fixed names.
- The root `mm_run.py` calls `mm_analyze.py` after the simulation, but that script is currently only available inside `1ZUG/`.
- Scripts are duplicated between the root directory and `1ZUG/`. For long-term maintenance, keeping a single source of truth would reduce drift.

## GitHub Publishing Notes

Before publishing or sharing the repository widely, consider whether generated heavy files should remain versioned, for example:

- `1ZUG/movie.dcd`
- `1ZUG/movie.pdb`
- `1ZUG/single_frags.npy`
- `1ZUG/checkpnt.chk`

Depending on the intended use, it may be better to keep only inputs, scripts, and workflow documentation in the repository.
