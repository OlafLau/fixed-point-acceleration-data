# Calculation drivers

These Julia files are the calculation-level drivers associated with the reported data. They document the model and solver calls; machine scheduling and recovery scripts are not included.

## Environment

Use a Julia environment containing Polyorder, Polymer, Scattering, FFTW, CSV, DataFrames, JLD2, JSON, and HDF5 as required by each driver. A typical invocation is:

```bash
julia --project=/path/to/polyorder-project script.jl [arguments]
```

For the Gyroid and HEX scripts, `POLYORDER_PROJECT=/path/to/polyorder-project` may be set instead. All script includes and default output locations are resolved relative to the script directory.

## Groups

- `benchmark/principal_worker.jl` reads a run directory containing `plan.json` and an `inputs/` directory, then writes traces and summaries under that run directory.
- `benchmark/timing_m_bcc.jl` performs the BCC history-length timing calculation. `SCFT_INPUT_DIR` selects the directory containing `bcc_seed.h5`; `SCFT_OUTPUT_DIR` selects the output directory.
- `discretization/` contains the shared worker and the contour/spatial scan entry point.
- `gyroid/` contains the neighboring-solution parameter scan.
- `hex/` contains the three warm-up batch drivers.
