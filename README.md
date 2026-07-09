# Data for "Accelerating Polymer Self-Consistent Field Theory with Robust Preconditioned Acceleration Algorithms"

This repository contains the processed numerical data and plotting notebooks used to generate the figures and tables in the manuscript:

**Accelerating Polymer Self-Consistent Field Theory with Robust Preconditioned Acceleration Algorithms**

The SCFT calculations were performed using Polyorder.jl, an in-house code under active development in the authors' group. The source code of the development repository is not included here. Documentation for the calculation workflow and software usage is available at:

https://www.yxliu.group/Polyorder.jl/dev/

## Repository Contents

- `data/timer/`: timing data used for the computational cost analysis.
- `data/residual_plot/bcc/`: convergence data for the BCC phase.
- `data/residual_plot/sigma/`: convergence data for the Frank-Kasper sigma phase.
- `data/gyroid/`: robustness data for Gyroid-phase calculations initialized from neighboring solutions.
- `data/success_rate/`: random-initialization success-rate data.
- `notebooks/`: Jupyter notebooks used to process the CSV files and generate plots.
- `figures/`: final figure files generated from the processed data.

## Notes

The files in this repository are processed data for reproducing the plotted results and tables in the manuscript. They are not a public release of the full SCFT simulation code.

Further inquiries regarding this work can be addressed to the corresponding author.
