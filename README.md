**BDQR Transmission Threshold Model**

This repository contains code and data supporting the manuscript: Estimating the relative contribution of transmission to bedaquiline resistance burden in tuberculosis. A transmission threshold framework. It includes R scripts for data cleaning, parameter preparation, and analysis, plus a C++ stochastic within-host birth–death simulation used to explore spontaneous BDQR dynamics.

**Repository Structure**
- **c_simulation/**: C++ stochastic within-host birth–death simulation for spontaneous BDQR. Build with CMake (see below).
- **bdqr_prev_data/**: BDQR prevalence data collected from a systematic review.
- **0_global_setup.R**: BDQR transmission threshold model and shared helper functions used throughout the project.
- **1_clean_data_sheets.R**: Scripts to clean WHO data and prevalence data when authors supplied data sheets.
- **1_data_setup.R**: Prepare input data for the threshold model.
- **1_parameter_setup.R**: Prepare input parameters (e.g., probability of acquisition and spontaneous resistance).
- **2_analysis.R**: Code that reproduces the analyses presented in the manuscript and appendix.
- **figures/**: Output figures produced by the analysis scripts.
- **bdqr_prev_data/** and **who_data/**: Supporting datasets used by the analysis.

**Quick start**
1. Requirements
   - R (compatible with the scripts in this repo; tested with R >= 4.0)
   - `cmdstanr`
   - CMake and a C++ compiler to build `c_simulation` (e.g., `gcc`/`clang`)

If this is the first time running `cmdstanr` you need to install the R package (see https://mc-stan.org/cmdstanr/articles/cmdstanr.html for more information):

```r
# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
check_cmdstan_toolchain()
# If your toolchain is configured correctly then CmdStan can be installed by calling the install_cmdstan() function:
install_cmdstan(cores = 2)
```

2. Build the simulation (from the repository root)

```bash
cd c_simulation
mkdir -p build && cd build
cmake ..
make
```

Running the simulation from the terminal produces `out.csv`, which contains the relevant data for calculating the probability of spontaneous resistance. 

3. Run the R analysis pipeline
   - Run `2_analysis.R` to reproduce figures and tables for the manuscript. Running `2_analysis.R` automatically calls:

```r
source('0_global_setup.R')
source('1_clean_data_sheets.R')
source('1_data_setup.R')
source('1_parameter_setup.R')
```

- Data files used by the scripts live in `bdqr_prev_data` and `who_data`.
