# ORCA1 Air-Sea Air-Sea Fluxes

`DOI:XXXXX.XXXXX`


## Context and Motivation

Purpose of this experiment is to compute the air-sea heat and momentum fluxes of the global ORCA1 config NEMO5 test case with different parameterizations. Air-sea fluxes and results are written in an output file with the NEMO output system (XIOS). 

#### Variations

- **W25_DET** : Mean of air-sea momentum and heat fluxes is computed with the Artifical Neural Network proposed by [Wu et al. 2025](https://doi.org/10.48550/arXiv.2503.03990). The mean is used alone as deterministic surface boundary condition. Results are compared with standard NEMO bulk formula.

#### Variations
- **W25_STO** : Mean and standard deviation of air-sea momentum and heat fluxes computed with the Artifical Neural Network proposed by [Wu et al. 2025](https://doi.org/10.48550/arXiv.2503.03990). Both are used to generate stochastic fluctuations of the surface boundary condition. The random map used in the stochastic process is generated with the NEMO STO module. Results are compared with standard NEMO bulk formula.

## Requirements

### Compilation

- NEMO version : [v4.2.1](https://forge.nemo-ocean.eu/nemo/nemo/-/releases/4.2.1) patched with [morays](https://github.com/morays-community/Patches-NEMO/tree/main/NEMO_v4.2.1) and local `CONFIG/my_src` sources.

- Code Compilation manager : none, use standard `makenemo` script


### Python

- Eophis version : [v1.1.0](https://github.com/alexis-barge/eophis/tree/v1.1.0)
- **[W25](https://github.com/jiarong-wu/mlflux)** dependencies:
  ```bash
  cd C1D_PAPA32_AirSeaFLuxes.W25/INFERENCES/W25ANN
  pip install -e .
  ```


### Run

- NEMO Production Manager : none, use submission script `job.ksh` in `RUN`


### Post-Process


- No Post-Process libraries

- Plotting : Python scripts `plots_res.py` and `plots_diff.py` in `POSTPROCESS`

