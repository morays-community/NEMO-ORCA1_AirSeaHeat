# ORCA1 Air-Sea Heat

`DOI:XXXXX.XXXXX`

## Context and Motivation

Purpose of this experiment is to correct the air-sea heat fluxes as a function of oceanic and atmospheric state predictors on a global ORCA1 config. More details about scientific context can be found in [Storto et al. 2024](https://doi.org/10.5194/gmd-2024-185). Corrected heat fluxes are written in an output file with the NEMO ouput system (XIOS).

#### Variations
- **S24** : Air-Sea heat fluxes correction computed with Artificial Neural Network proposed by [Storto et al. 2024](https://doi.org/10.5194/gmd-2024-185).

## Requirements

### Compilation

- NEMO version : [v4.0.7](https://forge.ipsl.fr/nemo/browser/NEMO/releases/r4.0/r4.0.7) patched with [morays](https://github.com/morays-community/Patches-NEMO/tree/main/NEMO_v4.0.7), local `CONFIG/my_src` sources and `ANNIF` module.
- **ANNIF** compilation:
	```bash
	cd ORCA1_AirSea_Heat.S24/CONFIG/ANNIF
	./Compile_gcc.ksh  # edit if necessary, should build annif.o
	```

- Code Compilation manager : none, use standard `makenemo` script
	```
	Don't forget to include annif.o module in NEMO arch file
	```

### Python

- Eophis version : [v1.0.1](https://github.com/meom-group/eophis/releases/tag/v1.0.1)
- **S24** dependencies :
	```bash
	pip install -f ORCA1_AirSea_Heat.S24/INFERENCES/requirements.txt
	```

### Run

- NEMO Production Manager : none, use submission script `job.ksh` in `RUN`


### Post-Process

- No Post-Process libraries

- Plotting : Python scripts in `POSTPROCESS`

