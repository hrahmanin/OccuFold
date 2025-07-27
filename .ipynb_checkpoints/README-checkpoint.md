# OccupancyInputCTCF



### Description
This GitHub repository provides tools for training machine learning models to predict 3D chromatin architecture from single-molecule footprinting data (e.g., methylation patterns). It integrates sequence features and occupancy profiles to infer genome folding and validates predictions by comparing Hi-C data with simulated chromatin loop extrusion, incorporating locus-specific occupancy rates and dynamic CTCF barriers.

![Workflow Figure](figures/workflowfigure.png)

### Structure of the repository
The structure of this repository follows as below:
- processing/ Scripts and pipelines for NGS data processing (e.g., handling SMF methylation footprint data and ChIP-seq data).
- models/ – Code for deep learning models (CNN architectures, training scripts, evaluation functions) used to predict CTCF occupancy or 3D contacts
- analysis/ – Notebooks or scripts for analyzing results (e.g. comparing predicted vs. actual Hi-C, generating figures).
- utils/ – Utility functions and tools (shared helper code for data I/O, metric calculations, etc.).
- outputs/ – Folder to store output files, such as processed data or model predictions (keeping them separate from code)

  

<!--### Requirements
- *Polychrom*: A toolkit for polymer simulations. (https://github.com/open2c/polychrom)
- *OpenMM*: A library for molecular simulations. (https://github.com/openmm/openmm)
- *Open2C* analysis packages (see https://github.com/open2c)-->

  
## Installation
First, 

```
git clone https://github.com/Fudenberg-Research-Group/OccupancyInputCTCF.git
```

### Workflow
#### Running simulations 
1. One-Dimensional Lattice Simulation: with running `workflow.py`


#### Processing simulation data
After running the simulations, the simulated trajectories can be processed to generate *in silico* ChIP-seq profiles, 1d contact maps, and 3d contact maps (optional). Scripts for data processing available in `processing`. Instructions are provided with the relevant python code.

#### Analysis
Once the data is processed, observable features can be quantified




