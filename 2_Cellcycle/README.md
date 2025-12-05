# Cell cycle analysis
## Introduction
This folder contains the scripts that are necessary to recreate the results shown in Figure 4 and 5. The scripts are separated into 3 discrete scripts:

1. FUCCI scRNA-seq analysis (Figure 4): ~~~1_scRNA-seq_analysis~~~
2. Deep-learning pipeline for cell shape, nuclear, and angular model (Figure 5): ~~~2_DeepLearnging_pipeline~~~
3. Analysis of deep-learning results (Figure5): ~~~3_DeepLearnging_analysis~~~

## Execution extructions
### Execution order
The scripts are partially dependent on eachother. We provide all temporary results necessary that make the scripts executable individually. The only file necessary for this is X.h5. For regeneration of all files we suggest the following order of execution:

1. FUCCI_scRNA-seq can be fully executed except for the last cell. Produces the following files:
- Requires:
  - a
  - b
- Generates:
  - a
  - b
2. Execution of the machine learning pipeline.
- Requires:
  - a
  - b
- Generates:
  - a
  - b
3. Analysis of deep-learning 
- Requires:
  - a
  - b
- Generates:
  - a
  - b
### Dependencies
For each notebook a requirements.txt is provided to create an environment compatible with the notebooks.
