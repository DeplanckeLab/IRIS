# scRNA-seq analyis

-> Describe what this is about
-> Describe which input files necessary


## General information
This folder contains: 
- the main script necessary to generate Figure 4 plots: Notebook1_Figure4_v3.ipynb
- the requirements file necessary to generate a conda evironment to run the script.

Additional necessary files are deposited in the following folders:
- `../input_files`: Contains pregenerated files to run the script.
- `../output_files`: Contains output files genrated by the script.
- `../utils`: Contains utility files/modules.

To reproduce Notebook1_Figure4_v3.ipynb:
1. First, recreate the environment (as explained below)
2. Then execute Notebook1_Figure4_v3.ipynb
3. Input files will be automatically retrieved from the input_files folder, and output files will be automatically saved in output_files folder. 


### Environment
Assumes a running jupyter-hub install and conda on the machine. Replace '<ENV_NAME> with environment name of choice.
~~~
conda create -n <ENV_NAME> python=3.12 -y
conda activate <ENV_NAME>
pip install -r requirements.txt
pip install ipykernel
python -m ipykernel install --user --name <ENV_NAME> --display-name "<ENV_NAME>"
~~~
The environment should now show up in jupyter-hub as kernel.


## Main script 
Script Notebook1_Figure4_v3.ipynb performs analyses on the single-cell transcriptomic data from (FUCCI)-expressing NIH-3T3 fibroblasts.
These include pre-processing analyses (filtering, integration, dimensionality reduction, etc.), analyses on the single-cell fucci intensities to compute the angular score and distance metric, cell cycle dynamics analyses, analyses on Dream Complex activity. 
- Anaylses are separated into chunks and introduced in each chunk.
- The script requires a series of input files. These are listed and described in chunck "Input data".
- The script generates a series of output files. These are listed and described in chunck "Output data"


