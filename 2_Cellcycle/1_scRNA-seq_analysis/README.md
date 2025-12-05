# scRNA-seq analyis
-> Describe what this is about
-> Describe which input files necessary



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
