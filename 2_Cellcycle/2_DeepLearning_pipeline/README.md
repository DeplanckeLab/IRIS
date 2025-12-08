# Image to transcriptome modelling

### Utilization
This folder containes all scripts to train the cell shape CNN, nucleus CNN, and fucci angle MLP. Scripts are run in the following order:

1. Training of models (see Traning the model):
    - CNN models: `image2transcriptome.py`
        - On Hoechst '405' for the nucleus shape model.
        - On the segmentation mask 'seg' for the cell shape model.
    - MLP models: `image2transcriptome_feature.py`
2. Results are then compiled in the following notebook: `compile_scores.ipynb`

The compiled data is analyzed/plotted in the 'DeepLearning_analysis' section. The specific calls of all scripts are listed below in 'Example calls'

### Installation

We recommend using Python 3.12. Replace `<ENV_NAME>` with name of choice.
~~~
conda create -n <ENV_NAME> python=3.12 -y
conda activate <ENV_NAME>
pip install -r requirements.txt
~~~
---

### Training the model
The `image2transcriptome.py` script supports 2 split types:
1. Random splitting: Use argument `--random_split`
2. Split by batch: Pass arguments `--train_experiments_prefix 'BATCH_PREFIX_TRAIN'` and `--val_experiments_prefix 'BATCH_PREFIX_VAL'`, _e.g._ `BATCH_PREFIX_TRAIN='JP'` indicates that JP's experimental batches will be used for training. The remaining batches that don't match `BATCH_PREFIX_TRAIN` or `BATCH_PREFIX_VAL` will be used for evaluation.

```
python image2transcriptome.py
       --data_dir <DATA_DIR_PATH>
       --out_dir <OUT_DIR_PATH>
       --dataset <DATASET_NAME>
       --channel <CHANNEL_NAME>
       --apply_segmentation_mask
       --model_name 'resnet'
       --save_model
       --species 'mouse'
       --device 0
       --seeds 0 1 2 3 4
```
This script trains multiple models in a leave-one-batch-out way and stores predictions for each test batch. 

Description of the main arguments:
* `data_dir`: Path of the directory containing the H5 data.
* `out_dir`: Path of the directory where the script outputs will be saved.
* `dataset`: Name of the compiled H5 file.
* `cell_types`: List of cell types to select. By default, all cell types will be included.
* `species`: List of species to select. By default, all species will be included.
* `channel`: Name of the channel in the H5 file that will be used to make predictions, e.g.: `405`, `seg`.
* `apply_segmentation_mask`: Whether to apply the segmentation mask to the images.
* `model_name`: Name of the model to train (`resnet` by default). Currently supported options: `resnet`.
* `save_model`: Whether to save the trained models.
* `discard_spikein`: Whether to discard spikein cells based on the average expression of a few pre-specified genes.
* `device`: GPU number.
* `seeds`: List of seeds used to initialize the parameters of the model. The script is run for each provided seed, storing results for each independent run.
* `n_folds`: Number of folds. The train data is split into `n_folds`. An independent model will be trained on `n_folds - 1` folds, using the remaining fold for early stopping. This step yields `n_folds` trained models. We run inference on the test set using each of these models and then average the predicted expression profiles.
* Other parameters include `bs` (batch size), `epochs` (maximum number of training epochs), `lr` (learning rate), and `wd` (weight decay).

**About the other scripts**
The script `loo_image2transcriptome_feature.py` similarly trains baseline models based on predefined features that can be specified via the `feature_name` argument (supported arguments: `size`, `roundness`, `eccentricity`, `axis_minor`, and `axis_major`).

---

### Example calls
Nucleus CNN model:
~~~
python image2transcriptome.py --data_dir ../input_files --out_dir ./results --dataset fucci_final --channel '405' --save_model --species 'mouse' --device 0 --seeds 0 1 2 3 4 --random_split
~~~
Cell shape CNN model:
~~~
python image2transcriptome.py --data_dir ../input_files --out_dir ./results --dataset fucci_final --channel 'seg' --train_augmentations 'random_rotation' --save_model --species 'mouse' --device 0 --seeds 0 1 2 3 4 --random_split
~~~
Angular MLP model
~~~
python image2transcriptome_feature.py --data_dir ../input_files --out_dir ./results --dataset fucci_final --channels 488 561 seg --fit_method mlp --species 'mouse' --device 0 --seeds 0 1 2 3 4 --random_split
~~~
---

### Visualizing the results

By default, for each training seed, the training script:
1. Computes batch-specific Pearson correlation scores for each gene.
2. Stores a scatter plot of the predicted expression profiles for the top 10 predicted genes.
3. Saves an AnnData object with the ground-truth and predictions. This object can later be loaded to visualize a UMAP representation of the predicted expression profiles.
4. Stores the results aggregated across multiple runs.
5. Saves the trained models for each training fold, if `save_model` is provided.

The notebook `comparison_loo.ipynb` and `comparison.ipynb` visualize the Pearson correlation scores for multiple channels and seeds.
