import h5py
import json
import pandas as pd
import numpy as np
import scanpy as sc

# Function to load gene expression data
def get_exp_data(gene_exp: bytes) -> dict:
    decode_data = gene_exp.decode("utf-8")
    decode_data = decode_data.replace("\'","\"")
    exp_dict = json.loads(decode_data)
    return exp_dict

def select_cell(cell_ds, channels, cell_types=None, species=None, exp_data=True, exclude_batches=None):
    if exp_data and "exp_data" not in cell_ds or 'exp_data/exp_data' not in cell_ds or not cell_ds.attrs['num_cells'] == 1:
        return False
    
    cell_type_check = cell_types is None
    if not cell_type_check:
        for c in cell_types:
            cell_type_check = cell_ds.attrs['Cell_type'] == c
            if cell_type_check:
                break

    species_check = species is None
    if not species_check:
        if 'species_annotation' not in cell_ds['exp_data'].attrs:  # Keep all cells of unknown species (probably mouse)
            species_check = True
        else:
            for s in species:
                species_check = cell_ds['exp_data'].attrs['species_annotation'] == s
                if species_check:
                    break

    batch_check = exclude_batches is None
    if not batch_check:
        for b in exclude_batches:
            batch_check = cell_ds.attrs['exp_id'] != b
            if not batch_check:
                break

    return (not exp_data or "exp_data" in cell_ds or 'exp_data/exp_data' in cell_ds) \
        and all([c in cell_ds for c in channels]) \
        and cell_ds.attrs['num_cells'] == 1 \
        and cell_type_check \
        and species_check \
        and batch_check

def load_data(ds_dir, channels=['405', 'bf', 'seg', '488', '561'], cell_types=None, species=None, exp_data=True, exclude_batches=None):
    list_dict = []
    channels_dict = {k: [] for k in channels}
    with h5py.File(ds_dir, 'r') as dataset:
        print(len(dataset.keys()))
        for cell_id in list(dataset.keys()):
            cell_ds = dataset[str(cell_id)]
            
            # Some cells don't have have a 405 channel. Filter single cells.
            if select_cell(cell_ds, channels, cell_types=cell_types, species=species, exp_data=exp_data, exclude_batches=exclude_batches):
                # Get the plane ids from the bf channel
                if 'bf' in cell_ds:
                    planes = list(cell_ds['bf'].keys())
                elif '405' in cell_ds:
                    planes = list(cell_ds['405'].keys())
                else:
                    raise ValueError('Cannot extract planes')
                # I extracted +/- 1 plane from focus. So the 2nd list index is the in focus plane
                # Careful: The ids have absolut index of their original planes.
                foc_plane_id = planes[1]

                if exp_data:
                    exp_data = cell_ds['exp_data/exp_data'][()]
                    exp_series = get_exp_data(exp_data)
                    exp_series['cell_id'] = cell_id
                    exp_series['cell_type'] = cell_ds.attrs['Cell_type']
                    exp_series['exp_id'] = cell_ds.attrs['exp_id']
                    if 'species_annotation' in cell_ds['exp_data'].attrs:
                        exp_series['Cell_species'] = cell_ds['exp_data'].attrs['species_annotation'] # cell_ds.attrs['Cell_species']
                    else:
                        exp_series['Cell_species'] = 'Unknown'
                    if 'well' in cell_ds.attrs:
                        split_id = cell_ds.attrs["exp_full_id"].split('_')
                        exp_series['Full_Cell_ID'] = f'{"_".join(split_id[:2])}_WC{cell_ds.attrs["well"]}_{split_id[2]}'
                    list_dict.append(exp_series)
                # print(cell_id, cell_ds.attrs['Cell_type'], cell_ds.attrs['exp_id'])

                # Get image, fucci, and segmentation mask
                for channel in channels:
                    img = cell_ds[channel][foc_plane_id][()]
                    """
                    if channel != 'seg':
                        min_val = np.min(img)
                        max_val = np.max(img)
                        img = (img - min_val) / (max_val - min_val) * 255
                    """
                    channels_dict[channel].append(img)

    # Store images in numpy array
    for k, v in channels_dict.items():
        channels_dict[k] = np.array(v)

    # Store expression data
    adata = None
    if exp_data:
        exp_df = pd.DataFrame(list_dict)
        if 'cell_id' in exp_df:
            exp_df = exp_df.set_index('cell_id')
        obs = exp_df[['cell_type', 'exp_id', 'Cell_species', 'Full_Cell_ID']]
        X = exp_df.drop(columns=['cell_type', 'exp_id', 'Cell_species', 'Full_Cell_ID'])
        print('Full_Cell_ID' in X.columns)
        adata = sc.AnnData(X, obs=obs)
        nan_mask = np.any(np.isnan(adata.X), axis=0)
        adata = adata[:, ~nan_mask]

    return adata, channels_dict