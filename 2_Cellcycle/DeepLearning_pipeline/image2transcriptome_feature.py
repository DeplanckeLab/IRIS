import torch
from torch.utils.data import DataLoader
import torch.nn.functional as F
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import json
import argparse
from pathlib import Path
import logging
import cv2
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
from skimage.measure import regionprops
from scipy.stats import kurtosis
from scipy.ndimage import generic_filter
from src.data import load_data
from src.eval import calculate_results_per_batch, adata_results
from src.plotting import scatter_top_genes, plot_example_images
sns.set_style('whitegrid')

parser = argparse.ArgumentParser()
parser.add_argument('--seeds', default=[0, 1, 2], nargs='+')
parser.add_argument('--suffix', default=None)
parser.add_argument('--data_dir')
parser.add_argument('--out_dir')
parser.add_argument('--feature_name', default='roundness')
parser.add_argument('--feature_names', default=None, nargs='+')
parser.add_argument('--fit_method', default='linear')
parser.add_argument('--dataset', default='fucci_3t3_221124')
parser.add_argument('--cell_types', default=None, nargs='+')
parser.add_argument('--species', default=None, nargs='+')
parser.add_argument('--exclude_batches', default=None, nargs='+')
parser.add_argument('--channel', default=None)
parser.add_argument('--channels', default=None, nargs='+')
parser.add_argument('--minmax_normalize', action='store_true')
parser.add_argument('--apply_segmentation_mask', action='store_true')
parser.add_argument('--min_genes', default=500, type=int)
parser.add_argument('--min_cells', default=5, type=int)
parser.add_argument('--n_top_genes', default=2000, type=int)
parser.add_argument('--gene_list_file', default=None)
parser.add_argument('--device', default=None)
parser.add_argument('--random_split', action='store_true')
parser.add_argument('--train_experiments_prefix', default='NG')
parser.add_argument('--val_experiments_prefix', default=None)
parser.add_argument('--train_augmentations', default='albumentations')
parser.add_argument('--test_augmentations', default=None)
parser.add_argument('--shuffle_size', default=None, type=int)
parser.add_argument('--gaussblur_size', default=None, type=int)
parser.add_argument('--bs', default=256, type=int)
parser.add_argument('--num_workers', default=8, type=int)
parser.add_argument('--emb_dim', default=32, type=int)
parser.add_argument('--backbone_size', default='18', type=str)
parser.add_argument('--epochs', default=2000, type=int)
parser.add_argument('--lr', default=1e-2, type=float)
parser.add_argument('--wd', default=1e-4, type=float)
parser.add_argument('--discard_spikein', action='store_true')
parser.set_defaults(random_split=False)
parser.set_defaults(minmax_normalize=False)
parser.set_defaults(apply_segmentation_mask=False)
parser.set_defaults(discard_spikein=False)
args = parser.parse_args()


def calculate_avg_mask_intensities(seg, im):
    seg_uint8 = cv2.threshold(seg, 2, 65535, cv2.THRESH_BINARY)[1] #seg is a binary matrix of 0 or 255 pixel values
    count=np.sum(seg_uint8)/255 #we count the total nb of pixels that belong to the cell
    if count!=0: #if segmentation exists
        masked_image = cv2.bitwise_and(im, im, mask=seg_uint8)
        sum_val = np.sum(masked_image)  #here we count the intensity
        mean = sum_val/count
    else:
        mean=0
    return mean


class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_dim, output_dim)
    
    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out

def fit_predict_mlp(X_train, X_val, X_test, y_train, y_val, y_test,
                    device=0,
                    hidden_dim=32,
                    num_epochs=500,
                    learning_rate=0.01,
                    patience=2):
    device = torch.device(f"cuda:{device}" if torch.cuda.is_available() else "cpu")
    input_dim = X_train.shape[1]
    output_dim = y_train.shape[1]
    
    # Convert to PyTorch tensors
    X_train = torch.tensor(X_train, dtype=torch.float32).to(device)
    X_val = torch.tensor(X_val, dtype=torch.float32).to(device)
    X_test = torch.tensor(X_test, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32).to(device)
    y_val = torch.tensor(y_val, dtype=torch.float32).to(device)
    y_test = torch.tensor(y_test, dtype=torch.float32)
    
    # Initialize model, loss function, and optimizer
    model = MLP(input_dim, hidden_dim, output_dim)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    model.to(device)
    
    # Training loop
    for epoch in range(num_epochs):
        model.train()
        optimizer.zero_grad()
        outputs = model(X_train)
        loss = criterion(outputs, y_train)
        loss.backward()
        optimizer.step()
        
        # Validation
        model.eval()
        early_stop_counter = 0
        with torch.no_grad():
            val_outputs = model(X_val)
            val_loss = criterion(val_outputs, y_val)
            if epoch == 0:
                best_val_loss = val_loss
            else:
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_model_state = model.state_dict()
                else:
                    early_stop_counter += 1
                    if early_stop_counter >= patience:
                        print(f"Early stopping at epoch {epoch+1}")
                        model.load_state_dict(best_model_state)
                        break
        
        if (epoch+1) % 10 == 0:
            print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, Val Loss: {val_loss.item():.4f}')
    
    # Evaluate on test set
    model.eval()
    with torch.no_grad():
        test_outputs = model(X_test.to(device)).cpu().numpy()

    return {'test_outputs': test_outputs}


def train_and_store_results(args, adata_subset, imgs, seg_masks, split_idxs, seed, out_dir, obsm_key = None):    
    # Calculate features
    features = []
    feature_names = args.feature_names if args.feature_names is not None else [args.feature_name]
    feature_names_str = '-'.join(feature_names)
    print(f'Calculating features {feature_names_str} ...')
    for i, (m, img) in tqdm(enumerate(zip(seg_masks, imgs))):
        props = regionprops((m > 0).astype(int))
        masked_pixels = img[m > 0]
        # Assuming there's only one region (the cell)
        if len(masked_pixels) == 0 or len(props) == 0:
            # Case where segmentation masks are all zeros
            feature_value = np.nan
        else:
            feature_values = []
            for feature_name in feature_names:
                if feature_name == 'roundness':
                    area = props[0].area 
                    perimeter = props[0].perimeter
                    feature_value = (4 * np.pi * area) / (perimeter ** 2)
                elif feature_name == 'eccentricity':
                    feature_value = props[0].eccentricity
                elif feature_name == 'axis_minor':
                    feature_value = props[0].axis_minor_length
                elif feature_name == 'axis_major':
                    feature_value = props[0].axis_major_length
                elif feature_name == 'solidity':
                    feature_value = props[0].solidity
                elif feature_name == 'size':
                    # feature_value = props[0].area
                    feature_value = np.mean((m > 0).astype(int))
                elif feature_name == 'fucci_angle':
                    feature_value = adata_subset[i].obs['fucci_angle'].values.ravel()[0]
                elif feature_name == 'kurtosis':
                    feature_value = kurtosis(masked_pixels, fisher=True, bias=False)
                elif feature_name == 'maes_max':
                    feature_value = np.max(masked_pixels) / 65535
                elif feature_name == 'maes_mean':
                    feature_value = np.mean(masked_pixels) / 65535
                elif feature_name == 'maes_var':
                    feature_value = np.var(masked_pixels / 65535)
                else:  
                    raise ValueError(f'Unknown feature name: {feature_name}')
                feature_values.append(feature_value)
        features.append(feature_values)
    features = np.array(features)
    
    if len(feature_names) == 1 and len(features.shape) == 1:
        features = features[..., None]
    
    # Create output dir and store arguments
    results_dir = f'{out_dir}/{args.fit_method}_{feature_names_str}_model/seed_{seed}'
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    with open(f'{results_dir}/args.txt', 'w') as f:
        json.dump(args.__dict__, f, indent=2)
    logging.info(f'Features shape: {features.shape}')
    assert features.shape[0] == adata_subset.shape[0], 'Features shape does not match adata shape'
    
    # Logging
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    logging.FileHandler(f'{out_dir}/outputs.log'),
    log.addHandler(logging.FileHandler(f'{results_dir}/outputs.log'))      # set the new handler
    log.addHandler(logging.StreamHandler())      # set the new handler

    logging.info('... instantiating data and model')
    train_idxs = split_idxs['train']
    val_idxs = split_idxs['val']
    test_idxs = split_idxs['test']
    X_train = features[train_idxs]
    X_val = features[val_idxs]
    X_test = features[test_idxs]
    train_gts = adata_subset[train_idxs].X
    val_gts = adata_subset[val_idxs].X
    test_gts = adata_subset[test_idxs].X

    # Discard NaN values
    nan_mask = np.isnan(X_train).any(axis=1)
    if np.any(nan_mask):
        logging.info(f'... discarding {np.sum(nan_mask)} NaN samples from training set')
        X_train = X_train[~nan_mask]
        train_gts = train_gts[~nan_mask]
    nan_mask = np.isnan(X_val).any(axis=1)
    if np.any(nan_mask):
        logging.info(f'... discarding {np.sum(nan_mask)} NaN samples from val set')
        X_val = X_val[~nan_mask]
        val_gts = val_gts[~nan_mask]
    nan_mask = np.isnan(X_test).any(axis=1)
    if np.any(nan_mask):
        logging.info(f'... discarding {np.sum(nan_mask)} NaN samples from test set')
        X_test = X_test[~nan_mask]
        test_gts = test_gts[~nan_mask]

    # Training model
    if args.fit_method == 'linear':
        # Fit linear regression
        logging.info('... fitting linear regression')
        reg = LinearRegression()
        reg.fit(X_train, train_gts)
        test_preds = reg.predict(X_test)
    elif args.fit_method == 'mlp':
        # Fit MLP
        logging.info('... fitting MLP')
        results = fit_predict_mlp(X_train, X_val, X_test, train_gts, val_gts, test_gts,
                                  device=args.device,
                                  hidden_dim=args.emb_dim,
                                  num_epochs=args.epochs,
                                  learning_rate=args.lr,
                                  patience=2)
        test_preds = results['test_outputs']

    # Run inference on test datasets (including different inference augmentations)
    df_check = calculate_results_per_batch(adata_subset[test_idxs], test_preds, test_gts)
    df_check.to_csv(f'{results_dir}/test_results.csv')

    # Plot results
    plt.figure(figsize=(16, 6))
    scatter_top_genes(test_preds, test_gts, df_check, adata_subset[test_idxs], k=10)
    plt.tight_layout()
    plt.savefig(f'{results_dir}/test_scatter_top_genes.pdf', bbox_inches='tight')

    # Store results for the main test set
    adata_pred = adata_results(adata_subset[test_idxs], test_preds, test_gts)
    adata_pred.write_h5ad(f'{results_dir}/test_preds.h5ad')

    """
    # Store train results
    logging.info('... storing train results')
    train_preds, train_gts, train_output_features = inference(model, dataloaders['train'], device)
    adata_pred = adata_results(adata_subset[train_idxs], train_preds, train_gts,
                               output_features=train_output_features)
    adata_pred.write_h5ad(f'{results_dir}/train_preds.h5ad')
    """

if __name__ == '__main__':
    # Gather channels
    channels = []
    if args.channels is not None:
        channels = args.channels
    if args.channel is not None:
        channels.append(args.channel)
    channels = list(sorted(list(set(channels))))  # Remove duplicates
    print(f'Using channels: {channels}')

    # Create output directory and logs
    channels_string = '-'.join(channels)
    ct_string = ''
    segm_string = ''
    split_string = ''
    suffix_string = ''
    if args.cell_types is not None:
        ct_string = '-' + '-'.join(args.cell_types)
    if not args.apply_segmentation_mask:
        segm_string = '-bg'
    if args.random_split:
        split_string = '-random'
    if args.suffix is not None:
        suffix_string = '_'+args.suffix
    out_dir = f'{args.out_dir}/{args.dataset}{ct_string}_{channels_string}{segm_string}{split_string}{suffix_string}'
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(f'{out_dir}/outputs.log'),
            logging.StreamHandler()
        ]
    )
    logging.info(f'Image2transcriptome model. Storing results in {out_dir}')
    logging.info('... loading data')
    channels_ = channels.copy()
    if args.apply_segmentation_mask and args.channel != 'seg':
        channels_.append('seg')
    
    # Load data
    adata, channels_dict = load_data(ds_dir=f'{args.data_dir}/{args.dataset}.h5',
                                     channels=channels_,
                                     cell_types=args.cell_types,
                                     species=args.species,
                                     exclude_batches=args.exclude_batches,
                                     exp_data=True)
    logging.info(f'Adata shape: {adata.shape}')

    # Detect low quality cells
    logging.info('... detecting low quality cells')
    filter_mask, _ = sc.pp.filter_cells(adata, min_genes=args.min_genes, inplace=False)
    adata = adata[filter_mask]
    for k, v in channels_dict.items():
        # Discard matching images
        channels_dict[k] = v[filter_mask]

    # Normalize
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Detect spike-in cells
    if args.discard_spikein:
        logging.info('... detecting spike-in cells')
        spike_in_genes = ['Gm13194', 'Gm5781', 'Gm16073', 'Gm45499', 'Gm16379', 'Gm7618', 'Gm18014', 'Gm12152']  # 'Gm4321'
        spike_in_genes = [g for g in spike_in_genes if g in adata.var.index]
        spike_in_mean = adata[:, spike_in_genes].X.mean(axis=-1)
        spike_in_mask = (spike_in_mean <= 2)
        logging.info(f'... discarding {(~spike_in_mask).sum()} spike-in cells')
        adata = adata[spike_in_mask]
        for k, v in channels_dict.items():
            # Discard matching images
            channels_dict[k] = v[spike_in_mask]
        
    # Select HVGs
    if args.gene_list_file is not None:
        logging.info(f'... loading gene list from {args.gene_list_file}')
        with open(args.gene_list_file, 'r') as f:
            gene_list = [line.strip().upper() for line in f.readlines()]
        n_genes = len(gene_list)
        genes = np.array([g.upper() for g in adata.var.index])
        gene_idxs = [np.argwhere(genes == g).ravel()[0] for g in gene_list if g in genes]
        logging.info(f'... selecting {len(gene_idxs)}/{n_genes} genes from the input list')
        adata_subset = adata[:, gene_idxs]
    else:
        logging.info(f'... selecting top {args.n_top_genes} HVGs')
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, batch_key="exp_id")
        g = adata.var.highly_variable
        adata_subset = adata[:, g.values]
        logging.info(f'Prc1 in data: {"Prc1" in adata_subset.var.index}')
        
    # Apply segmentation masks
    if args.apply_segmentation_mask:
        logging.info('... applying segmentation masks')
        seg_masks = channels_dict['seg'] / 255
        for k, v in channels_dict.items():
            if k != 'seg':
                channels_dict[k] = v * seg_masks
            # channels_dict[args.channel] = channels_dict[args.channel] * seg_masks
            logging.info(f'Min: {channels_dict[k].min()}. Max: {channels_dict[k].max()}')
        if not 'seg' in channels:
            del channels_dict['seg']

    # Splitting data by experiment ID
    logging.info('... splitting data')
    df = adata_subset.obs
    
    for seed in args.seeds:
        np.random.seed(int(seed))
        logging.info(f'Seed: {seed}')

        if args.random_split:
            # Split data at random, not by experiment ID
            n_total = df.shape[0]
            train_idxs = np.random.choice(n_total, int(0.6*n_total), replace=False)
            val_test_idxs = np.setdiff1d(np.arange(n_total), train_idxs)
            val_idxs = np.random.choice(val_test_idxs, int(0.5*len(val_test_idxs)), replace=False)
            test_idxs = np.setdiff1d(val_test_idxs, val_idxs)
        else:
            experiments = df['exp_id'].unique()
            train_experiments = [exp for exp in experiments if exp.startswith(args.train_experiments_prefix)]
            val_experiments = [] if args.val_experiments_prefix is None else [exp for exp in experiments if exp.startswith(args.val_experiments_prefix)]
            test_experiments = [exp for exp in experiments if exp not in train_experiments and exp not in val_experiments]  # [exp for exp in experiments if exp.startswith('JP241')]
            assert set(train_experiments + val_experiments + test_experiments) == set(experiments)

            train_df = df[df['exp_id'].isin(train_experiments)]
            val_df = df[df['exp_id'].isin(val_experiments)]
            test_df = df[df['exp_id'].isin(test_experiments)]
            train_idxs = np.argwhere(adata_subset.obs['exp_id'].isin(train_experiments).values).ravel()
            val_idxs = np.argwhere(adata_subset.obs['exp_id'].isin(val_experiments).values).ravel()
            test_idxs = np.argwhere(adata_subset.obs['exp_id'].isin(test_experiments).values).ravel()
            n_total = df.shape[0]
            logging.info(f'Train size: {train_df.shape[0]}/{n_total} ({100*train_df.shape[0]/n_total:.2f}%)')
            logging.info(f'Val size: {val_df.shape[0]}/{n_total} ({100*val_df.shape[0]/n_total:.2f}%)')
            logging.info(f'Test size: {test_df.shape[0]}/{n_total} ({100*test_df.shape[0]/n_total:.2f}%)')

        # Segmentation masks
        seg_masks = channels_dict['seg'] / 255
        split_idxs = {
            'train': train_idxs,
            'val': val_idxs,
            'test': test_idxs
        }

        # Normalize images
        logging.info('... normalizing images')
        imgs = []
        ch = args.channel if args.channel is not None else channels[0]
        for img in tqdm(channels_dict[ch]):
            min_val = np.min(img)
            max_val = np.max(img)
            img = (img - min_val) / (max_val - min_val+1e-8)
            imgs.append(img)
        imgs = np.array(imgs)

        # Calculate angles from FUCCI intensities
        if '488' in channels_dict and '561' in channels_dict:
            fucci1_data_norm = []
            fucci2_data_norm = []
            for i in range(adata.shape[0]):
                fucci1 = channels_dict['488'][i]
                fucci2 = channels_dict['561'][i]
                seg = channels_dict['seg'][i]
                fucci1_data_norm.append(calculate_avg_mask_intensities(seg, fucci1))
                fucci2_data_norm.append(calculate_avg_mask_intensities(seg, fucci2))
            adata_subset.obs['fucci1'] = fucci1_data_norm
            adata_subset.obs['fucci2'] = fucci2_data_norm

            logging.info('... calculating angles from FUCCI intensities')
            pca = PCA(n_components=2)
            pcs = pca.fit_transform(np.log1p(adata_subset.obs[['fucci1', 'fucci2']]))
            angles = np.arctan2(pcs[:, 1], pcs[:, 0])
            adata_subset.obs['fucci_angle'] = angles

        # Train and store results for each seed
        logging.info('... training and storing results')
        train_and_store_results(args, adata_subset, imgs, seg_masks, split_idxs, seed, out_dir)
