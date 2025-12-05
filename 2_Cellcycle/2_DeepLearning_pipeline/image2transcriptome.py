import torch
from torch.utils.data import DataLoader
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import json
import argparse
from pathlib import Path
import logging
from src.data import load_data
from src.dataset import IndexedImage2TranscriptomeDataset
from src.models import PretrainedResNet, PretrainedEfficientNet
from src.training import train, inference
from src.eval import calculate_results_per_batch, adata_results
from src.plotting import scatter_top_genes, plot_example_images
sns.set_style('whitegrid')

parser = argparse.ArgumentParser()
parser.add_argument('--seeds', default=[0, 1, 2], nargs='+')
parser.add_argument('--suffix', default=None)
parser.add_argument('--data_dir', default='/mlbio_scratch/vinas/iris-data')
parser.add_argument('--out_dir', default='/mlbio_scratch/vinas/image2transcriptome/results')
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
parser.add_argument('--device', default=0, type=int)
parser.add_argument('--random_split', action='store_true')
parser.add_argument('--train_experiments_prefix', default='NG')
parser.add_argument('--val_experiments_prefix', default=None)
parser.add_argument('--train_augmentations', default='albumentations')
parser.add_argument('--test_augmentations', default=None)
parser.add_argument('--shuffle_size', default=None, type=int)
parser.add_argument('--gaussblur_size', default=None, type=int)
parser.add_argument('--bs', default=256, type=int)
parser.add_argument('--num_workers', default=8, type=int)
parser.add_argument('--model_name', default='resnet')
parser.add_argument('--emb_dim', default=512, type=int)
parser.add_argument('--backbone_size', default='18', type=str)
parser.add_argument('--epochs', default=50, type=int)
parser.add_argument('--lr', default=1e-3, type=float)
parser.add_argument('--wd', default=1e-4, type=float)
parser.add_argument('--save_model', action='store_true')
parser.add_argument('--discard_spikein', action='store_true')
parser.set_defaults(random_split=False)
parser.set_defaults(minmax_normalize=False)
parser.set_defaults(apply_segmentation_mask=False)
parser.set_defaults(save_model=False)
parser.set_defaults(discard_spikein=False)
args = parser.parse_args()


def train_and_store_results(args, adata_subset, dataloaders, seed, out_dir, obsm_key = None):    
    # Instantiate model
    logging.info('... instantiating model')
    if obsm_key is None:
        out_dim = adata_subset.shape[-1]
    else:
        out_dim = adata_subset.obsm[obsm_key].shape[-1]

    # torch.cuda.empty_cache()
    device = torch.device(f'cuda:{args.device}' if torch.cuda.is_available() else 'cpu')

    if args.model_name == 'EfficientNet':
        model =  PretrainedEfficientNet(out_dim=out_dim, cnn_channels_adapter=len(channels_dict))
    else:
        model = PretrainedResNet(args.backbone_size, embed_dim=args.emb_dim, out_dim=out_dim, cnn_channels_adapter=len(channels_dict))
    
    model = model.to(device)
    
    # Create output dir and store arguments
    results_dir = f'{out_dir}/seed_{seed}'
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    with open(f'{results_dir}/args.txt', 'w') as f:
        json.dump(args.__dict__, f, indent=2)

    # Logging
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    logging.FileHandler(f'{out_dir}/outputs.log'),
    log.addHandler(logging.FileHandler(f'{results_dir}/outputs.log'))      # set the new handler
    log.addHandler(logging.StreamHandler())      # set the new handler

    # Train model
    logging.info('... training model')
    model = train(model,
                  dataloaders,
                  device=device,
                  num_epochs=args.epochs,
                  lr=args.lr,
                  wd=args.wd,
                  seed=seed)
    if args.save_model:
        torch.save(model.state_dict(), f'{results_dir}/model.pth')
    
    # Run inference on test datasets (including different inference augmentations)
    logging.info('... running inference on test sets')
    for key, dataloader in dataloaders.items():
        if not key.startswith('test'):
            continue
        logging.info(f'Running inference checks for {key}')
        preds, gts, _ = inference(model, dataloader, device)
        df_check = calculate_results_per_batch(adata_subset[test_idxs], preds, gts)
        df_check.to_csv(f'{results_dir}/{key}_results.csv')

        # Plot results
        plt.figure(figsize=(16, 6))
        scatter_top_genes(preds, gts, df_check, adata_subset[test_idxs], k=10)
        plt.tight_layout()
        plt.savefig(f'{results_dir}/{key}_scatter_top_genes.pdf', bbox_inches='tight')

        if key == 'test':
            # Store results for the main test set
            adata_pred = adata_results(adata_subset[test_idxs], preds, gts)
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
    model_string = args.model_name
    out_dir = f'{args.out_dir}/{args.dataset}{ct_string}_{channels_string}{segm_string}{split_string}_{model_string}{suffix_string}'
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

    # Choose resize size based on model
    resize_size = 256
    if args.model_name == 'EfficientNet':
        resize_size = 224
    logging.info(f'Resize size: {resize_size}')

    # Splitting data by experiment ID
    logging.info('... splitting data')
    df = adata_subset.obs

    for seed in args.seeds:
        logging.info(f'Setting random seed: {seed}')
        np.random.seed(int(seed))

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

        # Create datasets
        logging.info(f'Image2transcriptome model. Channels: {",".join(channels)}')
        logging.info('... creating Pytorch datasets')
        obsm_key = None
        # multiple_channels = len(channels_dict) > 1
        # img_list = channels_dict[args.channel]
        train_dataset = IndexedImage2TranscriptomeDataset(adata_subset,
                                                        # img_list=img_list,
                                                        channels_dict=channels_dict,
                                                        idxs=train_idxs,
                                                        obsm_key=obsm_key,
                                                        augmentation=args.train_augmentations,
                                                        gaussblur_size=args.gaussblur_size if args.train_augmentations == 'gaussblur' else None,
                                                        shuffle_size=args.shuffle_size if args.train_augmentations == 'shuffle' else None,
                                                        repeat_channels=False)
        val_dataset = IndexedImage2TranscriptomeDataset(adata_subset,
                                                        # img_list=img_list,
                                                        channels_dict=channels_dict,
                                                        idxs=val_idxs,
                                                        obsm_key=obsm_key,
                                                        repeat_channels=False)
        test_dataset = IndexedImage2TranscriptomeDataset(adata_subset,
                                                        # img_list=img_list,
                                                        channels_dict=channels_dict,
                                                        idxs=test_idxs,
                                                        obsm_key=obsm_key,
                                                        augmentation=args.test_augmentations,
                                                        gaussblur_size=args.gaussblur_size if args.test_augmentations == 'gaussblur' else None,
                                                        shuffle_size=args.shuffle_size if args.test_augmentations == 'shuffle' else None,
                                                        repeat_channels=False)
        datasets = {"train": train_dataset,
                    "val": val_dataset,
                    "test": test_dataset}
        dataloaders = {
            "train": DataLoader(
                datasets["train"],
                batch_size=args.bs,
                shuffle=True,
                num_workers = args.num_workers
            ),
            "val": DataLoader(
                datasets["val"],
                batch_size=args.bs,
                shuffle=False,
                num_workers = args.num_workers
            ),
            "test": DataLoader(
                datasets["test"],
                batch_size=args.bs,
                shuffle=False,
                num_workers = args.num_workers
            ),
        }

        # Include test datasets for inference checks
        for shuffle_size in [1, 4, 8, 16, 32, 64]:
            shuffle_dataset = IndexedImage2TranscriptomeDataset(adata_subset,
                                                        channels_dict=channels_dict,
                                                        idxs=test_idxs,
                                                        obsm_key=obsm_key,
                                                        augmentation='shuffle',
                                                        shuffle_size=shuffle_size,
                                                        repeat_channels=False)
            datasets[f'test_shuffle_{shuffle_size}'] = shuffle_dataset
            dataloaders[f'test_shuffle_{shuffle_size}'] = DataLoader(shuffle_dataset,
                batch_size=args.bs,
                shuffle=False,
                num_workers = args.num_workers
            )
        for gaussblur_size in  [3, 9, 27, 81, 151]:
            gaussblur_dataset = IndexedImage2TranscriptomeDataset(adata_subset,
                                                        channels_dict=channels_dict,
                                                        idxs=test_idxs,
                                                        obsm_key=obsm_key,
                                                        augmentation='gaussblur',
                                                        gaussblur_size=gaussblur_size,
                                                        repeat_channels=False)
            datasets[f'test_gaussblur_{gaussblur_size}'] = gaussblur_dataset
            dataloaders[f'test_gaussblur_{gaussblur_size}'] = DataLoader(gaussblur_dataset,
                batch_size=args.bs,
                shuffle=False,
                num_workers = args.num_workers
            )
        
        # Plot some example images
        logging.info('... plotting example images')
        train_dataset_noaug = IndexedImage2TranscriptomeDataset(adata_subset,
                                                        # img_list=img_list,
                                                        channels_dict=channels_dict,
                                                        idxs=train_idxs,
                                                        obsm_key=obsm_key,
                                                        augmentation=None,
                                                        repeat_channels=False)
        
        # Plot some example images with and without augmentation
        for channel_idx in range(len(channels_dict)):
            plot_example_images(train_dataset_noaug, train_dataset, out_dir, key='train', channel_idx=0)
        
            for key, dataset in datasets.items():
                if not key.startswith('test_'):
                    continue
                for channel_idx in range(len(channels_dict)):
                    plot_example_images(test_dataset, dataset, out_dir, key=key, channel_idx=channel_idx)

        # Train and store results for each seed
        train_and_store_results(args, adata_subset, dataloaders, seed, out_dir)