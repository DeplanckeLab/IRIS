import h5py
import json
import pandas as pd
import anndata as ad
import skimage as ski
from skimage.util import montage
import numpy as np
import cv2


def get_exp_data(gene_exp: bytes) -> dict:
    decode_data = gene_exp.decode("utf-8")
    decode_data = decode_data.replace("\'","\"")
    exp_dict = json.loads(decode_data)
    return exp_dict
    
def expdf_from_h5(h5_path):
    exp_list = []
    with h5py.File(h5_path, 'r') as dataset:
        cell_ids = list(dataset.keys())
        for cell_id in cell_ids:
            cell_ds = dataset[cell_id]
            if "exp_data" in cell_ds and cell_ds.attrs['num_cells'] == 1:
                try:
                    exp_data = cell_ds['exp_data/exp_data'][()]
                except:
                    print("no exp dataset:", cell_id)
                else:
                    exp_dict = get_exp_data(exp_data)
                    exp_series = pd.Series(exp_dict)
                    exp_series.name = cell_id
                    exp_list.append(exp_series)

    exp_df = pd.DataFrame(exp_list)
    return exp_df

def metadf_from_h5(h5_path):
    meta_list = []
    with h5py.File(h5_path, 'r') as dataset:
        cell_ids = list(dataset.keys())
        for cell_id in cell_ids:
            cell_ds = dataset[cell_id]

            if "exp_data/exp_data" in cell_ds and cell_ds.attrs['num_cells'] == 1:
                attr_list = [('cell_id', cell_id)]

                for attr_name in cell_ds.attrs.keys():
                    attr_value = cell_ds.attrs[attr_name] 
                    attr_list.append((attr_name, attr_value))

                for attr_name in cell_ds['exp_data'].attrs.keys():
                    attr_value = cell_ds['exp_data'].attrs[attr_name] 
                    attr_list.append((attr_name, attr_value))

                attr_dict = {at_k: at_v for at_k, at_v in attr_list}
                meta_list.append(pd.Series(attr_dict))
    meta_df = pd.DataFrame(meta_list) 
    meta_df.set_index('cell_id', inplace=True)
    return meta_df

def annd_from_h5(h5_path):
    exp_df = expdf_from_h5(h5_path)
    meta_df = metadf_from_h5(h5_path)
    exp_df = exp_df.fillna(0)
    adata = ad.AnnData(exp_df, exp_df.index.to_frame(), exp_df.columns.to_frame())
    adata.obs = meta_df
    print(adata)
    return adata


def get_image_by_cid(cell_id_list, channel, h5_path):
    img_dict = {}
    with h5py.File(h5_path, 'a') as dataset:
        for cell_id in cell_id_list:
            cell_ds = dataset[str(cell_id)]
            foc_plane = str((cell_ds.attrs['focal_plane']))
            img = (cell_ds[channel][foc_plane][()]/255)
            #seg = cell_ds['seg'][foc_plane][()]
            img = img.astype('uint8')
            #img = img[75:175,75:175].copy()
            #seg = seg[75:175,75:175].copy()
            img = img - np.min(img)
            img = img/np.max(img)

            rgba_img = np.zeros((img.shape[0], img.shape[1], 4))
            rgba_img[..., 0] = img
            rgba_img[..., 1] = img
            rgba_img[..., 2] = img
            rgba_img[..., 3] = np.ones(img.shape)#seg

            img_dict[cell_id] = rgba_img
    return img_dict

def gen_montage_multich(cell_id_list, channel, h5_path, anno_list = None, ch_sep=False, verbose=True):
    img_dict = {ch: [] for ch in channel.keys()}
    with h5py.File(h5_path, 'r') as dataset:
        if anno_list == 'cid':
            anno_list = cell_id_list
        if anno_list and len(anno_list) != len(cell_id_list):
            return None

        for l_id, cell_id in enumerate(cell_id_list):
            for ch in channel:
                cell_ds = dataset[str(cell_id)]
                foc_plane = str((cell_ds.attrs['focal_plane']))
                img = cell_ds[ch][foc_plane][()]/256
                img = img.astype('uint8')

                if channel[ch]["norm"] and not channel[ch]["col"] == "grey":
                    img = norm_minmax(img, verbose=verbose)
                if channel[ch]["col"] == "grey":
                    c_mean = np.mean(img)
                    adj_fact = 150 / c_mean if c_mean > 0 else 1
                    img = np.clip(img*adj_fact, 0, 255)/255#.astype(np.uint8)

                img_dict[ch].append(img)

    for ch in img_dict.keys():
        img_mon = montage(img_dict[ch], grid_shape = (1, len(cell_id_list)))
        if not channel[ch]["norm"]:
            img_mon = norm_minmax(img_mon, verbose=verbose)
        img_dict[ch] = img_mon

    if ch_sep:
        for ch in img_dict.keys():
            img_mon = img_dict[ch]
            rgba_img = np.zeros((img_mon.shape[0], img_mon.shape[1], 3))
            img_color = channel[ch]["col"]
            if img_color == "red":
                rgba_img[..., 0] = img_mon
                rgba_img[..., 1] = img_mon * 0.2
                rgba_img[..., 2] = img_mon * 0.2
            
            if img_color == "green":
                rgba_img[..., 1] = img_mon
                rgba_img[..., 0] = img_mon * 0.2
                rgba_img[..., 2] = img_mon * 0.2
            
            if img_color == "blue":
                # Mix in some other colors to make the blue brighter
                rgba_img[..., 0] = img_mon * 0.2
                rgba_img[..., 1] = img_mon * 0.2
                rgba_img[..., 2] = img_mon

            if img_color == "grey":
                rgba_img[..., 0] = img_mon
                rgba_img[..., 1] = img_mon
                rgba_img[..., 2] = img_mon
                
            img_dict[ch] = rgba_img

        img_ar_all = [img for img in img_dict.values()]
        img_montage = montage(img_ar_all, grid_shape = (len(channel),1), channel_axis=3)
        return img_montage

    else:
        shape_x, shape_y = list(img_dict.values())[0].shape
        rgba_img = np.zeros((shape_x, shape_y, 3))
        for ch, _ in channel.items():
            img_color = channel[ch]["col"]
            img_mon = img_dict[ch]
            if img_color == "red":
                rgba_img[..., 0] = img_mon
            
            if img_color == "green":
                rgba_img[..., 1] = img_mon
            
            if img_color == "blue":
                rgba_img[..., 2] = img_mon

        return rgba_img

def norm_minmax(img, verbose=True):
    if verbose:
        print(f'Image min {np.min(img)} max {np.max(img)}')
    img_max = np.max(img)
    if img_max > 0:
        if verbose:
            print('correcting')
        img = (img - np.min(img))/(np.max(img)-np.min(img))
    return img

print(f"[sc_utils] Using scikit-image version: {ski.__version__}")
print(f"[sc_utils] Using cv2 version: {cv2.__version__}")

