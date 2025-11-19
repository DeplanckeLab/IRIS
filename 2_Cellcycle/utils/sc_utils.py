import h5py
import json
import pandas as pd
import anndata as ad
import skimage as ski
from skimage.util import montage
from skimage import exposure
import numpy as np
import cv2
from pandas.api.types import is_numeric_dtype
import adjustText
from adjustText import adjust_text
from matplotlib.patches import FancyArrowPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


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

def add_obs_as_attr(obs_df, h5_path):
    with h5py.File(h5_path, 'a') as dataset:
        for index, row in obs_df.iterrows():
            cell_ds = dataset[index]
            for meta_item in row.items():
                attr_name, attr_value = meta_item
                cell_ds.attrs[attr_name] = attr_value


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
                    #clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
                    #img = clahe.apply(img)
                    c_mean = np.mean(img)
                    adj_fact = 150 / c_mean if c_mean > 0 else 1
                    img = np.clip(img*adj_fact, 0, 255)/255#.astype(np.uint8)
                    #img = cv2.equalizeHist(img)/255
                    #img = norm_minmax(img)

#                if anno_list:
#                    cell_anno = anno_list[l_id]
#                    cv2.putText(img,cell_anno, (0,15), cv2.FONT_HERSHEY_PLAIN, 1, (1,1,1),1, 1)
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

def mask_img(img, mask):
    mask = (mask > 0).astype(np.uint8)
    masked_img = img*mask
    return masked_img

def mask_img_wdil(img, mask):
    mask = (mask > 0).astype(np.uint8)
    #print("Size before:", np.sum(mask))
    kernel = np.ones((3, 3), dtype=np.uint8)
    mask = cv2.dilate(mask, kernel, iterations=2)
    #print("Size after:", np.sum(mask))
    masked_img = img*mask
    return masked_img

#def bump_detect(mask):
#    mask = (mask > 0).astype(np.uint8)
#    kernel = np.ones((5,5), np.uint8)
#    eroded = cv2.erode(mask, kernel, iterations=5)
#    dilated = cv2.dilate(eroded, kernel, iterations=5)
#    difference  = mask - dilated
#    bump_px = np.sum(difference)
#    return bump_px


def bump_detect(mask):
    mask = (mask > 0).astype(np.uint8)
    sk = ski.morphology.skeletonize(mask == 1)
    sk_len = np.sum(sk)
    return sk_len

def meas_bump(h5_path):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            foc_plane = str((cell_ds.attrs['focal_plane']))
            seg = cell_ds['seg'][foc_plane][()]
            bump_px = bump_detect(seg)
            out_dict[cell_id] = bump_px
    return out_dict

def circ_func(mask):
    mask = (mask > 0).astype(np.uint8)
    perim = ski.measure.perimeter(mask, neighborhood=4)
    circ = ((4*3.1415)*np.sum(mask))/(perim**2)
    return circ
    
def meas_circ(h5_path):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            foc_plane = str((cell_ds.attrs['focal_plane']))
            seg = cell_ds['seg'][foc_plane][()]
            circ = circ_func(seg)
            out_dict[cell_id] = circ
    return out_dict
    

def measure_channel(h5_path, channel, mask_func=False):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            now_break = False
            if channel in cell_ds:
                foc_plane = str((cell_ds.attrs['focal_plane']))
                img = cell_ds[channel][foc_plane][()]
                seg = cell_ds['seg'][foc_plane][()]/255
                
                if mask_func:
                    img = mask_func(img, seg)
             
                cell_area = np.sum(seg)
                if cell_area > 0:
                    intensity = np.sum(img) / cell_area
                    out_dict[cell_id] = intensity
                #now_break = True
            else:
                out_dict[cell_id] = 0
            if now_break:
                break
    return out_dict

def measure_channel_max(h5_path, channel):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            now_break = False
            if channel in cell_ds:
                foc_plane = str((cell_ds.attrs['focal_plane']))
                image = cell_ds[channel][foc_plane][()]
                seg = cell_ds['seg'][foc_plane][()]
                intensity = np.max(image)
                out_dict[cell_id] = intensity
                #now_break = True
            else:
                out_dict[cell_id] = 0
            if now_break:
                break
    return out_dict

def measure_size(h5_path):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            now_break = False
            if 'seg' in cell_ds:
                foc_plane = str((cell_ds.attrs['focal_plane']))
                seg = cell_ds['seg'][foc_plane][()]/255
                intensity = np.sum(seg)
                out_dict[cell_id] = intensity
            else:
                out_dict[cell_id] = 0
            if now_break:
                break
    return out_dict

def bump_detect(mask):
    mask = (mask > 0).astype(np.uint8)
    sk = ski.morphology.skeletonize(mask == 1)
    sk_len = np.sum(sk)
    return sk_len

def meas_bump(h5_path):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            foc_plane = str((cell_ds.attrs['focal_plane']))
            seg = cell_ds['seg'][foc_plane][()]
            bump_px = bump_detect(seg)
            out_dict[cell_id] = bump_px
    return out_dict


def measure_cyto_mem_int(h5_path, channel):
    out_dict = {}
    with h5py.File(h5_path, 'r') as dataset:
        now_break = 0
        #for cell_id in ['16144']:#dataset.keys():
        for cell_id in dataset.keys():
            cell_ds = dataset[str(cell_id)]
            if channel in cell_ds:
                kernel_out = np.ones((11,11), np.uint8)
                kernel_in = np.ones((21,21), np.uint8)
                foc_plane = str((cell_ds.attrs['focal_plane']))
                foc_638 = str((cell_ds.attrs['focal_plane'] + 1))
                if foc_638 not in cell_ds[channel]:
                    foc_638 = foc_plane
                image = cell_ds[channel][foc_plane][()]
                seg = cell_ds['seg'][foc_638][()]/255
                seg_cyto = cv2.erode(seg,kernel_in,iterations=1)
                seg_ero = cv2.dilate(seg,kernel_out,iterations=1)

                seg_mem = seg_ero - seg_cyto

                image_mem = image * seg_mem
                image_cyto = image * seg_cyto

                int_mem = np.sum(image_mem) / np.sum(seg_mem)
                int_cyto = np.sum(image_cyto) / np.sum(seg_cyto)

                int_ratio = int_mem/int_cyto

                out_dict[cell_id] = int_ratio
                #now_break = now_break + 1
            else:
                out_dict[cell_id] = 0
            if now_break == 10:
                break
    #return image, seg, seg_cyto, seg_ero, image_mem, image_cyto
    return out_dict

def measure_cyto_mem2(h5_path, channel):
    out_dict = {}
    with h5py.File(h5_path, 'a') as dataset:
        now_break = 0
        #for cell_id in ['16144']:#dataset.keys():
        for cell_id in ['16144']:
            cell_ds = dataset[str(cell_id)]
            if channel in cell_ds:
                foc_plane = str((cell_ds.attrs['focal_plane']))
                foc_638 = str((cell_ds.attrs['focal_plane']))
                if foc_638 in cell_ds[channel]:
                    foc_plane = foc_638
                image = cell_ds[channel][foc_plane][()]
                x_len, y_len = image.shape
                y, x = np.indices(image.shape)
                #print(x,y)
                center = (int(x_len/2),int(y_len/2))
                distance_from_center = np.sqrt((x - center[1])**2 + (y - center[0])**2)
                radial_bins = np.arange(0,distance_from_center.max(), 10)
                radial_profile = [
                    image[(distance_from_center >= r) & (distance_from_center < r + 1)].mean()
                    for r in radial_bins
                        ]
                return image, distance_from_center, radial_profile


def add_image_by_highlight(ax, highlight_dots, high_df, channel, ds_dir):
    texts = []
    x_coords = []
    y_coords = []
    #for collection in ax.collections:
    #    offsets = collection.get_offsets()
    #    [x_coords.append(x) for x in offsets[:,0]] 
    #    [y_coords.append(y) for y in offsets[:,1]] 
        
    for cell_id, x, y in highlight_dots:
        texts.append(ax.text(x, y, cell_id, ha='center', va='center', fontsize = 'small')) # explode_radius 

    adh_texts = adjust_text(texts, ax=ax, x=x_coords, y=y_coords, expand=(2, 2), arrowprops=dict(arrowstyle='->', color='red'), 
                            time_lim=10, max_move=(50,50), only_move={"text": "xy", "static": "xy", "explode": "xy", "pull": "xy"})
    image_dict = get_image_by_cid(high_df.index.to_list(), channel, ds_dir)

    if adh_texts:
        for adj_text in adh_texts:
            for label in adj_text:
                if isinstance(label, FancyArrowPatch): continue
                x, y = label.get_position()
                index = int(label.get_text())
                label.set_text(None)
                image = image_dict[str(index)]
                offset_image = OffsetImage(image, zoom=0.1)
                ab = AnnotationBbox(offset_image, (x, y), frameon=False, alpha = 0.5)
                ax.add_artist(ab)

def highligh_dots(ax, df, offsets, ds_dir):
    highlight_dots = []
    
    for off_index in offsets:
        for x_coord, y_coord in ax.collections[off_index].get_offsets():
            tmp_frame = df
    
            # X axis
            axis_title = ax.get_xlabel()
            if is_numeric_dtype(tmp_frame[axis_title].dtype):
                tmp_frame = tmp_frame[np.isclose(tmp_frame[axis_title], x_coord)]
            else:
                ticks = [tick.get_text() for tick in ax.get_xticklabels()]
                axis_lookup = {num:label for num, label in enumerate(ticks)}
                closest_value = min(axis_lookup.keys(), key=lambda x:abs(x-x_coord))
                tmp_frame = tmp_frame[tmp_frame[axis_title] == axis_lookup[closest_value]]
                        
            # Y axis
            axis_title = ax.get_ylabel()
            if is_numeric_dtype(tmp_frame[axis_title].dtype):
                tmp_frame = tmp_frame[np.isclose(tmp_frame[axis_title], y_coord)]
            else:
                ticks = [tick.get_text() for tick in ax.get_yticklabels()]
                axis_lookup = {num:label for num, label in enumerate(ticks)}
                closest_value = min(axis_lookup.keys(), key=lambda x:abs(x-y_coord))
                tmp_frame = tmp_frame[tmp_frame[axis_title] == axis_lookup[closest_value]]                     
        
            if len(tmp_frame.index.tolist()) != 0:
                cell_id = tmp_frame.index.tolist()[0]
                highlight_dots.append((cell_id, x_coord, y_coord))
    
    add_image_by_highlight(ax, highlight_dots, df, '638', ds_dir)



print(f"[sc_utils] Using scikit-image version: {ski.__version__}")
print(f"[sc_utils] Using cv2 version: {cv2.__version__}")
print(f"[sc_utils] Using adjust_text version: {adjustText.__version__}")

