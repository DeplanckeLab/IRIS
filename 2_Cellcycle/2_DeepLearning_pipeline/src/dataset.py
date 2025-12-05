import torch
from torch.utils.data import Dataset
# from torchvision import transforms as T
import albumentations as A
import numpy as np
import random
import cv2


class MinMaxNormalize(A.ImageOnlyTransform):
    def __init__(self, always_apply=True, p=1.0):
        super(MinMaxNormalize, self).__init__(always_apply=always_apply, p=p)

    def apply(self, img, **params):
        min_val = img.min()
        max_val = img.max()
        if max_val > min_val:
            img = (img - min_val) / (max_val - min_val)
        else:
            img = img * 0.0
        return img.astype(np.float32)

class RandomShuffle(A.ImageOnlyTransform):
    def __init__(self, always_apply=True, p=1.0, patch_size=1):
        super(RandomShuffle, self).__init__(always_apply=always_apply, p=p)
        self.patch_size = (patch_size, patch_size)

    def apply(self, img, **params):
        image = img
        h, w, n_channels = image.shape

        # Determine number of patches
        patch_h, patch_w = self.patch_size
        patches = []

        # Crop image to be a multiple of patch size
        h_crop = h - h % patch_h
        w_crop = w - w % patch_w
        image = image[:h_crop, :w_crop]

        # Extract patches
        for y in range(0, h_crop, patch_h):
            for x in range(0, w_crop, patch_w):
                patch = image[y:y+patch_h, x:x+patch_w]
                patches.append(patch)

        # Shuffle patches
        random.shuffle(patches)

        # Reconstruct the image
        shuffled_image = np.zeros_like(image)
        i = 0
        for y in range(0, h_crop, patch_h):
            for x in range(0, w_crop, patch_w):
                shuffled_image[y:y+patch_h, x:x+patch_w] = patches[i]
                i += 1

        return shuffled_image.astype(np.float32)

class GaussBlur(A.ImageOnlyTransform):
    def __init__(self, always_apply=True, p=1.0, kernel_size=5):
        super(GaussBlur, self).__init__(always_apply=always_apply, p=p)
        self.kernel_size = (kernel_size, kernel_size)

    def apply(self, img, **params):
        blur = cv2.blur(img, self.kernel_size)
        return blur

class IndexedImage2TranscriptomeDataset(Dataset):
    def __init__(self, adata, channels_dict, idxs, obsm_key=None, resize_size=256, repeat_channels=True, augmentation=None, shuffle_size=None, gaussblur_size=None):
        self.adata = adata
        self.obsm_key = obsm_key
        self.channels_dict = channels_dict
        self.channel_names = sorted(list(channels_dict.keys()))
        self.idxs = idxs
        transforms = [A.Resize(height=resize_size, width=resize_size, p=1.0)]
        if augmentation == 'albumentations':
            # Define transform
            augmentations = [
                A.HorizontalFlip(p=0.5),
                A.RandomBrightnessContrast(p=0.2, brightness_limit = (-0.01, 0.01), contrast_limit =  (-0.01, 0.01)),
                A.GaussNoise(std_range=(0.005, 0.01), p=1.0),
                A.Rotate(p=0.5),
            ]
            transforms.extend(augmentations)
        elif augmentation == 'random_rotation':
            transforms.append(A.Rotate(p=0.5))
        elif augmentation == 'shuffle':
            assert shuffle_size is not None
            transforms.append(RandomShuffle(patch_size=shuffle_size, p=1.0))
        elif augmentation == 'gaussblur':
            assert gaussblur_size is not None
            transforms.append(GaussBlur(kernel_size=gaussblur_size, p=1.0))

        # MinMax normalization and toTensor always applied
        transforms.append(MinMaxNormalize())
        transforms.append(A.ToTensorV2())
        self.transform = A.Compose(transforms)

        self.repeat_channels = repeat_channels

    def __len__(self):
        return len(self.idxs)
        
    def __getitem__(self, i):
        idx = self.idxs[i]
        image = np.array([self.channels_dict[channel][idx] for channel in self.channel_names])  # Shape: (C, H, W)
        # image = self.img_list[idx]
        image = image / 65535.0
        image = image.astype(np.float32)
        image = self.transform(image=image.transpose(1, 2, 0))  # Transpose to (H, W, C) for albumentations

        # image = image / 255.0
        """
        if self.repeat_channels and image.shape[0] == 1:
            image = image.repeat((3, 1, 1))
        """
        # image = image.repeat((3, 1, 1)).float()
        # image = image.float()
        if self.obsm_key is None:
            y = self.adata[idx].X
        else:
            y = self.adata[idx].obsm[self.obsm_key]
        y = torch.tensor(y).ravel().float()
        return image['image'], y