import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torchvision.models import resnet18, resnet34, resnet50, resnet101, resnet152

class PretrainedResNet(nn.Module):
    """
    ResNet backbone
    """
    def __init__(self, backbone_size, embed_dim, out_dim, cnn_channels_adapter=None) -> None:
        super().__init__()

        self.embed_dim = embed_dim
        self.backbone_size = backbone_size
        self.out_dim = out_dim
        if self.backbone_size == "18":
            self.backbone = resnet18(weights="DEFAULT")
        elif self.backbone_size == "34":
            self.backbone = resnet34(weights="DEFAULT")
        elif self.backbone_size == "50":
            self.backbone = resnet50(weights="DEFAULT")
        elif self.backbone_size == "101":
            self.backbone = resnet101(weights="DEFAULT")
        elif self.backbone_size == "152":
            self.backbone = resnet152(weights="DEFAULT")
        else:
            raise ValueError(
                "Invalid backbone size. Must be one of ['18', '34', '50', '101', '152']"
            )
        if cnn_channels_adapter is not None:
            # Convolutional layer to adapt the number of input channels
            self.backbone.conv1 = nn.Conv2d(
                cnn_channels_adapter, 64, kernel_size=(7, 7), stride=(2, 2), padding=(3, 3), bias=False
            )
        self.backbone.fc = (
            nn.Identity()
        )  # remove the last fully connected layer by setting it to an identity function
        self.linear = nn.Linear(self.embed_dim, self.out_dim)

    def forward(self, x):
        x = self.backbone(x)
        x = self.linear(x)
        return x

class PretrainedEfficientNet(nn.Module):
    def __init__(self, out_dim, embed_dim=1000, cnn_channels_adapter=None):
        super().__init__()

        self.embed_dim = embed_dim
        self.out_dim = out_dim
        self.backbone = torch.hub.load('NVIDIA/DeepLearningExamples:torchhub', 'nvidia_efficientnet_b0', pretrained=True)
        for param in self.backbone.parameters():  # Freeze the backbone
            param.requires_grad = False

        if cnn_channels_adapter is not None:
            # Convolutional layer to adapt the number of input channels
            self.backbone.stem.conv = nn.Conv2d(
                cnn_channels_adapter, 32, kernel_size=(3, 3), stride=(2, 2), padding=(1, 1), bias=False
            )
            self.backbone.stem.conv.requires_grad = True
        self.linear = nn.Linear(self.embed_dim, self.out_dim)

    def forward(self, x):
        x = self.backbone(x)
        x = self.linear(x)
        return x

"""
from lightly.models import utils
from lightly.models.modules import MAEDecoderTIMM, MaskedVisionTransformerTIMM, MAEBackbone
from lightly.transforms import MAETransform

class MAE(nn.Module):
    def __init__(self, vit, out_dim=2000):
        super().__init__()
        decoder_dim = 512
        self.mask_ratio = 0.0
        self.patch_size = vit.patch_embed.patch_size[0]
        self.backbone = MaskedVisionTransformerTIMM(vit=vit)
        self.sequence_length = self.backbone.sequence_length
        self.decoder = MAEDecoderTIMM(
            num_patches=vit.patch_embed.num_patches,
            patch_size=self.patch_size,
            embed_dim=vit.embed_dim,
            decoder_embed_dim=decoder_dim,
            decoder_depth=1,
            decoder_num_heads=16,
            mlp_ratio=4.0,
            proj_drop_rate=0.0,
            attn_drop_rate=0.0,
        )
        self.head = nn.Linear(1024, out_dim)

    def forward_encoder(self, images, idx_keep=None):
        return self.backbone.encode(images=images, idx_keep=idx_keep)

    def forward_decoder(self, x_encoded, idx_keep, idx_mask):
        # build decoder input
        batch_size = x_encoded.shape[0]
        x_decode = self.decoder.embed(x_encoded)
        x_masked = utils.repeat_token(
            self.decoder.mask_token, (batch_size, self.sequence_length)
        )
        x_masked = utils.set_at_index(x_masked, idx_keep, x_decode.type_as(x_masked))

        # decoder forward pass
        x_decoded = self.decoder.decode(x_masked)

        # predict pixel values for masked tokens
        x_pred = utils.get_at_index(x_decoded, idx_mask)
        x_pred = self.decoder.predict(x_pred)
        return x_pred

    def freeze_encoder(self):
        for param in self.backbone.parameters():
            param.requires_grad = False
        print("Encoder has been frozen.")
        
    def freeze_encoder_partial(self):
        for param in self.backbone.parameters():
            param.requires_grad = False
        for name, param in self.backbone.named_parameters():
            if 'vit.blocks.23' in name:
                param.requires_grad = True
        print("Encoder has been partially frozen.")

    def forward(self, images):
        batch_size = images.shape[0]
        idx_keep, idx_mask = utils.random_token_mask(
            size=(batch_size, self.sequence_length),
            mask_ratio=self.mask_ratio,
            device=images.device,
        )
        x_encoded = self.forward_encoder(images=images, idx_keep=idx_keep)
        feature = torch.mean(x_encoded, 1)
        return self.head(feature)
"""