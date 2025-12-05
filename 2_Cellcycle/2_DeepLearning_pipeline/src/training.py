import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms as T
from torchvision.models import resnet18, resnet34, resnet50, resnet101, resnet152
from collections import OrderedDict
from tqdm import tqdm
import numpy as np
import logging

class EarlyStopper:    
    def __init__(self, patience, delta=0.0):
        self.patience = patience
        self.counter = 0
        self.min_val_loss = float("inf")
        self.delta = delta
        self.best_model_state = None

    def early_stop(self, val_loss, model):
        if self.best_model_state is None:
            self.best_model_state = model.state_dict()
        if val_loss < self.min_val_loss:
            self.min_val_loss = val_loss
            self.counter = 0
            self.best_model_state = model.state_dict()
        elif val_loss > self.delta + self.min_val_loss:
            self.counter += 1
            if self.counter >= self.patience:
                return True
        return False
    

def train(model, dataloaders, device, num_epochs=30, lr=1e-3, wd=1e-4, seed=0, patience=10, step_size=10):
    torch.manual_seed(seed)
    
    # Train
    early_stopper = EarlyStopper(patience=patience)
    optimizer = optim.Adam(
        lr=lr,
        params=model.parameters(),
        weight_decay=wd
    )
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=step_size, gamma=0.1)

    pbar_epoch = tqdm(range(num_epochs), position=0, leave=True)
    phases = ['train']
    if len(dataloaders['val'].dataset) > 0:
        phases.append('val')
    for epoch in pbar_epoch:
        for phase in phases:
            if phase == "train":
                model.train()
            else:
                model.eval()

            running_loss = 0.0
            pbar = tqdm(dataloaders[phase], position=0, leave=True)
            for batch in pbar:
                images, ys = batch
                images = images.to(device)
                ys = ys.to(device)
                optimizer.zero_grad()

                with torch.set_grad_enabled(phase == "train"):
                    y_pred = model(images)
                    loss = torch.nn.MSELoss()(ys, y_pred)

                    if phase == "train":
                        loss.backward()
                        optimizer.step()
                pbar.set_postfix(loss=f"{phase} batch loss: {loss.item():.4f}")
                running_loss += loss.item() * images.size(0)

            if phase == "train":
                if scheduler is not None:
                    scheduler.step()

            epoch_loss = running_loss / len(dataloaders[phase].dataset)
            logging.info(f"{phase} loss: {epoch_loss:.4f}")

            if phase == "val" and early_stopper.early_stop(epoch_loss, model):
                logging.info(f"Early stopping at epoch {epoch}")
                model.load_state_dict(early_stopper.best_model_state)
                break
        else:
            continue
        break
    return model

def inference(model, dataloader, device, eval_on = 'test', store_intermediate=False):
    if store_intermediate:
        # Utilities to get latent features (prior to output layer) as per Johannes' request
        def remove_all_forward_hooks(model: torch.nn.Module) -> None:
            for name, child in model._modules.items():
                if child is not None:
                    if hasattr(child, "_forward_hooks"):
                        child._forward_hooks = OrderedDict()
                    remove_all_forward_hooks(child)
        activation = {'fc': []}
        def get_activation(name):
            def hook(model, input, output):
                activation[name].extend(output.detach().cpu().numpy())
            return hook
        remove_all_forward_hooks(model)
        model.backbone.fc.register_forward_hook(get_activation('fc'))
        
    # Inference
    model.eval()
    pbar = tqdm(dataloader, position=0, leave=True)
    preds = []
    gts = []
    with torch.no_grad():
        for batch in pbar:
            images, ys = batch
            images = images.to(device)
            ys = ys.to(device)
            y_pred = model(images)
            gts.extend(ys.cpu().numpy())
            preds.extend(y_pred.cpu().numpy())
    preds = np.array(preds)
    gts = np.array(gts)
    output_features = None
    if store_intermediate:
        output_features = np.stack(activation['fc'])
    return preds, gts, output_features