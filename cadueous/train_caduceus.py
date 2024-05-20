from data import MFDataset
from transformers import AutoModelForMaskedLM, AutoTokenizer
from torch.utils.data import DataLoader
import torch
from torch import nn
from tqdm import tqdm
import wandb
import numpy as np
from models import *

wandb.init(project="corn")
# Define a larger dictionary containing various settings
config = {
    "model": "caduceous",
    "learning_rate": 1e-5,
    "epochs": 5,
    "batch_size": 4,
    "classes": 2,
    "grad_acc": 16,
    "class_weighting": "na",
    "other": "genes, abundance"
    # Add more settings as needed
}

# Save the configuration
wandb.config.update(config)

# initialize datasets and dataloader
dataset = MFDataset()
generator = torch.Generator().manual_seed(42)
train, test = torch.utils.data.random_split(dataset, [0.9, 0.1], generator=generator)

train_loader = DataLoader(train, batch_size=config["batch_size"], shuffle=True)
test_loader = DataLoader(test, batch_size=1, shuffle=False)

model = MFNetCaduceus().cuda()
optimizer = torch.optim.AdamW(model.parameters(), lr=config["learning_rate"], weight_decay=0.0)

def evaluate():
    model.eval()

    test_losses = []
    for sequences, labels in test_loader:
        labels = labels.cuda()
        with torch.no_grad():
            preds = model(sequences)
        loss = (labels - preds)**2
        loss = loss.mean()    
            
        test_losses.append(loss.item())
    wandb.log({"test_loss": np.mean(test_losses)})
    return np.mean(test_losses)

n_steps = 0
losses = []

for epoch in range(1, 6):
    for sequences, labels, saeff in tqdm(train_loader):
        model.train()
        labels = labels.cuda()

        preds = model(sequences)

        loss = (labels - preds)**2
        loss = loss.mean() / config["grad_acc"]
        loss.backward()
        losses.append(loss.item())

        if n_steps % config["grad_acc"] == 0 and n_steps > 0: 
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            optimizer.zero_grad()

            wandb.log({"train_loss": np.sum(losses)}, step=n_steps)
            losses = []
        
        n_steps += config["batch_size"]

    evaluate()
