from data import MFDataset
from transformers import AutoModelForMaskedLM, AutoTokenizer
from torch.utils.data import DataLoader
import torch
from torch import nn
from tqdm import tqdm
import wandb
import numpy as np
from models import MFNetCaduceusClassifier
import torch.nn.functional as F
from sklearn.metrics import confusion_matrix, f1_score
from transformers import get_cosine_schedule_with_warmup

wandb.init(project="corn")
# Define a larger dictionary containing various settings
config = {
    "model": "caduceous",
    "learning_rate": 1e-4,
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
dataset = MFDataset(classes = config["classes"])
generator = torch.Generator().manual_seed(42)
train, test = torch.utils.data.random_split(dataset, [0.9, 0.1], generator=generator)

train_loader = DataLoader(train, batch_size=config["batch_size"], shuffle=True)
test_loader = DataLoader(test, batch_size=1, shuffle=False)

model = MFNetCaduceusClassifier(num_classes=config["classes"]).cuda()
optimizer = torch.optim.AdamW(model.parameters(), lr=config["learning_rate"], weight_decay=0.0)
schedule = get_cosine_schedule_with_warmup(optimizer, num_warmup_steps=100, num_training_steps=(len(train_loader) // config["grad_acc"]) * config["epochs"])
criterion = nn.CrossEntropyLoss()
def evaluate():
    model.eval()
    all_preds = []
    all_labels = []
    test_losses = []
    with torch.no_grad():
        for sequences, labels, sa_eff in test_loader:
            labels = labels.cuda()
            preds = model(sequences)
            loss = criterion(preds, labels)
            test_losses.append(loss.item())
            _, predicted = torch.max(preds, 1)
            all_preds.extend(predicted.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    accuracy = (np.array(all_preds) == np.array(all_labels)).mean()
    
    avg_loss = np.mean(test_losses)
    f1 = f1_score(all_labels, all_preds)
    wandb.log({"test_loss": avg_loss, "test_accuracy": accuracy, "f1_score": f1,
              "confusion_matrix": wandb.plot.confusion_matrix(probs=None, y_true=all_labels, preds=all_preds)},
              step=n_steps)    
n_steps = 0
losses = []
scaler = torch.cuda.amp.GradScaler()
for epoch in range(1, 6):
    for sequences, labels, saeff in tqdm(train_loader):
        model.train()
        labels = labels.cuda()
        with torch.cuda.amp.autocast():
            preds = model(sequences)
            one_hot_labels = F.one_hot(labels, num_classes = config["classes"]).float()
            #print("preds: ",preds.shape)
            #print("one_hot: ",one_hot_labels.shape)
            loss = criterion(preds, one_hot_labels)
        loss = loss / config["grad_acc"]
        scaler.scale(loss).backward()
        losses.append(loss.item() * config["grad_acc"])

        if n_steps % config["grad_acc"]*config["batch_size"] == 0 and n_steps > 0: 
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()
            optimizer.zero_grad()
            schedule.step()            

            wandb.log({"train_loss": np.mean(losses)}, step=n_steps)
            losses = []
        
        n_steps += config["batch_size"]
    evaluate()
