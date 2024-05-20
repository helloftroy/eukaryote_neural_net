from data import MFDataset
from transformers import AutoModelForMaskedLM, AutoTokenizer, AutoModelForSequenceClassification
from torch.utils.data import DataLoader
import torch
from torch import nn
from tqdm import tqdm
import wandb
import numpy as np
from models import MFNetRegressor, MFNetClassifier, MFNetClassifierWithStart
import torch.nn.functional as F
from scipy.ndimage import gaussian_filter1d
from transformers import get_cosine_schedule_with_warmup
from sklearn.metrics import confusion_matrix, f1_score, precision_recall_curve

wandb.init(project="corn")
# Define a larger dictionary containing various settings
config = {
    "model": "100m NT",
    "learning_rate": 1e-5,
    "epochs": 5,
    "batch_size": 4,
    "classes": 5,
    "grad_acc": 16,
    "class_weighting": "na",
    "other": "gene with their abundance"
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
model = MFNetClassifierWithStart(num_classes=config["classes"]).cuda()
#print(model)

"""
# just takes the x bottom layers of the model and the classifier, so its smaller for faster training
for param in model.parameters():
    param.requires_grad = False

for param in model.model.classifier.parameters():
    param.requires_grad = True

for param in model.model.esm.encoder.layer[-7:].parameters():
    param.requires_grad = True

model_parameters = filter(lambda p: p.requires_grad, model.parameters())
params = sum([np.prod(p.size()) for p in model_parameters])
#print("Parameters: ", params)
"""
optimizer = torch.optim.Adam(model.parameters(), lr=config["learning_rate"])
schedule = get_cosine_schedule_with_warmup(optimizer, num_warmup_steps=100, num_training_steps=(len(train_loader) // config["grad_acc"]) * config["epochs"])
# use when weighting
#class_weights = torch.tensor([1.0, 9.0]).cuda()  # Assigning a weight ratio of 9:1 based on the class distribution
#criterion = nn.CrossEntropyLoss(weight=class_weights)
criterion = nn.CrossEntropyLoss()

def evaluate():
    model.eval()

    test_losses = []
    all_preds = []
    all_labels = []
    correct = 0
    total = 0
    #all_sa = []
    with torch.no_grad():
        for sequences, labels, sa_eff in test_loader:
            labels = labels.cuda()
            preds = model(sequences)
            loss = criterion(preds, labels)  # Calculate the loss using the specified criterion
            test_losses.append(loss.item())
            _, predicted = torch.max(preds, 1)
            all_preds.extend(predicted.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())
            #print(sa_eff.item(), predicted.item(), labels.item())
        #    all_sa.append([sa_eff,labels,predicted])

            
    accuracy = (np.array(all_preds) == np.array(all_labels)).mean()#correct / total
    avg_loss = np.mean(test_losses)
    f1 = f1_score(all_labels, all_preds, average='macro')
    print("Test Loss:", avg_loss,"\nTest Accuracy:", accuracy, "\nF1 Score:", f1, "\nConfusion Matrix:" )
    cm = confusion_matrix(all_labels, all_preds)
    for i in range(len(cm)):
        print(f"Actual {i}\t" + "\t".join(str(x) for x in cm[i]))
    #wandb.log({"my_lineplot_id" : wandb.plot.line(table, "recall_micro", 
     #      "precision_micro", stroke=None, title="Average Precision")})

    wandb.log({"test_loss": avg_loss, "test_accuracy": accuracy, "f1_score": f1,
            "confusion_matrix": wandb.plot.confusion_matrix(probs=None, y_true=all_labels, preds=all_preds)},
            step=n_steps)

n_steps = 0
losses = []
scaler = torch.cuda.amp.GradScaler()

for epoch in range(1, 6):
    evaluate()
    for sequences, labels, saeff in tqdm(train_loader):
        model.train()
        labels = labels.cuda()
        
        with torch.cuda.amp.autocast():
            
            preds = model(sequences)
            #apply smoothing to the targets
            one_hot_labels = F.one_hot(labels, num_classes=config["classes"]).float()
            #smoothed_labels = torch.from_numpy(gaussian_filter1d(one_hot_labels.cpu().numpy(), 0.5)).cuda()
            #filterdData = gaussian_filter1d( labels, 0.47 )
            loss = criterion(preds, one_hot_labels)
            #loss = criterion(preds, smoothed_labels)  # Calculate the loss using the specified criterion
        loss = loss / config["grad_acc"]  # Update the loss calculation for gradient accumulation

        scaler.scale(loss).backward()
 
        losses.append(loss.item())

        if n_steps % (config["grad_acc"]*config["batch_size"]) == 0 and n_steps > 0: 
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()
            optimizer.zero_grad()
            schedule.step()
            wandb.log({"train_loss": np.sum(losses),
                    "lr": optimizer.param_groups[0]["lr"]}, step=n_steps)
            losses = []
        
        n_steps += config["batch_size"]
    evaluate()
