from data import MFDataset
from transformers import AutoModelForMaskedLM, AutoTokenizer, AutoModelForSequenceClassification, AutoModel, AutoTokenizer
from torch.utils.data import DataLoader
import torch
from torch import nn
from tqdm import tqdm
import wandb
import numpy as np
#from peft import get_peft_config, get_peft_model, LoraConfig, TaskType

# initialize model
#model_name = 'agro-nucleotide-transformer-1b'
model_name = 'nucleotide-transformer-v2-100m-multi-species'
class FeedForward(nn.Module):
    def __init__(self, dim, hidden_dim, dropout = 0.):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(dim),
            nn.Linear(dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, dim),
            nn.Dropout(dropout)
        )

    def forward(self, x):
        return self.net(x)

class MFNetCaduceus(nn.Module):

    def __init__(self):
        super().__init__()

        self.tokenizer = AutoTokenizer.from_pretrained('kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16', trust_remote_code=True)
        self.caduceus = AutoModel.from_pretrained('kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16', trust_remote_code=True)
        self.ff = FeedForward(dim=512, hidden_dim=1024*2, dropout=0.2)
        self.head = nn.Linear(512, 1)

    def forward(self, sequences):

        tokens = torch.tensor(self.tokenizer.encode(sequences[0]), dtype=torch.int).cuda().unsqueeze(0).long()
        x = self.caduceus(tokens).last_hidden_state.mean(dim=1)
        x = self.ff(x)
        x = self.head(x)
        return x

class MFNetCaduceusClassifier(nn.Module):
    def __init__(self, num_classes):
        super().__init__()
        self.tokenizer = AutoTokenizer.from_pretrained('kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16', trust_remote_code=True)
        self.caduceus = AutoModel.from_pretrained('kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16', num_labels=num_classes, trust_remote_code=True)
        self.ff = FeedForward(dim=512, hidden_dim=1024*2, dropout=0.1)
        self.head = nn.Linear(512, num_classes)
    def forward(self, sequences):
        tokens = self.tokenizer(sequences, padding=True, return_tensors="pt")["input_ids"].cuda().long()
        #tokens = torch.tensor(tokens).cuda()
        outs = self.caduceus(tokens).last_hidden_state
        outs = self.ff(outs).mean(dim=1)
        outs = self.head(outs)
        return outs
  
