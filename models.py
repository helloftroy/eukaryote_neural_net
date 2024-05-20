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


class MFNetClassifierWithStart(nn.Module):

    def __init__(self, num_classes):
        super().__init__()

        self.tokenizer = AutoTokenizer.from_pretrained(f'InstaDeepAI/{model_name}',  trust_remote_code=True)
        self.model = AutoModelForMaskedLM.from_pretrained(
            'InstaDeepAI/nucleotide-transformer-v2-100m-multi-species',  trust_remote_code=True)
        self.ff = FeedForward(dim=512, hidden_dim=1024*2, dropout=0.1)        #+2 if also include the length
        self.head = nn.Linear(512, num_classes) # +2 if also include the length # 1024 if 500m model

    def forward(self, sequences):#, lengths):

        tokens = self.tokenizer(sequences, padding="max_length", max_length=1024, truncation=True)["input_ids"]
        tokens = torch.tensor(tokens).cuda()
        attention_mask = (tokens != self.tokenizer.pad_token_id)
        outs = self.model(tokens,
            attention_mask=attention_mask,
            encoder_attention_mask=attention_mask,
            output_hidden_states=True
        )
        x = outs['hidden_states'][-1]
        #x[:,20:-20] = 0 # set between 20 and 20 == 0 # IF YOU SET MAKE SURE SEQUENCE IS NOT TOO SHORT
        #front_mean = x[:,:5].mean(dim=1)
        #end_mean = x[:,-5:].mean(dim=1)
        #x = torch.cat([front_mean, end_mean], dim=-1)
        x = x.mean(dim=1)
        #if using 2 inputs 
        #x = torch.cat([x, starts.unsqueeze(1).float(), length.unsqueeze(1).float()], dim=-1)
    #    x = torch.cat([x, starts.unsqueeze(1).float()], dim=-1)
        x = self.ff(x)
        x = self.head(x)
        #print(outs.logits.shape)
        return x

class MFNetClassifier(nn.Module):

    def __init__(self, num_classes):
        super().__init__()

        self.tokenizer = AutoTokenizer.from_pretrained(f'InstaDeepAI/{model_name}',  trust_remote_code=True)
        self.model = AutoModelForSequenceClassification.from_pretrained(
            'InstaDeepAI/nucleotide-transformer-v2-100m-multi-species', num_labels=num_classes,  trust_remote_code=True)

    def forward(self, sequences):

        tokens = self.tokenizer(sequences, padding="max_length", max_length=1024, truncation=True)["input_ids"]
        tokens = torch.tensor(tokens).cuda()
        attention_mask = (tokens != self.tokenizer.pad_token_id)
        outs = self.model(tokens,
            attention_mask=attention_mask,
            #encoder_attention_mask=attention_mask,
            #output_hidden_states=True
        )
        #x = outs['hidden_states'][-1]
        #x = x.mean(dim=1)
        #x = self.head(x)
        #print(outs.logits.shape)
        return outs.logits


class MFNetRegressor(nn.Module):
    
    def __init__(self):
        super().__init__()    
        
        self.tokenizer = AutoTokenizer.from_pretrained(f'InstaDeepAI/{model_name}', trust_remote_code=True)
        self.embedder = AutoModelForMaskedLM.from_pretrained(f'InstaDeepAI/{model_name}', trust_remote_code=True)
        self.head = nn.Linear(1500, 1)

    def forward(self, sequences):
 
        tokens = self.tokenizer(sequences, padding="max_length", max_length=1024, truncation=True)["input_ids"]
        tokens = torch.tensor(tokens).cuda()

        attention_mask = (tokens != self.tokenizer.pad_token_id)

        outs = self.embedder(
            tokens,
            attention_mask=attention_mask,
            encoder_attention_mask=attention_mask,
            output_hidden_states=True
        )
        x = outs['hidden_states'][-1]
        x = x.mean(dim=1)
        x = self.head(x)
        return x



class MFNetCaduceus(nn.Module):

    def __init__(self):
        super().__init__()

        self.tokenizer = AutoTokenizer.from_pretrained('kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16', trust_remote_code=True)
        self.caduceus = AutoModel.from_pretrained('kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16', trust_remote_code=True)
        self.ff = FeedForward(dim=512, hidden_dim=1024*2, dropout=0.2)
        self.head = nn.Linear(512, 1)

    def forward(self, sequences):

        tokens = torch.tensor(self.tokenizer.encode(sequences[0]), dtype=torch.int).cuda().unsqueeze(0).long()

        x = self.caduceus(tokens).mean(dim=1)
        x = self.ff(x)
        x = self.head(x)
        return x

