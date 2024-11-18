from transformers import BertModel, BertConfig
from torch import nn
from torch.nn import functional as F
import constants
import torch

class MLM(nn.Module):
    def __init__(self, config):
        super(MLM, self).__init__()
        config = BertConfig().from_json_file(config)
        self.encoder = BertModel(config)
        self.classifier = nn.Linear(self.encoder.config.hidden_size,len(constants.tokens_list))
    def forward(self, batched_tokens):
        attention_mask = (batched_tokens != constants.tokens_to_id[constants.pad_token]).int()
        x = self.encoder(batched_tokens, attention_mask).last_hidden_state.relu()
        return self.classifier(x)


class RNAEncoder(nn.Module):
    def __init__(self, encoder, dimension = 256):
        super(RNAEncoder, self).__init__()
        self.encoder = encoder
        self.linear = nn.Linear(self.encoder.config.hidden_size,dimension)
        #self.register_buffer('eps', torch.tensor(1e-12))
    def forward(self, batched_tokens):
        attention_mask = (batched_tokens != constants.tokens_to_id[constants.pad_token]).int()
        # last_hidden_state: batch size, lenght, hidden size
        x = self.encoder(batched_tokens, attention_mask).last_hidden_state[:,0,:].relu() # consider embedding of CLS token
        x = self.linear(x)
        scaler = x.norm(dim=1, p=2, keepdim=True)
        #scaler = torch.where(scaler > self.eps, scaler, self.eps)
        x = x/scaler # force the embedding to reside on a unit sphere
        return x
