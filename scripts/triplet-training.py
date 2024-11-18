#!/usr/bin/env python
import torch
import os
import numpy as np
from torch.utils.data import DataLoader, SequentialSampler, RandomSampler
from torch.nn.utils.rnn import pad_sequence
from model import MLM, RNAEncoder
import sys
from torch.nn import functional as F
from torch.optim import AdamW
from dataset import RNASet2
import constants
from transformers import BertModel, BertConfig, get_linear_schedule_with_warmup
from sklearn.metrics import precision_recall_curve,auc,roc_auc_score
constants.init(1)
import argparse

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("training with triplet loss")

eps = 1e-8

def cosine_distance(embedding_1,embedding_2):
    return 1 - (embedding_1 * embedding_2).sum(axis=-1,keepdim=True)

def euclidean_distance(embedding_1,embedding_2):
    return (embedding_1 - embedding_2).pow(2).sum(1)


def pairwise_cosine_distance(embedding):
    dp = torch.mm(embedding, embedding.t()) # dot product / pairwise cosine distance
    dg = torch.diag(dp) # diagonal elements, equals to 1 
    distance = (dg.unsqueeze(0) - 2.0*dp + dg.unsqueeze(1))/2 # convert similarity to distance
    distance = F.relu(distance) # make sure distance >= 0
    return distance

def get_anchor_positive_triplet_mask(labels):
    """Return a 2D mask where mask[a, p] is True iff a and p are distinct and have same label.
    Args:
        labels: tf.int32 `Tensor` with shape [batch_size]
    Returns:
        mask: tf.bool `Tensor` with shape [batch_size, batch_size]
    """
    # Check that i and j are distinct
    indices_equal = torch.eye(labels.size(0), device=labels.device).bool()
    indices_not_equal = ~indices_equal

    # Check if labels[i] == labels[j]
    # Uses broadcasting where the 1st argument has shape (1, batch_size) and the 2nd (batch_size, 1)
    labels_equal = (labels.unsqueeze(0) == labels.unsqueeze(1)) &  (labels.unsqueeze(0) != -1)
    return labels_equal & indices_not_equal 


def get_anchor_negative_triplet_mask(labels):
    """Return a 2D mask where mask[a, n] is True iff a and n have distinct labels.
    Args:
        labels: (batch size,)
    Returns:
        mask: bool (batch size, batch size)
    """
    # Check if labels[i] != labels[k]
    # Uses broadcasting where the 1st argument has shape (1, batch_size) and the 2nd (batch_size, 1)

    return (~(labels.unsqueeze(0) == labels.unsqueeze(1))) | (labels.unsqueeze(1) == -1)


def triplet_loss(embeddings, labels,margin = 0.32):
    # calculate within batch pairwise distance
    distance = pairwise_cosine_distance(embeddings)
    mask_anchor_positive = get_anchor_positive_triplet_mask(labels).float()
    mask_anchor_negative = get_anchor_negative_triplet_mask(labels).float()
    instance_mask = (mask_anchor_positive.sum(axis=1,keepdim=True) > 0) & (mask_anchor_negative.sum(axis=1,keepdim=True) > 0)
    anchor_positive_dist = mask_anchor_positive * distance # (batch size, batch size)
    hardest_positive_dist, _ = anchor_positive_dist.max(1, keepdim=True) 
    max_distance, _ = distance.max(1, keepdim=True) # for each anchor, max distance with in this mini batch
    distance2 = distance + max_distance*(1.0 - mask_anchor_negative) # set unmasked value (not negative pair) as the max distance
    hardest_negative_dist, _ = distance2.min(1, keepdim=True) # get the hardest negative instance of each anchor
    loss = hardest_positive_dist - hardest_negative_dist + margin
    loss = loss[instance_mask]
    return F.relu(loss).mean()



def evaluate(embedding_1,embedding_2,label_1,label_2):
    distances = cosine_distance(embedding_1, embedding_2)
    labels = (label_1 == label_2).int()
    scores = -distances.detach().cpu().numpy()
    labels = labels.detach().cpu().numpy()
    AUROC = roc_auc_score(labels, scores)
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    AUPRC = auc(recall, precision)
    return AUROC, AUPRC 


def collate_pairs(examples):
    batched_tokens_1, batched_tokens_2 = [], []
    labels_1, labels_2, levels = [], [], []
    for tokens_1, tokens_2, label_1, label_2, level in examples:
        batched_tokens_1.append(tokens_1)
        batched_tokens_2.append(tokens_2)
        labels_1.append(label_1)
        labels_2.append(label_2)
        levels.append(level)
    batched_tokens_1 = pad_sequence(batched_tokens_1, batch_first=True, padding_value=constants.tokens_to_id[constants.pad_token])
    batched_tokens_2 = pad_sequence(batched_tokens_2, batch_first=True, padding_value=constants.tokens_to_id[constants.pad_token])
    labels_1 = torch.tensor(labels_1).long()
    labels_2 = torch.tensor(labels_2).long()
    levels = torch.tensor(levels).long()
    return batched_tokens_1, batched_tokens_2, labels_1, labels_2, levels


def collate_fn(examples):
    batched_tokens, labels  = [], []
    for tokens, label  in examples:
        batched_tokens.append(tokens)
        labels.append(label)
    batched_tokens = pad_sequence(batched_tokens, batch_first=True, padding_value=constants.tokens_to_id[constants.pad_token])
    labels = torch.tensor(labels).long()
    return batched_tokens, labels


def main():
    parser = argparse.ArgumentParser(description='train RNA encoder with contrastive loss')
    parser.add_argument('--encoder-config','-ec',required=True, help="model configuration")
    parser.add_argument('--dimension','-dim',type=int,default=256, help="dimenion of embedding vector")
    parser.add_argument('--train-sequence','-ts',type=str,required=True,help="sequence for training")
    parser.add_argument('--train-group','-tg',type=str,required=True,help="groupping of training instance")
    parser.add_argument('--val-sequence','-vs',type=str,help="sequence for validation")
    parser.add_argument('--val-group','-vg',type=str,required=True,help="groupping of validation instance")
    parser.add_argument('--positive-fraction','-pf',type=float,default=0.5,help="positive pair fraction in training")
    parser.add_argument('--shuffled-fraction','-sf',type=float,default=0.5,help="shuffling this fraction in negative instances")
    parser.add_argument('--batch-size', '-bs', type=int,default=128, help="batch size for scanning")
    parser.add_argument('--weight-decay', '-wd', type=float, default=0.01,help="weight decay to use")
    parser.add_argument('--warmup-steps',default=0,type=int,help="number of steps for learning rate warm up")
    parser.add_argument('--device', '-d', default="cuda:0", choices=["cuda:0","cuda:1","cpu"],help="Device to run the model")
    parser.add_argument('--models', '-m', required=True, type=str, help="directory to save model check points")
    parser.add_argument('--pretrained-model', '-pm', help="whether start from pretrained model .")
    parser.add_argument('--max-grad-norm',default=1,type=float,help="gradient clipping")
    parser.add_argument('--metrics', help="where to save training metrics .")
    parser.add_argument('--epoches', '-e', type=int, default=256, help="Number of epoches to train")
    parser.add_argument('--learning-rate', '-lr', type=float, default=5e-5, help="Learning rate ro use")
    parser.add_argument('--margin', '-mg', type=float, default=0.32, help="margin for contrastive loss")
    args = parser.parse_args()


    if args.pretrained_model is not None:
        logger.info("Load pretrained weights ...")
        pretranied_model = MLM(args.encoder_config)
        state_dict = torch.load(args.pretrained_model)
        pretranied_model.load_state_dict(state_dict)
        logger.info("Intialize RNA encoder with pretrained weights ...")
        encoder = pretranied_model.encoder
        del pretranied_model
    else:
        logger.info("Intialize RNA encoder with random weights ...")
        config = BertConfig().from_json_file(args.encoder_config)
        encoder = BertModel(config)
    model = RNAEncoder(encoder,dimension=args.dimension).to(args.device)

    logger.info(f"Model checkpoints will be saved to {args.models} .")

    if not os.path.exists(args.models):
        logger.warn(f"the directory {args.models} does not exists, create it .")
        os.makedirs(args.models)

    logger.info("load training instances ...")
    train_set = RNASet2(args.train_sequence, args.train_group, 
                           positive_fraction = args.positive_fraction, 
                           shuffled_fraction = args.shuffled_fraction, crop = 0.9, paired = False)
    train_sampler = RandomSampler(train_set)

    logger.info("load validation instances ...")
    val_set = RNASet2(args.val_sequence, args.val_group,
                           positive_fraction = 0.2, shuffle_kmer=2,
                           shuffled_fraction = 0.0, crop = 0.9, paired = True)
    val_sampler = RandomSampler(val_set)

    train_loader = DataLoader(train_set, sampler=train_sampler, batch_size=args.batch_size, collate_fn=collate_fn)
    val_loader = DataLoader(val_set, sampler=val_sampler, batch_size=args.batch_size, collate_fn=collate_pairs)


    # not apply weight decay to bias and layerNorm weight
    no_decay = ["bias", "LayerNorm.weight"]
    optimizer_grouped_parameters = [{
        "params": [p for n, p in model.named_parameters() if not any(nd in n for nd in no_decay)],
        "weight_decay": args.weight_decay,
        },{
        "params": [p for n, p in model.named_parameters() if any(nd in n for nd in no_decay)],
        "weight_decay": 0},]
    
    optimizer = AdamW(optimizer_grouped_parameters,lr=args.learning_rate, eps=1e-8, betas=(0.9,0.999))
    total_steps = len(train_loader)*args.epoches
    scheduler = get_linear_schedule_with_warmup(optimizer, 
                                                num_warmup_steps=args.warmup_steps, 
                                                num_training_steps=total_steps)

    
    logger.info(f"training metrics will be saved to {args.metrics}")
    fm = open(args.metrics,"w")

    logger.info("start training ...")
    for e in range(args.epoches):
        i = 0
        train_losses = []
        for batched_tokens, labels in train_loader:
            batched_tokens, labels = batched_tokens.to(args.device), labels.to(args.device)
            optimizer.zero_grad()
            embeddings = model(batched_tokens) # (batch size, hidden dimension)
            loss = triplet_loss(embeddings, labels, margin = args.margin)
            loss.backward()
            train_losses.append(loss.item())
            torch.nn.utils.clip_grad_norm_(model.parameters(), args.max_grad_norm) # apply gradient clipping
            optimizer.step() # update weights
            scheduler.step() # update learning rate
            train_set.update_cached_classes()
            if i%256 == 0:
                # summarize training metrics 
                print("train", e, i, round(np.mean(train_losses),5), file = fm, sep="\t")
                train_losses = []
                # get validation performance 
                j = 0
                val_AUROCs, val_AUPRCs = [], []
                val_hard_AUROCs, val_hard_AUPRCs = [], []
                model = model.eval()
                for batched_tokens_1, batched_tokens_2, labels_1, labels_2, levels in val_loader:
                    batched_tokens_1, batched_tokens_2, labels_1, labels_2, levels = batched_tokens_1.to(args.device), batched_tokens_2.to(args.device), labels_1.to(args.device), labels_2.to(args.device), levels.to(args.device)
                    embeddings_1 = model(batched_tokens_1) # (batch size, hidden dimension)
                    embeddings_2 = model(batched_tokens_2)
                    AUROC, AUPRC = evaluate(embeddings_1, embeddings_2, labels_1, labels_2)
                    hard_mask = (levels == 2) | (labels_1 != labels_2)
                    hard_AUROC, hard_AUPRC = evaluate(embeddings_1[hard_mask,:], embeddings_2[hard_mask,:], labels_1[hard_mask], labels_2[hard_mask])
                    val_AUROCs.append(AUROC)
                    val_AUPRCs.append(AUPRC)
                    val_hard_AUROCs.append(hard_AUROC)
                    val_hard_AUPRCs.append(hard_AUPRC)
                    if j >= 100:
                        break
                    j += 1
                # summarize validation metrics 
                print("val", e, i, round(np.mean(val_AUROCs),3), round(np.mean(val_AUPRCs),3),
                                   round(np.mean(val_hard_AUROCs),3), round(np.mean(val_hard_AUPRCs),3),
                        file= fm, sep = "\t")
                fm.flush()
                model = model.train()
            i += 1
        torch.save(model.state_dict(),f"{args.models}/{e}.{i}.pt")
    logger.info("all done .")

if __name__ == "__main__":
    main()
