import os
import torch
import numpy as np
from torch.utils.data import Dataset
from torch.nn.utils.rnn import pad_sequence
from collections import defaultdict
import logging
from tqdm import tqdm
import sys
import constants
import re
from ushuffle import shuffle
import random
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("load sequence instances")


def tokenize(sequence,k=1):
    """
    input: RNA sequence of length L
    output: integer list of shape  L-k+1 + 2
            first element is cls token
            last token is seq token
    """
    tokens = [constants.tokens_to_id[constants.cls_token]]
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        if kmer in constants.kmer_tokens_set:
            tokens.append(constants.tokens_to_id[kmer])
        else:
            tokens.append(constants.tokens_to_id[constants.unk_token])
    tokens.append(constants.tokens_to_id[constants.sep_token])
    return torch.tensor(tokens)

def collate(examples):
    return pad_sequence(examples, batch_first=True, padding_value=constants.tokens_to_id[constants.pad_token])

def collate2(examples):
    batched_tokens = []
    batched_tags = []
    for tokens, tags in examples:
        batched_tokens.append(tokens)
        batched_tags.append(tags)
    batched_tokens = pad_sequence(batched_tokens, batch_first=True, padding_value=constants.tokens_to_id[constants.pad_token])
    batched_tags = pad_sequence(batched_tags, batch_first=True, padding_value=0)
    #for token in constants.special_tokens_list:
    #    batched_tags[batched_tokens==token] = 0
    batched_tags = batched_tags.long()
    return batched_tokens, batched_tags

def get_masked_tokens(batched_tokens, mlm_probability=0.2):
    device = batched_tokens.device
    probability_matrix = torch.full(batched_tokens.shape, mlm_probability)
    special_token_mask = torch.full(batched_tokens.shape,False)
    for token in constants.special_tokens_list:
        special_token_mask[batched_tokens==constants.tokens_to_id[token]] = True
    probability_matrix[special_token_mask] = 0
    mlm_mask = torch.bernoulli(probability_matrix).bool()
    labels = batched_tokens.detach().clone()
    labels[~mlm_mask] = -100
    # 80% of the time, we replace masked input tokens with tokenizer.mask_token ([MASK])
    masking_mask = torch.bernoulli(torch.full(batched_tokens.shape, 0.8)).bool() & mlm_mask
    batched_tokens[masking_mask] = constants.tokens_to_id[constants.mask_token]
    # 10% of the time, we replace masked input tokens with random word
    random_mask = torch.bernoulli(torch.full(batched_tokens.shape, 0.5)).bool() & mlm_mask & (~masking_mask)
    random_words = torch.randint(len(constants.special_tokens_list),len(constants.tokens_list), batched_tokens.shape, dtype=torch.long,device=device)
    batched_tokens[random_mask] = random_words[random_mask]
    # 10% of the time, left as is
    return batched_tokens, labels

def load_fasta(path):
    """
    Load fasta file into an sequence dict
    Each sequence records could span multiple lines
    """
    sequences = {}
    attrs = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if len(line)==0:
                continue
            if line.startswith(">"):
                seq_id0= line.replace(">","").strip()
                if " " in seq_id0:
                    p = seq_id0.find(" ")
                    seq_id = seq_id0[:p]
                    attr = seq_id0[p+1:].split(" ")[0]
                    attrs[seq_id] = attr
                else:
                    seq_id = seq_id0
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line.upper().replace("T","U")
    return sequences, attrs

def load_group(path):
    grouped_seq_ids = {}
    with open(path) as f:
        for line in f:
            seq_id, group_id = line.strip().split("\t")
            if group_id not in grouped_seq_ids:
                grouped_seq_ids[group_id] = []
            grouped_seq_ids[group_id].append(seq_id)
    return grouped_seq_ids

class RNASet(Dataset):
    def __init__(self, fasta, group, max_length=254, weighted=False, min_cluster_size=10):
        self.max_length =  max_length
        logger.info("load sequence ...")
        sequences, _  = load_fasta(fasta) 
        logger.info("load groupping ...")
        grouped_seq_ids = load_group(group)
        cluster_ids = sorted(list(grouped_seq_ids.keys()))
        self.sequences = []
        n_sequence = 0
        self.weights = None
        if weighted:
            self.weights = []
        logger.info("summarize instances ...")
        logger.info(f"consider clusters with size larger than {min_cluster_size} .")
        for cluster_id in cluster_ids:
            sequences_by_cluster = []
            for seq_id in grouped_seq_ids[cluster_id]:
                if seq_id not in sequences:
                    continue
                sequence = sequences[seq_id]
                sequences_by_cluster.append(sequence)
            if len(sequences_by_cluster) < min_cluster_size:
                continue
            n_sequence += len(sequences_by_cluster)
            if self.weights is not None:
                self.weights.append(np.log2(len(sequences_by_cluster)))
            self.sequences.append(sequences_by_cluster)
        n_clusters = len(self.sequences)
        if self.weights is not None:
            self.weights = np.array(self.weights)
            self.weights = self.weights/self.weights.sum()
            self.indices = np.arange(len(self.sequences))
        logger.info(f"{n_sequence} sequences are loaded into {n_clusters} clusters .")
    def __len__(self):
        return len(self.sequences)
    def __getitem__(self,idx):
        if self.weights is not None:
            idx = np.random.choice(self.indices,p = self.weights) 
        cluster_size = len(self.sequences[idx])
        cidx = np.random.randint(cluster_size)
        sequence = self.sequences[idx][cidx]
        if len(sequence) > self.max_length:
            start = np.random.randint(len(sequence)-self.max_length)
            sequence = sequence[start:start+self.max_length]
        tokens = tokenize(sequence)
        return tokens

class RNASet2(Dataset):
    def __init__(self, fasta, group, positive_fraction=0.5, 
                 shuffled_fraction=0.2, max_length=254, crop = 0.9, paired = True,
                 weighted = None, return_sequence = False, shuffle_kmer = 2,
                 n_cached_classes = 50, multiview = 0):
        """
        fasta: input RNA sequences
        group: RNA clustered to 3 level hierarchical, follows seq_id,cluster0_id,cluster1_id,cluster2_id
          cluster0_id: defined by Rfam cm model
          cluster1_id: defined by blastn with max evalue = 0.01 and min qcov = 0.6 followed by leiden clustering
          cluster2_id: defined by cd-hit-est -r 0 -c 0.8
          all sequence with same cluster0_id share a common label
        shuffled_fraction: 
          negative pairs are either sampled from sequence with different cluster0_id (with a probability of 1-shuffled_fraction), or random select a sequence and shuffle it (preserve k-mer frequency, with a probability of shuffled_fraction)
        positive_fraction:
          fraction of positive sequence pair
        max_length: max length of sequence to consider
        crop: crop a subsequence as input
        paired: sample a instance or a instance pair
        weighted: if not specified, sample each RNA family with uniform probability, else sampled with specified probability
        """
        self.max_length =  max_length
        self.paired = paired
        self.return_sequence = return_sequence
        self.n_cached_classes = n_cached_classes
        self.multiview = multiview
        logger.info("load sequence ...")
        sequences, _  = load_fasta(fasta)        
        
        logger.info("load groupping ...")
        groupped_sid = defaultdict(set)
        groupped_c2id = defaultdict(set)
        groupped_c1id = defaultdict(set)
        with open(group) as f:
            for line in f:
                sid, c0id, c1id, c2id = line.strip().split("\t")
                groupped_sid[c2id].add(sid)
                groupped_c2id[c1id].add(c2id)
                groupped_c1id[c0id].add(c1id)

        self.positive_fraction = positive_fraction
        self.crop = crop
        self.shuffled_fraction = shuffled_fraction
        self.max_length = max_length
        self.max_length_full = int(max_length/crop)
        self.n_sequence = 0
        self.n_cluster0 = 0
        self.n_cluster1 = 0
        self.n_cluster2 = 0

        self.sequences = []
        self.seq_ids = []
        c0ids = sorted(list(groupped_c1id.keys()))
        n_too_long = 0
        n_singlton = 0
        n_absent = 0
        self.shuffle_kmer = shuffle_kmer
        self.weights = None
        if weighted is not None:
            self.weights = []
            weights = {}
            with open(weighted) as f:
                for line in f:
                    c0id, count = line.strip().split("\t")
                    count = float(count)
                    weights[c0id] = count

        for c0id in c0ids:
            cluster0 = [] # a level 0 cluster
            cluster0_ids = []
            c1ids = sorted(list(groupped_c1id[c0id]))
            for c1id in c1ids:
                cluster1 = [] # a level 1 cluster
                cluster1_ids = []
                c2ids = sorted(list(groupped_c2id[c1id]))
                for c2id in c2ids:
                    cluster2 = [] # a level 2 cluster
                    cluster2_ids = []
                    for sid in groupped_sid[c2id]:
                        if sid not in sequences:
                            n_absent += 1
                            continue
                        if len(sequences[sid]) > self.max_length_full:
                            n_too_long += 1
                            continue
                        cluster2.append(sequences[sid])
                        cluster2_ids.append(sid)
                        self.n_sequence += 1
                    if len(cluster2) > 1: # not consider singleton sequence
                        cluster1.append(cluster2)
                        cluster1_ids.append(cluster2_ids)
                        self.n_cluster2 += 1
                    else:
                        n_singlton += 1
                if len(cluster1) > 0:
                    cluster0.append(cluster1)
                    cluster0_ids.append(cluster1_ids)
                    self.n_cluster1 += 1
            if len(cluster0) > 0:
                if self.weights is not None:
                    assert c0id in weights
                    self.weights.append(weights[c0id])
                self.sequences.append(cluster0)
                self.seq_ids.append(cluster0_ids)
                self.n_cluster0 += 1
        if self.weights is not None:
            self.weights = np.array(self.weights)/sum(self.weights)
        del sequences 
        del _       
        if not self.paired:
            self.update_cached_classes()
        logger.info(f"{n_absent} sequence are not present in fasta file")
        logger.info(f"{n_too_long} sequences which are too long were discarded .")
        logger.info(f"{n_singlton} singleton sequence were discarded .")
        logger.info("we have:")
        logger.info(f"{self.n_cluster0} level 0 clusters,")
        logger.info(f"{self.n_cluster1} level 1 clusters,")
        logger.info(f"{self.n_cluster2} level 2 clusters,")
        logger.info(f"{self.n_sequence} sequences in total .")
        
    def update_cached_classes(self):
        if self.weights is None:
            self.cached_classes = np.random.choice(np.arange(self.n_cluster0),size=self.n_cached_classes,replace=False)
        else:
            self.cached_classes = np.random.choice(np.arange(self.n_cluster0),size=self.n_cached_classes,p = self.weights, replace=False)
    def __len__(self):
        return self.n_sequence

    def get_instance_pair(self, idx):
        if self.weights is None:
            idx0_1 = idx%self.n_cluster0
        else:
            idx0_1 = np.random.choice(np.arange(self.n_cluster0),p=self.weights)
        idx0_2 = idx0_1
        if np.random.rand() < self.positive_fraction:
            # sample a positive pair
            idx1_1 = np.random.randint(len(self.sequences[idx0_1]))
            idx1_2 = np.random.randint(len(self.sequences[idx0_2])) # may sample from a different level 1 cluster
            idx2_1 = np.random.randint(len(self.sequences[idx0_1][idx1_1]))
            idx2_2 = np.random.randint(len(self.sequences[idx0_2][idx1_2])) 
            idx3_1 = np.random.randint(len(self.sequences[idx0_1][idx1_1][idx2_1]))
            idx3_2 = np.random.randint(len(self.sequences[idx0_2][idx1_2][idx2_2]))
            sequence_1 = self.sequences[idx0_1][idx1_1][idx2_1][idx3_1]
            sequence_2 = self.sequences[idx0_2][idx1_2][idx2_2][idx3_2]
            if idx1_1 != idx1_2:
                level = 2
            elif idx2_1 != idx2_2:
                level = 1
            else:
                level = 0
        else:
            # sample a negative pair
            idx1_1 = np.random.randint(len(self.sequences[idx0_1]))
            idx2_1 = np.random.randint(len(self.sequences[idx0_1][idx1_1]))
            idx3_1 = np.random.randint(len(self.sequences[idx0_1][idx1_1][idx2_1]))
            sequence_1 = self.sequences[idx0_1][idx1_1][idx2_1][idx3_1]
            if np.random.rand() < self.shuffled_fraction:
                # shuffle sequence as negative instance
                idx0_2 = -1
                sequence_2 = shuffle(sequence_1.encode(), self.shuffle_kmer).decode()
                level = 1
            else:
                while idx0_2 == idx0_1:
                    idx0_2 = np.random.randint(len(self.sequences))
                idx1_2 = np.random.randint(len(self.sequences[idx0_2]))
                idx2_2 = np.random.randint(len(self.sequences[idx0_2][idx1_2]))
                idx3_2 = np.random.randint(len(self.sequences[idx0_2][idx1_2][idx2_2]))
                sequence_2 = self.sequences[idx0_2][idx1_2][idx2_2][idx3_2]
                level = 0
        length_1 = int(self.crop*len(sequence_1))
        length_2 = int(self.crop*len(sequence_2))

        start_1 = np.random.randint(len(sequence_1)-length_1)
        start_2 = np.random.randint(len(sequence_2)-length_2)

        sequence_1 = sequence_1[start_1:start_1+length_1]
        sequence_2 = sequence_2[start_2:start_2+length_2]
        
        if self.return_sequence:
            seq_id_1 = self.seq_ids[idx0_1][idx1_1][idx2_1][idx3_1]
            seq_id_1 = f"{seq_id_1}:{start_1}-{start_1+length_1}"
            if idx0_2 == -1:
                seq_id_2 = seq_id_1 + ":shuffled"
            else:
                seq_id_2 = self.seq_ids[idx0_2][idx1_2][idx2_2][idx3_2]
            seq_id_2 = f"{seq_id_2}:{start_2}-{start_2+length_2}"
            return seq_id_1, seq_id_2, sequence_1, sequence_2, idx0_1, idx0_2, level

        tokens_1 = tokenize(sequence_1)
        tokens_2 = tokenize(sequence_2)
        return tokens_1, tokens_2, idx0_1, idx0_2, level
        
    
    def get_instance(self, idx):
        to_shuffle = False
        if np.random.rand() < self.positive_fraction:
            idx0 = self.cached_classes[idx%self.n_cached_classes]
        else:
            idx0 = idx%self.n_cluster0
            if np.random.rand() < self.shuffled_fraction:
                to_shuffle = True
        idx1 = np.random.randint(len(self.sequences[idx0])) 
        idx2 = np.random.randint(len(self.sequences[idx0][idx1]))
        idx3 = np.random.randint(len(self.sequences[idx0][idx1][idx2]))
        sequence = self.sequences[idx0][idx1][idx2][idx3]
        if to_shuffle:
            sequence = shuffle(sequence.encode(), self.shuffle_kmer).decode()
            label = -1
        else:
            label = idx0
        length = int(self.crop*len(sequence))
        start = np.random.randint(len(sequence)-length)
        sequence = sequence[start:start+length]
        if self.return_sequence:
            seq_id = self.seq_ids[idx0][idx1][idx2][idx3]
            if label == -1:
                seq_id = seq_id + ":shuffled"
            seq_id = f"{seq_id}:{start}-{start+length}"
            return seq_id, sequence, label
        tokens = tokenize(sequence)
        return tokens, label
    
    
    def get_multiviewed_instance(self, idx):
        if self.weights is None:
            idx0 = idx%self.n_cluster0
        else:
            idx0 = np.random.choice(np.arange(self.n_cluster0),p=self.weights)
        instances = []
        if np.random.rand() < self.positive_fraction:
            label = idx0
            for i in range(self.multiview):
                idx1 = np.random.randint(len(self.sequences[idx0])) 
                idx2 = np.random.randint(len(self.sequences[idx0][idx1]))
                idx3 = np.random.randint(len(self.sequences[idx0][idx1][idx2]))
                sequence = self.sequences[idx0][idx1][idx2][idx3]
                length = int(self.crop*len(sequence))
                start = np.random.randint(len(sequence)-length)
                sequence = sequence[start:start+length]
                tokens = tokenize(sequence)
                instances.append(tokens)
        else:
            label = -1
            idx1 = np.random.randint(len(self.sequences[idx0])) 
            idx2 = np.random.randint(len(self.sequences[idx0][idx1]))
            idx3 = np.random.randint(len(self.sequences[idx0][idx1][idx2]))
            sequence = self.sequences[idx0][idx1][idx2][idx3] 
            length = int(self.crop*len(sequence))
            start = np.random.randint(len(sequence)-length)
            sequence = sequence[start:start+length]
            for i in range(self.multiview):
                sequence = shuffle(sequence.encode(), self.shuffle_kmer).decode()
                tokens = tokenize(sequence)
                instances.append(tokens)
        return instances, label
    def __getitem__(self, idx):
        if self.paired:
            return self.get_instance_pair(idx)
        elif self.multiview > 1:
            return self.get_multiviewed_instance(idx)
        else:
            #print("cached instance:", *self.cached_classes, sep="\t")
            return self.get_instance(idx)




class NucleotideMapping(Dataset):
    def __init__(self, alignments, group, max_length=254,return_sequence = False,weighted=None):
        """
        alignments: directory contains alignments in stockholm format
        group: seq_id, rfam_id, subcluster id
        max_length: max sequence length
        """

        logger.info(f"load grouping from {group} ...")
        groupped_sid = defaultdict(set)
        groupped_c1id = defaultdict(set) 

        self.max_length = max_length
        self.return_sequence = return_sequence
        seq_ids = set()
        with open(group) as f:
            for line in f:
                sid, c0id, c1id = line.strip().split("\t")
                groupped_sid[c1id].add(sid) 
                groupped_c1id[c0id].add(c1id)
                seq_ids.add(sid)
                
        self.weights = None
        if weighted is not None:
            self.weights = []
            weights = {}
            with open(weighted) as f:
                for line in f:
                    c0id, count = line.strip().split("\t")
                    count = float(count)
                    weights[c0id] = count
                
                
        logger.info(f"load alignments from {alignments} ...")
        from Bio import AlignIO
        c0ids = set()
        for stk in os.listdir(alignments):
            c0id = stk[:stk.rfind(".")]
            if c0id not in groupped_c1id:
                continue
            c0ids.add(c0id)
        c0ids = sorted(list(c0ids))
        logger.info(f"MSA of {len(c0ids)} RNA families found .")

        self.alignments = []
        self.seq_ids = []
        self.cluster0_ids = []
        self.n_sequence = 0 
        for c0id in c0ids:
            path = os.path.join(alignments,f"{c0id}.stk")
            logger.info(f"load {path} ...")
            align = AlignIO.read(path, "stockholm")
            c0dict = {}
            for record in align:
                seq_id = record.id
                sequence = str(record.seq)
                sequence = sequence.upper().replace("T","U")
                if seq_id not in seq_ids:
                    continue
                c0dict[seq_id] = sequence
            c0_alignments = []
            c0_ids = []
            for c1id in sorted(list(groupped_c1id[c0id])):
                c1_alignments = []
                c1_ids = []
                for sid in groupped_sid[c1id]:
                    if sid not in c0dict:
                        continue
                    c1_alignments.append(c0dict[sid])
                    c1_ids.append(sid)
                if len(c1_alignments) > 1:
                    c0_alignments.append(c1_alignments)
                    c0_ids.append(c1_ids)
                    self.n_sequence += len(c1_ids)
            if len(c0_alignments) > 0:
                self.alignments.append(c0_alignments)
                self.seq_ids.append(c1_ids)
                self.cluster0_ids.append(c0id)
                if self.weights is not None:
                    assert c0id in weights
                    self.weights.append(weights[c0id])
        if self.weights is not None:
            self.weights = np.array(self.weights)/sum(self.weights)

    def __len__(self):
        return self.n_sequence

    def __getitem__(self, idx):
        if self.weights is None:
            idx0 = idx%len(self.alignments)
        else:
            idx0 = np.random.choice(np.arange(len(self.alignments)),p=self.weights)
            
        cluster0_id = self.cluster0_ids[idx0]
        idx1_1 = np.random.randint(len(self.alignments[idx0]))
        idx2_1 = np.random.randint(len(self.alignments[idx0][idx1_1]))

        idx1_2 = np.random.randint(len(self.alignments[idx0]))
        idx2_2 = np.random.randint(len(self.alignments[idx0][idx1_2]))

        if idx1_1 == idx1_2:
            while idx2_1 == idx2_2:
                idx2_2 = np.random.randint(len(self.alignments[idx0][idx1_2]))
        hardness = int(idx1_1 != idx1_2)
        sequence_1 = self.alignments[idx0][idx1_1][idx2_1]
        sequence_2 = self.alignments[idx0][idx1_2][idx2_2]
        assert len(sequence_1) == len(sequence_2)
        nucs_1 =  np.array(list(sequence_1))
        nucs_2 =  np.array(list(sequence_2))
        mask_1 = (nucs_1!=".") & (nucs_1!="-")
        mask_2 = (nucs_2!=".") & (nucs_2!="-")
        mask = mask_1 & mask_2
        similarity = (nucs_1 == nucs_2)[mask].sum()/mask.sum()
        aligned_positions_1 = torch.tensor((mask_1.cumsum() - 1)[mask])
        aligned_positions_2 = torch.tensor((mask_2.cumsum() - 1)[mask])
        # get ungapped sequence and positions of aligned nucleotides
        sequence_1 = sequence_1.replace(".","").replace("-","")
        sequence_2 = sequence_2.replace(".","").replace("-","")
        sequence_1 = sequence_1[:self.max_length]
        sequence_2 = sequence_2[:self.max_length]
        inside = (aligned_positions_1<self.max_length) & (aligned_positions_2<self.max_length)
        aligned_positions_1 = aligned_positions_1[inside]
        aligned_positions_2 = aligned_positions_2[inside]
        if self.return_sequence:
            return cluster0_id, hardness, similarity, sequence_1, sequence_2, aligned_positions_1, aligned_positions_2
        tokens_1 = tokenize(sequence_1)
        tokens_2 = tokenize(sequence_2)
        #return cluster0_id, hardness, similarity, tokens_1, tokens_2, aligned_positions_1 + 1, aligned_positions_2 + 1 # consider the prepended <CLS> token
        aligned_positions_1, aligned_positions_2 = aligned_positions_1 + 1, aligned_positions_2 + 1
        #p_sep_1, p_sep_2 = tokens_1.shape[0] - 1, tokens_2.shape[0] - 1
        ##unaligned_mask_1 = torch.full(tokens_1.shape[0],True)
        #unaligned_mask_1[aligned_positions_1] = False
        #unaligned_mask_2 = torch.full(tokens_2.shape[0])
        #unaligned_mask_2[aligned_positions_2] = False
        correspondance = torch.full((self.max_length+2,self.max_length+2),False)
        correspondance[aligned_positions_1, aligned_positions_2] = True
        #correspondance[unaligned_mask_1,p_sep_2] = True
        #correspondance[p_sep_1,unaligned_mask_2] = True
        return tokens_1, tokens_2, correspondance
