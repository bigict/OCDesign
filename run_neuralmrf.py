import argparse
import time
import torch
import torch.nn as nn
import torch.nn.functional as F
import json
import torch.utils.data as Data
import random
from protein_mpnn_utils import loss_nll, loss_smoothed, gather_edges, gather_nodes, gather_nodes_t, cat_neighbors_nodes, _scores, _S_to_seq, tied_featurize, parse_PDB
from protein_mpnn_utils import StructureDataset, StructureLoader, StructureDatasetPDB, ProteinMPNN, processData, reduce_mean
import numpy as np
from torch.utils.data.distributed import DistributedSampler
import os
import copy

def parse_to_json(pdb_path):
    aminodict = {"GLY": "G", "ALA": "A", "VAL": "V", "LEU": "L", "ILE": "I", "PRO": "P", "PHE": "F", "TYR": "Y", "TRP": "W", "SER": "S",
                 "THR": "T", "CYS": "C", "MET": "M", "ASN": "N", "GLN": "Q", "ASP": "D", "GLU": "E", "LYS": "K", "ARG": "R", "HIS": "H"}
    file = open(pdb_path)
    line = file.readline()
    num_of_chains = 0
    seq = ""
    resSeq = 0
    N_chain = []
    CA_chain = []
    C_chain = []
    O_chain = []
    switch = False
    
    
    name = os.path.basename(pdb_path).split('.')[0]
    while line:
        if line[0:6] == "HEADER":
            max = len(line) - 1
            while(line[max - 1: max] == " "):
                max -= 1
            min = max
            while(line[min - 1: min] != " "):
                min -= 1
            name = line[min: max] + "_" + args.chain
            num_of_chains += 1
        
        if line[0:4] == "ATOM":
            if int(line[23:26]) != resSeq and line[21:22] == args.chain:
                seq += aminodict[line[17:20]]
                resSeq = int(line[23:26])
            if line[13:16].strip() == "N" and line[21:22] == args.chain:
                N_chain.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
            if line[13:16].strip() == "CA" and line[21:22] == args.chain:
                CA_chain.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
            if line[13:16].strip() == "C" and line[21:22] == args.chain:
                C_chain.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
            if line[13:16].strip() == "O" and line[21:22] == args.chain:
                O_chain.append([float(line[31:38]), float(line[39:46]), float(line[47:54])])
            
        line = file.readline()
    
    file.close()
    
    '''
    if args.fix_native_pos is not None and args.fix_native_val is not None:
        for i, pos in enumerate(args.fix_native_pos):
            test = list(seq)
            test[int(pos)] = args.fix_native_val[int(i)]
            seq = ''.join(test)
    '''
    
    f = open("test.jsonl", 'w')
    f.write('{"seq_chain_' + args.chain + '": "' + seq + '", "coords_chain_' + args.chain + '": {"N_chain_' + args.chain + '": [')

    for n in N_chain:
        if switch:
            f.write(', ')
        f.write('[' + str(n[0]) + ', ' + str(n[1]) + ', ' + str(n[2]) + ']')
        switch = True
    switch = False
    f.write('], "CA_chain_' + args.chain + '": [')

    for ca in CA_chain:
        if switch:
            f.write(', ')
        f.write('[' + str(ca[0]) + ', ' + str(ca[1]) + ', ' + str(ca[2]) + ']')
        switch = True
    switch = False
    f.write('], "C_chain_' + args.chain + '": [')

    for c in C_chain:
        if switch:
            f.write(', ')
        f.write('[' + str(c[0]) + ', ' + str(c[1]) + ', ' + str(c[2]) + ']')
        switch = True
    switch = False
    f.write('], "O_chain_' + args.chain + '": [')   
        
    for o in O_chain:
        if switch:
            f.write(', ')
        f.write('[' + str(o[0]) + ', ' + str(o[1]) + ', ' + str(o[2]) + ']')
        switch = True
    switch = False
        
    f.write(']}, "name": "' + name + '", "num_of_chains": ' + str(num_of_chains) + ', "seq": "' + seq + '"}')
    f.close()
    
    

def test(args):
    time_1 = time.time()
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '5678'
    torch.distributed.init_process_group('nccl', world_size=args.num, rank=args.num-1)
    
    if args.seed != 0:
        seed=args.seed
    else:
        seed=int(np.random.randint(0, high=999, size=1, dtype=int)[0])
    
    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)
    
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
    bias_AAs_np = np.zeros(len(alphabet))
    omit_AAs_np = np.array([AA in "X" for AA in alphabet]).astype(np.float32)

    if args.CASP14_test_set:
        test_set_path = "/home/wangxinglong/MPNN_dataset/CASP14.jsonl"
    else:
        test_set_path = "/home/wangxinglong/MPNN_dataset/CAMEO2022_601.jsonl"
    # test_dataset = processData(test_set_path)
    # test_dataset = processData(args.json_path)
    
    parse_to_json(args.pdb_path)
    test_dataset = processData("test.jsonl")
    
    checkpoint_path = args.checkpoint_path 
    device = torch.device(f"cuda:{args.device}")
    model = torch.load(checkpoint_path, map_location=device)
    model = model.module
    
    # config the fix residue for fix position
    fix_dict, fix_positions = parse_fix_dict(args.fix_file)
    if len(fix_dict) == 0: fix = False
    else: fix = True

    identity_cont_relax_list, identity_design_list, dataSet, target_value_list = [], [], [], []
    D = []

    for ix, protein in enumerate(test_dataset):
        entry = {}
        entry["name"] = protein["name"]
        batch_clones = [copy.deepcopy(protein) for i in range(1)]
        X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list, visible_list_list, masked_list_list, masked_chain_length_list_list, chain_M_pos, residue_idx, dihedral_mask, tied_pos_list_of_lists_list, pssm_coef, pssm_bias, pssm_log_odds_all, bias_by_res_all, tied_beta = tied_featurize(batch_clones, device, None)
        if args.ContGreedy:
            # S_sample[0] = model.sample(X, chain_encoding_all, residue_idx, mask, S, args.fix, args.ContGreedy, args.fix_positions)
            S_sample[0] = model.sample(X, chain_encoding_all, residue_idx, mask, S, fix, args.ContGreedy, fix_dict = fix_dict, fixed_positions = fix_positions)
            mask_ = torch.ones(1, S.size(1)).squeeze().numpy()
            native_seq = _S_to_seq(S[0].cpu().numpy(), mask_)
            cont_relax_seq = _S_to_seq(cont_relax_seqs[0].cpu().numpy(), mask_)
            designed_seq = _S_to_seq(S_sample[0].cpu().numpy(), mask_)
            entry["Native_sequence"] = native_seq
            entry["Cont_relax_sequence"] = cont_relax_seq
            entry["Designed_sequence"] = designed_seq
            """
            print("Native_sequence:",native_seq)
            print("Cont_relax_sequence:",cont_relax_seq)
            print("Designed_sequence:",designed_seq)
            """
            identity_cont_relax = cont_relax_seqs[0].eq(S[0]).cpu().numpy().sum() / S.size(1)
            identity_design = S_sample[0].eq(S[0]).cpu().numpy().sum() / S.size(1)
            identity_cont_relax_list.append(identity_cont_relax)
            identity_design_list.append(identity_design)
            entry["Identity_cont_relax"] = identity_cont_relax
            entry["Identity_design"] = identity_design
            entry["Target_value_ContRelax"] = target_values_ContRelax[0]
            entry["Target_value_ContGreedy"] = target_values_ContGreedy[0]
            """
            print("ix:",ix)
            print("identity_cont_relax:",identity_cont_relax)
            print("identity_design:",identity_design)
            print("avg_identity_cont_relax:",torch.mean(torch.tensor(identity_cont_relax_list)).item())
            print("avg_identity_design:",torch.mean(torch.tensor(identity_design_list)).item())
            """
            dataSet.append(entry)

        else:
            S_sample = model.sample(X, chain_encoding_all, residue_idx, mask, S, fix, args.ContGreedy, fix_dict = fix_dict, fixed_positions = fix_positions)
            #print("S_sample:",S_sample)
            time_2 = time.time()
            mask_ = torch.ones(1, S.size(1)).squeeze().numpy()
            native_seq = _S_to_seq(S[0].cpu().numpy(), mask_)
            designed_seq = _S_to_seq(S_sample[0].cpu().numpy(), mask_)
            entry["Native_sequence"] = native_seq
            entry["Designed_sequence"] = designed_seq
            """
            print("Native_sequence:",native_seq)
            print("Designed_sequence:",designed_seq)
            """
            identity = S_sample[0].eq(S[0]).cpu().numpy().sum() / S.size(1)
            identity_design_list.append(identity)
            entry["Identity"] = identity
            """
            print("ix:",ix)
            print("identity:",identity)
            print("avg_identity:",torch.mean(torch.tensor(identity_design_list)).item())
            """
            dataSet.append(entry)
    
    # output to the stdout/write to the file
    for entry in dataSet:
        output_info = ">" + str(entry["name"]) + " Identity:" + str(entry["Identity"]) + "\n" + str(entry["Designed_sequence"])
        if args.output_file is not None:
            with open(args.output_file, "w") as f:f.write(output_info)
        print(f'=' * 80)
        print(output_info)
        print(f'=' * 80)

def parse_fix_dict(fix_file):
    with open(fix_file) as f:lines = [line.strip() for line in f.readlines()]
    fix_dict = {int(line.split(' ')[0]): line.split(' ')[1] for line in lines}
    fix_positions = [i - 1 for i in list(fix_dict.keys())]
    return fix_dict, fix_positions
            
if __name__ == "__main__":
    # using argparser parse arguments
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # base config 
    argparser.add_argument("--num", type=int, default=1, help="World size")
    argparser.add_argument("--seed", type=int, default=0, help="If set to 0 then a random seed will be picked")
    argparser.add_argument("--device", type=str, default="1", help="GPU")

    # config model path, sample method,
    argparser.add_argument("--CASP14_test_set", type=int, default=1, help="Whether to use CASP14 as the test set")
    argparser.add_argument("--checkpoint_path", type=str, default="./params/model_110.pth", help="Path of model")
    argparser.add_argument("--ContGreedy", type=int, default=0, help="ContGreedy")

    # config the input/output process information
    argparser.add_argument("--pdb_path", type = str, required = True, default = "example/pdb/1acf.pdb",  help = "Input pdb path")
    argparser.add_argument("--chain", type=str, default="C", help= "Assigned chain")
    argparser.add_argument("--output_file",  default = None, help="Specify the output file")

    # config fix residues for some position
    # argparser.add_argument("--fix", type=int, default=0, help="Fixed position")
    # argparser.add_argument("--fix_positions", nargs = '+', help="Fixed poaition")
    # argparser.add_argument("--fix_native_pos", nargs='+', default=[], help="Assigned fixed position in native sequence")
    # argparser.add_argument("--fix_native_val", nargs='+', default=[], help="Assigned fixed value correspond to the position in native sequence")
    argparser.add_argument("--fix_file", type = str, help="Fixed position")

    
    
    # It is necessary to assign one and only one exist chain
    
    '''argparser.add_argument("--json_path", type = str, 
                          default = "/home/wangxinglong/MPNN_dataset/5U4Y_C.jsonl", 
                          help    = "Path of proteinA structure information file with the json format")'''
                          


    args = argparser.parse_args()
    
    test(args)
