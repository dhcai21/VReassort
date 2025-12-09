import os
from Bio import SeqIO
import pandas as pd
import platform
import dendropy
import numpy as np
import random

def data_checking_processing(args):
    if args.type == 'seq':
        """ sequence rename """
        seg_name,head_all,seq_all = load_seq(args)
        print("Checking the numbers and names of strains")
        num_flag, name_flag = check_name(head_all, seg_name)
        if num_flag:
            print("Number checking: Pass")
            seg_name,fasta_path = rename(args,name_flag,seg_name,head_all,seq_all)
        else:
            print("Error: strain numbers are inconsistent across segments")
            exit(1)
        print(seg_name)
        """ MSA, TrimAl and FastTree """
        tree_path = msa_and_tree(args,seg_name,fasta_path)
    elif args.type == 'tree':
        """ sequence rename """
        seg_name,head_all = load_tree(args)
        tree_path = os.path.join(args.out,'tree')
        print("Checking the numbers and names of strains")
        num_flag, name_flag = check_name(head_all, seg_name)
        if num_flag:
            print("Number checking: Pass")
            if not name_flag:
                print("Name checking: Failed")
                print("Error: The leaves' names on the segment trees are not consistent")
                exit(1)
        else:
            print("Error: strain numbers are inconsistent across segments")
            exit(1)
    else:
        print("Error: Please specify the type of the input data (seq or tree)")
        exit(1)
    return seg_name,tree_path


def load_seq(args):
    """ data loading """
    df_input = pd.read_csv(args.input,keep_default_na=False)
    seg_name = df_input['segment'].values
    seg_path = df_input['path'].values
    seg_num = len(seg_name)
    """ sequences processing """
    head_all = []
    seq_all = []
    for i in range(seg_num):
        head = []
        seq = []
        for s in SeqIO.parse(seg_path[i],'fasta'):
            head.append(s.description)
            seq.append(s.seq)
        head_all.append(head)
        seq_all.append(seq)
    return seg_name,head_all,seq_all

def check_name(head_all, seg_name):
    """
    head_all: list of lists of headers (one inner list per file)
    seg_name: list of file/segment names (same length as head_all)
    Returns: num_flag, name_flag
    """
    num_flag = True
    name_flag = True
    lengths = [len(h) for h in head_all]
    if not lengths:
        return num_flag, name_flag
    first_len = lengths[0]
    # check numbers
    if any(l != first_len for l in lengths[1:]):
        num_flag = False
        print("Strain count per file:")
        for i, l in enumerate(lengths):
            label = seg_name[i] if seg_name is not None else f"file_{i+1}"
            print(f"{label}: {l}")
        return num_flag, name_flag
    # check headers (ignoring order)
    ref_set = set(head_all[0])
    for i in range(1, len(head_all)):
        cur_set = set(head_all[i])
        if cur_set != ref_set:
            name_flag = False
    return num_flag, name_flag


def rename(args,name_flag,seg_name,head_all,seq_all):
    """
    name_flag: consistent names across samples, True or False
    head_all: the headers
    seq_all: the sequences
    """
    """ rename processing """
    seq_num = len(head_all[0])
    seg_num = len(head_all)
    fasta_path = os.path.join(args.out,'fasta') 
    os.system(f'mkdir -p {fasta_path}')
    print("Synchronizing sequence names")
    if name_flag:
        """ mapping pre-processing """
        print("Consistent: mapping by names")
        mapping = {head_all[0][i]:f"S{i+1}" for i in range(seq_num)}
        mapping2 = {f"S{i+1}":head_all[0][i] for i in range(seq_num)}
        for i in range(seg_num):
            f = open(f"{fasta_path}/{seg_name[i]}.fasta",'w')
            index = [j for j in range(seq_num)]
            if args.reshuffle:
                random.shuffle(index)
            for j in index:
                head = mapping[head_all[i][j]]
                seq = seq_all[i][j]
                f.write(f">{head}\n{seq}\n")
    else:
        print("Inconsistent: mapping by sequence orders")
        for i in range(seg_num):
            f = open(f"{fasta_path}/{seg_name[i]}.fasta",'w')
            head_list = []
            seq_list = []
            index = [j for j in range(seq_num)]
            if args.reshuffle:
                random.shuffle(index)
            for j in index:
                head = f"S{j+1}"
                seq = seq_all[i][j]
                f.write(f">{head}\n{seq}\n")
    """ save the mapping files """
    df_head = pd.DataFrame()
    df_head['Name'] = [f"S{i+1}" for i in range(seq_num)]
    if name_flag:
        head = [mapping2[s] for s in df_head['Name'].values]
        for i in range(seg_num):
            df_head[seg_name[i]] = head
    else:
        for i in range(seg_num):
            df_head[seg_name[i]] = head_all[i]
    df_head.to_csv(os.path.join(args.out,'name_mapping.csv'),index=False)
    return seg_name,fasta_path

def load_tree(args):
    """ data loading """
    df_input = pd.read_csv(args.input,keep_default_na=False)
    seg_name = df_input['segment'].values
    seg_path = df_input['path'].values
    seg_num = len(seg_name)
    """ copy the trees """
    tree_path = os.path.join(args.out,'tree')
    os.system(f"mkdir -p {tree_path}")
    """ sequences processing """
    head_all = []
    for i in range(seg_num):
        os.system(f"cp {seg_path[i]} {tree_path}/{seg_name[i]}.treefile")
        t = dendropy.Tree.get(path=seg_path[i],schema='newick')
        leaves = np.array([i for i in t.leaf_nodes() if i.is_leaf])
        head = [i.taxon.label for i in leaves]
        head_all.append(head)
    return seg_name,head_all


def msa_and_tree(args,seg_name,fasta_path):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    system_name = platform.system()
    tree_path = os.path.join(args.out,'tree')
    os.system(f"mkdir -p {tree_path}")
    for i in range(len(seg_name)):
        print(f"{i+1}. Segment:  {seg_name[i]}")
        seq_file = f"{fasta_path}/{seg_name[i]}.fasta"
        msa_file = f"{fasta_path}/{seg_name[i]}_msa.fasta"
        trim_file = f"{fasta_path}/{seg_name[i]}_trimed.fasta"
        """ MAFFT """
        print(f"{seg_name[i]}: Multiple Sequence Alignment (MAFFT)")
        os.system(f"mafft --quiet --thread {args.thread} {seq_file} > {msa_file}")
        """ Trimal """
        print(f"{seg_name[i]}: MSA Trimming (TrimAl)")
        os.system(f"trimal -in {msa_file} -out {trim_file} -gt 0.9 -cons 60") # potential para
        """ FastTree """
        print(f"{seg_name[i]}: Phylogenetic Tree Construction (FastTree)")
        os.system(f"export OMP_NUM_THREADS={args.thread}")
        if system_name == 'Darwin':
            os.system(f"{current_dir}/FastTreeMP_mac -quiet -cat 4 -nt -gtr -gamma {trim_file} > {tree_path}/{seg_name[i]}.treefile")
        else:
            os.system(f"{current_dir}/FastTreeMP_linux -quiet -cat 4 -nt -gtr -gamma {trim_file} > {tree_path}/{seg_name[i]}.treefile")
    return tree_path
