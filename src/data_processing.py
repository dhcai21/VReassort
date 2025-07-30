import os
from Bio import SeqIO
import pandas as pd
import platform

def rename(args):
    """ data loading """
    df_input = pd.read_csv(args.input,keep_default_na=False)
    seg_name = df_input['segment'].values
    seg_path = df_input['path'].values
    seg_num = len(seg_name)
    """ rename processing """
    fasta_path = os.path.join(args.out,'fasta') 
    os.system(f'mkdir -p {fasta_path}')
    head_all = []
    for i in range(seg_num):
        f = open(f"{fasta_path}/{seg_name[i]}.fasta",'w')
        n = 0
        head = []
        for s in SeqIO.parse(seg_path[i],'fasta'):
            n = n + 1
            head.append(s.description)
            f.write(f">S{n}\n{s.seq}\n")
        head_all.append(head)
    """ save the mapping files """
    df_head = pd.DataFrame()
    df_head['Name'] = [f"S{i+1}" for i in range(len(head_all[0]))]
    for i in range(seg_num):
        df_head[seg_name[i]] = head_all[i]
    df_head.to_csv(os.path.join(args.out,'name_mapping.csv'),index=False)
    return seg_name,fasta_path

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
