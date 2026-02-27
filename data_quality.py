#!/usr/bin/env python3
import os
import sys
import random
import itertools
import statistics
import argparse
import platform
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
import dendropy
from dendropy.calculate import treecompare

def subsample_once(alignment, fraction):
    aln_length = alignment.get_alignment_length()
    n_cols_to_keep = int(aln_length * fraction)
    if n_cols_to_keep == 0:
        raise ValueError("Fraction too small — no columns selected.")
    selected_columns = sorted(random.sample(range(aln_length), n_cols_to_keep))
    new_records = []
    for record in alignment:
        new_seq = "".join(record.seq[i] for i in selected_columns)
        new_record = record[:]
        new_record.seq = Seq(new_seq)
        new_records.append(new_record)
    return MultipleSeqAlignment(new_records)

def compute_rf_statistics(tree_files):
    if len(tree_files) < 2:
        return None
    tns = dendropy.TaxonNamespace()
    trees = [dendropy.Tree.get(path=path, schema="newick", taxon_namespace=tns) for path in tree_files]
    rf_distances = []
    for t1, t2 in itertools.combinations(trees, 2):
        dist = treecompare.unweighted_robinson_foulds_distance(t1, t2)
        rf_distances.append(dist)
    n_leaves = len(trees[0].leaf_nodes())
    
    # Maximum possible RF distance for unrooted trees is 2 * (n_leaves - 3)
    max_rf = 2 * (n_leaves - 3)
    mean_rf = statistics.mean(rf_distances)
    median_rf = statistics.median(rf_distances)
    normalized_mean_rf = mean_rf / max_rf if max_rf > 0 else 0
    normalized_median_rf = median_rf / max_rf if max_rf > 0 else 0
    return mean_rf, normalized_mean_rf, median_rf, normalized_median_rf

def process_msa(msa_path, threads, fraction=0.8, n_reps=5, seed=21):
    if seed is not None:
        random.seed(seed)
        
    if not os.path.exists(msa_path):
        print(f"Error: File not found: {msa_path}")
        sys.exit(1)
        
    folder = os.path.dirname(os.path.abspath(msa_path))
    msa_name = os.path.basename(msa_path)
    
    print(f"Reading alignment: {msa_path}")
    alignment = AlignIO.read(msa_path, "fasta")
    
    output_dir = os.path.join(folder, "subsample")
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(msa_name)[0]
    
    tree_files = []
    
    # Set OpenMP threads for FastTreeMP
    os.system(f"export OMP_NUM_THREADS={threads}")
    
    # Get current directory of the script to locate FastTreeMP binaries
    current_dir = os.path.dirname(os.path.abspath(__file__))
    system_name = platform.system()
    
    for i in range(1, n_reps + 1):
        print(f"Generating subsample {i}/{n_reps}...")
        new_alignment = subsample_once(alignment, fraction)
        sub_msa = os.path.join(output_dir, f"{base_name}_sub{i}.fa")
        AlignIO.write(new_alignment, sub_msa, "fasta")
        
        tree_file = os.path.join(output_dir, f"{base_name}_sub{i}.tree")
        
        # OS-specific FastTreeMP command
        if system_name == 'Darwin':
            cmd = f"{current_dir}/src/FastTreeMP_mac -quiet -cat 4 -nt -gtr -gamma {sub_msa} > {tree_file}"
        else:
            cmd = f"{current_dir}/src/FastTreeMP_linux -quiet -cat 4 -nt -gtr -gamma {sub_msa} > {tree_file}"
        os.system(cmd)
        
        if os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
            tree_files.append(tree_file)
        else:
            print(f"Warning: Tree file {tree_file} was not generated correctly.")
            
    print("Computing RF statistics...")
    rf_result = compute_rf_statistics(tree_files)
    return rf_result

def main():
    parser = argparse.ArgumentParser(description="Calculate normalized RF distance for a given MSA file.")
    parser.add_argument("-msa", required=True, help="Path to the input MSA file (FASTA format)")
    parser.add_argument("-thread", type=int, default=8, help="Number of threads for FastTreeMP (default: 8)")
    parser.add_argument("-fraction", type=float, default=0.8, help="Fraction of columns to subsample (default: 0.8)")
    parser.add_argument("-n_reps", type=int, default=5, help="Number of subsample replicates (default: 5)")
    parser.add_argument("-seed", type=int, default=21, help="Random seed for reproducibility (default: 21)")
    
    args = parser.parse_args()
    
    results = process_msa(
        msa_path=args.msa, 
        threads=args.thread, 
        fraction=args.fraction, 
        n_reps=args.n_reps, 
        seed=args.seed
    )
    
    if results is not None:
        mean_rf, norm_mean_rf, median_rf, norm_median_rf = results
        print("\n=== Results ===")
        print(f"MSA File             : {args.msa}")
        print(f"Mean RF              : {mean_rf:.4f}")
        print(f"Normalized Mean RF   : {norm_mean_rf:.4f}")
        print("=================\n")
    else:
        print("Failed to compute RF statistics. Ensure that FastTreeMP ran successfully and generated at least 2 trees.")

if __name__ == "__main__":
    main()
