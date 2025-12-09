#!/usr/bin/env python
import argparse
import os
import time
import sys
import multiprocessing as mp
mp.set_start_method("fork")

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_dir}/src")

from data_processing import *
from feature_processing import *
from identify_processing import *

import argparse
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="Script for processing segment sequences or trees.")

# ==== Basic arguments ====
parser.add_argument(
    '--input', type=str, required=True,
    help="CSV file indicating the input files of segment sequences/trees and their paths."
)
parser.add_argument(
    '--type', type=str, required=True, default='seq',
    help="Indicator of the type of input data: 'tree' or 'seq' (default: seq)."
)
parser.add_argument(
    '--out', type=str, required=False, default="result",
    help="Output folder (default: result)."
)
parser.add_argument(
    '--model', type=str, required=False, default=f"{current_dir}/src/trained.pt",
    help="Path to trained model file (default: src/trained.pt)."
)
parser.add_argument(
    '--thread', type=int, required=False, default=8,
    help="Number of threads to use (default: 8)."
)
parser.add_argument(
    '--reshuffle', action='store_true',
    help="Reshuffle sequence order to generate trees with slight topology variation."
)

# ==== Feature parameters ====
parser.add_argument(
    '--n_thres', type=int, required=False, default=20,
    help="Minimum normalization cap boundary (default: 20)."
)
parser.add_argument(
    '--b_thres1', type=float, required=False, default=0.0003,
    help="Lower branch collapsing threshold (default: 0.0003)."
)
parser.add_argument(
    '--b_thres2', type=float, required=False, default=0.001,
    help="Upper branch collapsing threshold (default: 0.001)."
)
parser.add_argument(
    '--theta', type=float, required=False, default=0.01,
    help="Branch length distance threshold for increment calculation (default: 0.01)."
)
parser.add_argument(
    '--B', type=float, required=False, default=5,
    help="Balance factor for the branch length distance component (default: 5)."
)

# ==== Clade-related parameters ====
parser.add_argument(
    '--min_clade_size', type=int, required=False, default=5,
    help="Minimum clade size (default: 5)."
)
parser.add_argument(
    '--max_clade_size', type=int, required=False, default=150,
    help="Maximum clade size (default: 150)."
)
parser.add_argument(
    '--jd', type=float, required=False, default=0.8,
    help="Minimum Jaccard index between clades for similarity consideration (default: 0.8)."
)
parser.add_argument(
    '--max_clade_diff', type=int, required=False, default=10,
    help="Maximum allowed number of different strains between clades (default: 10)."
)

args = parser.parse_args()

# args = parser.parse_args(['--input', 'input.csv', '--out', 'result', '--thread', '8'])

print("**** Parameters ****")
print(args)
######### 1. Data Processing (data_processing.py) ###############
print("**** Data Processing ****")
seg_name,tree_path = data_checking_processing(args)

######### 2. Feature Processing (feature_processing.py) ###############
print("**** Feature Processing ****")
""" Basic Feature """
print("Basic feature")
feature_path = call_single(args,seg_name,tree_path)
""" Pair feature """
print("Generating feature")
call_pair(args,seg_name,tree_path,feature_path)

# """ for testing """
# seg_name = df_input = pd.read_csv(args.input,keep_default_na=False)
# seg_name = df_input['segment'].values
# feature_path = os.path.join(args.out,'feature')


######## 3. Prediction Processing (identification.py model.py) ############
print("**** Identification Processing ****")
""" Identification """
call_identify(args,seg_name,feature_path)
print("**** Finish ****")





