import argparse
import os
import time
import sys
import multiprocessing as mp
mp.set_start_method("fork")

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from data_processing import rename, msa_and_tree
from feature_processing import *
from model import *
from identify_processing import *

# if __name__ == '__main__':
    #parse args
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True, help = "a csv file indicateing the files of segment sequences and their path.")
parser.add_argument('--out', type=str, required=False, default = "result" , help = "output folder (default: result).")
parser.add_argument('--thread', type=int, required=False, default = 8, help = "number of threads (default: 8).")

args = parser.parse_args()

# args = parser.parse_args(['--input', 'input.csv', '--out', 'result', '--thread', '8'])

print("**** Parameters ****")
print(args)
######### 1. Data Processing (data_processing.py) ###############

""" sequence rename """
print("**** Data Processing ****")
seg_name,fasta_path = rename(args)
print(seg_name)

""" MSA, TrimAl and FastTree """
tree_path = msa_and_tree(args,seg_name,fasta_path,current_dir)


######### 2. Feature Processing (feature_processing.py) ###############
print("**** Feature Processing ****")
""" Basic Feature """
print("Basic feature")
feature_path = call_single(args,seg_name,tree_path)
""" Pair feature """
print("Generating feature")
call_pair(args,seg_name,tree_path,feature_path)


######## 3. Prediction Processing (identification.py model.py) ############
print("**** Identification Processing ****")
""" Identification """
call_identify(args,seg_name,feature_path,current_dir)
print("**** Finish ****")





