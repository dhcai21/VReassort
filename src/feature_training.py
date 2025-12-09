import argparse
import os
import time
import sys
import multiprocessing as mp
mp.set_start_method("fork")

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_dir}")

# from data_processing import *
from feature_processing import *



def parse_args():
    parser = argparse.ArgumentParser(
        description="Process subfolders and generate query list for data processing."
    )

    # --- Folder and file parameters ---
    parser.add_argument("--main_folder", type=str, required=True,
                        help="Path to the main folder (main_folder/subfolder/samplex/seg1.treefile).")
    parser.add_argument("--output", type=str, default="train_feature.csv",
                        help="Output file name. (default: train_feature.csv)")
    parser.add_argument("--thread", type=int, default=8,
                        help="Number of threads to use.")

    # --- Processing parameters ---
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
    parser.add_argument("--vec_size", type=int,default=40,
                        help="Feature vector size. (default: 40)")

    # --- Training or task parameters ---
    parser.add_argument("--times", type=int, default=5,
                        help="scale of Negative strains. (default 5, Positive : Negative = 1:5)")

    return parser.parse_args()



if __name__ == '__main__':
    """ loading parameters """
    # rng = np.random.default_rng()
    config = parse_args()
    print(config)
    file_out = config.output
    s1,s2,s3,s4,label,name,na1,na2,na3,na4,r1,r2,r3,r4 = run_train(config)
    """ save the leaves feature """ 
    data = pd.DataFrame()
    data['na1'] = na1
    data['na2'] = na2
    data['na3'] = na3
    data['na4'] = na4
    data['s1'] = s1
    data['s2'] = s2
    data['s3'] = s3
    data['s4'] = s4
    data['r1'] = r1
    data['r2'] = r2
    data['r3'] = r3
    data['r4'] = r4
    data['label'] = label
    data['name'] = name
    data.to_csv(f'{file_out}',index=False)


