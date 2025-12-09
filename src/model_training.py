import argparse
import os
import time
import sys
import random
import numpy as np
import pandas as pd
from types import SimpleNamespace
import torch
from torch.utils.data import DataLoader
import concurrent.futures
from model import Siamese_cnn, myDS_train

# import multiprocessing as mp
# mp.set_start_method("fork")

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_dir}")



def feature_processing(data):
    res1 = pd.DataFrame()
    res2 = pd.DataFrame()
    res1['na1'] = data['na1'].values
    res1['na2'] = data['na2'].values
    res1['s1'] = data['s1'].values
    res1['s2'] = data['s2'].values
    res1['r1'] = data['r1'].values
    res1['r2'] = data['r2'].values
    res1['label'] = data['label'].values
    res1['name'] = data['name'].values
    res2['na1'] = data['na3'].values
    res2['na2'] = data['na4'].values
    res2['s1'] = data['s3'].values
    res2['s2'] = data['s4'].values
    res2['r1'] = data['r3'].values
    res2['r2'] = data['r4'].values
    res2['label'] = data['label'].values
    res2['name'] = data['name'].values
    res = pd.concat([res1,res2])
    return res


def training(config): 

    """ Cuda Check """
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(f"Using {device}")
    """ Loading training dataset """
    train_data = pd.read_csv(config.feature)
    train_data = feature_processing(train_data)
    
    trainDS=myDS_train(train_data,config)
    print(f"Training data size: {len(trainDS)}")
    
    """ Model setting """
    siamese = Siamese_cnn(config).to(device)
    print(f"Network structure:")
    print(siamese)
    
    criterion = torch.nn.MSELoss().to(device)
    
    # optimizer
    learning_rate = config.lr
    optimizer = torch.optim.Adam(filter(lambda x: x.requires_grad, siamese.parameters()), lr=learning_rate)
    
    # save the model
    ckpt_path = os.path.join(config.output)
    """ Train the model """
    train_loss = []
    train_s_time = time.time()
    epoch_len = len(str(config.epochs))
    for epoch in range(config.epochs):
        # loss
        train_loss = []
        
        # dataloader
        train_dataloader = DataLoader(dataset=trainDS, shuffle=True, num_workers=0, batch_size=config.batch_size)

        for idx, data in enumerate(train_dataloader, 0):

            # Load the data
            s1, s2,label,_ = data
            s1 = s1.to(device)
            s2 = s2.to(device)
            label = label.to(device)
            
            # clear gradients
            optimizer.zero_grad()

            # Calculate the probability
            pro = siamese(s1, s2)
            pro = pro.squeeze(0)

            # loss backward
            loss = criterion(pro, label)
            loss.backward()
            optimizer.step()
            train_loss.append(loss.data.cpu())
        
        # Record the loss at each epoch
        ave_loss = np.average(train_loss)
        print(f'Epoch [{epoch + 1:>{epoch_len}}/{config.epochs:>{epoch_len}}] ' + f'train_loss: {ave_loss:.5f} ')
    """ training time """
    train_e_time = time.time()
    print(f"Training time: {(train_e_time - train_s_time) / 60 :.2f} min.")
    return siamese




def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a neural network model."
    )

    # --- input/output files ---
    parser.add_argument(
        "--feature",
        type=str,
        required=True,
        help="Path to the feature file (e.g., train_feature.csv)."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="trained.pt",
        help="Output file name for the trained model (default: trained.pt)."
    )

    # --- network parameters ---
    parser.add_argument(
        "--input_size",
        type=int,
        default=40,
        help="Input feature dimension (default: 40)."
    )
    parser.add_argument(
        "--h1_dim",
        type=int,
        default=128,
        help="Number of hidden units in the first fully connected layer (default: 128)."
    )
    parser.add_argument(
        "--h2_dim",
        type=int,
        default=32,
        help="Number of hidden units in the second fully connected layer (default: 32)."
    )
    parser.add_argument(
        "--cnn_size",
        type=int,
        default=5,
        help="Kernel size for CNN layers (default: 5)."
    )
    parser.add_argument(
        "--cnn_step",
        type=int,
        default=1,
        help="Stride for CNN layers (default: 1)."
    )

    # --- training hyperparameters ---
    parser.add_argument(
        "--epochs",
        type=int,
        default=100,
        help="Number of training epochs (default: 100)."
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=64,
        help="Mini-batch size (default: 64)."
    )
    parser.add_argument(
        "--lr",
        type=float,
        default=0.001,
        help="Learning rate (default: 0.001)."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=21,
        help="Random seed for reproducibility (default: 21)."
    )

    # ==== Feature parameters ====
    parser.add_argument(
        '--n_thres', type=int, required=False, default=20,
        help="Minimum normalization cap boundary (default: 20)."
    )
    parser.add_argument(
        '--theta', type=float, required=False, default=0.01,
        help="Branch length distance threshold for increment calculation (default: 0.01)."
    )
    parser.add_argument(
        '--B', type=float, required=False, default=5,
        help="Balance factor for the branch length distance component (default: 5)."
    )
    parser.add_argument("--vec_size", type=int,default=40,
                        help="Feature vector size. (default: 40)")

    return parser.parse_args()



if __name__ == '__main__':
    """ loading parameters """
    config = parse_args()
    seed = config.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    """ network """
    network = SimpleNamespace()
    network.h1_dim = config.h1_dim
    network.h2_dim = config.h2_dim
    network.cnn_size = config.cnn_size
    network.cnn_step = config.cnn_step
    network.input_size = config.input_size
    """ training """
    print(config)
    siamese = training(config)
    siamese = siamese.cpu()
    train_dict = {
        'model': siamese.state_dict(),
        'network': network
    }
    torch.save(train_dict, config.output)
    exit()

