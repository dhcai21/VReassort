#coding=utf-8
import os
import random
import numpy as np
import pandas as pd
import sys
import copy
from collections import Counter
from types import SimpleNamespace
import torch
from torch.utils.data import DataLoader
from ete3 import Tree
import multiprocessing as mp
import concurrent.futures
from tqdm import tqdm


current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_dir}")

from model import Siamese_cnn, myDS


""" main: identification """

def call_identify(args,seg_name,feature_path):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    result_path = os.path.join(args.out,'identification')
    os.system(f"mkdir -p {result_path}")
    query_list = []
    df_result_summary = pd.DataFrame()
    pair1 = []
    pair2 = []
    leaf_list = []
    clade_list = []
    for i in range(len(seg_name)-1):
        for j in range(i+1,len(seg_name)):
            leaf_file = f"{feature_path}/{seg_name[i]}_{seg_name[j]}_leaf.csv"
            clade_file = f"{feature_path}/{seg_name[i]}_{seg_name[j]}_clade.csv"
            para_model = SimpleNamespace(leaf_file=leaf_file,clade_file=clade_file,out=result_path)
            para_model.prefix=f"{seg_name[i]}_{seg_name[j]}"
            para_model.n_thres = args.n_thres
            para_model.B = args.B
            para_model.theta = args.theta
            para_model.ckpt = args.model
            query_list.append(para_model)
            out_leaf = os.path.join(result_path, f"{para_model.prefix}_leaf.csv")
            out_clade = os.path.join(result_path, f"{para_model.prefix}_clade.csv")
            pair1.append(seg_name[i])
            pair2.append(seg_name[j])
            leaf_list.append(out_leaf)
            clade_list.append(out_clade)
    df_result_summary['segi'] = pair1
    df_result_summary['segj'] = pair2
    df_result_summary['leaf'] = leaf_list
    df_result_summary['clade'] = clade_list
    out_summary = os.path.join(args.out, f"result_files.csv")
    df_result_summary.to_csv(out_summary,index=False)
    if len(query_list) == 1:
        identify(query_list[0])
    else:
        with concurrent.futures.ProcessPoolExecutor(min(len(query_list),args.thread)) as executor:
            list(tqdm(executor.map(identify, query_list), total=len(query_list)))
    return len(query_list)



def feature_processing(data,mode):
    if mode ==1:
        res = pd.DataFrame()
        res['na1'] = data['na1'].values
        res['na2'] = data['na2'].values
        res['s1'] = data['s1'].values
        res['s2'] = data['s2'].values
        res['r1'] = data['r1'].values
        res['r2'] = data['r2'].values
        res['name'] = data['name'].values
    elif mode == 2:
        res = pd.DataFrame()
        res['na1'] = data['na3'].values
        res['na2'] = data['na4'].values
        res['s1'] = data['s3'].values
        res['s2'] = data['s4'].values
        res['r1'] = data['r3'].values
        res['r2'] = data['r4'].values
        res['name'] = data['name'].values
    return res

def prediction(data, siamese, vec_size, config, mask_strain=[]):
    device = next(siamese.parameters()).device
    testDS = myDS(data, vec_size,config, mask_strain)
    test_dataloader = DataLoader(dataset=testDS, num_workers=0, batch_size=1)
    result, pro, test_name = [], [], []
    # siamese.eval()
    with torch.no_grad():
        for idx, (s1, s2, name) in enumerate(test_dataloader):
            s1 = s1.to(device)
            s2 = s2.to(device)
            output = siamese(s1, s2)
            output = output.squeeze(0)
            pro.append(output.cpu().numpy())
            test_name.append(name[0])
    pro = np.array(pro)
    test_name = np.array(test_name)
    result = pd.DataFrame({'score': pro, 'name': test_name})
    return result


def two_pro(data1,data2):
    name1 = data1['name'].values
    name2 = data2['name'].values
    pro1 = data1['score'].values
    pro2 = data2['score'].values
    pro = {n:[0,0] for n in name1}
    for i in range(len(name1)):
        n1 = name1[i]
        n2 = name2[i]
        p1 = pro1[i]
        p2 = pro2[i]
        pro[n1][0] = p1
        pro[n2][1] = p2
    final_pro = []
    for key in pro.keys():
        p1 = pro[key][0]
        p2 = pro[key][1]
        final_pro.append(p1+p2)
    p2 = [pro[key][1] for key in pro.keys()]
    data = pd.DataFrame()
    data['score'] = final_pro
    data['s1'] = pro1
    data['s2'] = p2
    data['name'] = name1
    return data


def clade_neighbors_match(feature,reassort_index):
    """ load the basic data of the predicted reassortant clades """
    clades = [i.split(',') for i in feature['share_leaves'].values[reassort_index]]
    na1 = [i.split(',') for i in feature['na1'].values[reassort_index]]
    na2 = [i.split(',') for i in feature['na2'].values[reassort_index]]
    r1 = [i.split(',') for i in feature['r1'].values[reassort_index]]
    r2 = [i.split(',') for i in feature['r2'].values[reassort_index]]
    r3 = [i.split(',') for i in feature['r3'].values[reassort_index]]
    r4 = [i.split(',') for i in feature['r4'].values[reassort_index]]
    s1 = [i.split(',') for i in feature['s1'].values[reassort_index]]
    s2 = [i.split(',') for i in feature['s2'].values[reassort_index]]
    s3 = [i.split(',') for i in feature['s3'].values[reassort_index]]
    s4 = [i.split(',') for i in feature['s4'].values[reassort_index]]
    s_vec = [s1,s2,s3,s4]
    r_vec = [r1,r2,r3,r4]
    num_reassort = len(reassort_index)
    """ check if the clade i is in the top-20 neighbors of j (only one of the segments) """
    match_info = {i:{} for i in range(num_reassort)}
    for i in range(num_reassort):
        c = clades[i]
        for j in range(num_reassort):
            cj = clades[j]
            if i ==j:
                continue
            if len (set(cj).intersection(set(c))) !=0:
                match_info[i][j] = 0
                continue
            n1 = na1[j]
            n2 = na2[j]
            m1 = len(set(c).intersection(set(n1[:20])))
            m2 = len(set(c).intersection(set(n2[:20])))
            if m1 ==0 and m2 ==0:
                match_info[i][j] = 0
            elif m1 >= m2:
                match_info[i][j] = 1
            else:
                match_info[i][j] = 2
    return clades,match_info,s_vec,r_vec


def compare_length(idx,m,s_vec):
    num = int(len(s_vec[0][idx])/2)
    v1 = np.array([float(i) for i in s_vec[0][idx]])[num:]
    v2 = np.array([float(i) for i in s_vec[1][idx]])[num:]
    v3 = np.array([float(i) for i in s_vec[2][idx]])[num:]
    v4 = np.array([float(i) for i in s_vec[3][idx]])[num:]
    v1 = v1/max(v1)
    v2 = v2/max(v2)
    v3 = v3/max(v3)
    v4 = v4/max(v4)
    if m == 1:
        val = abs(v2-v4)
    else:
        val = abs(v1-v3)
    return val



def real_mask_clades(match_info,s_vec,r_vec):
    candidate_pair = [] # for checking the comapre score
    mask_clades = [] # the 'real' reassortant one 
    num_reassort = len(match_info)
    fake_clades = {i:[] for i in range(num_reassort)} # the corresponding fp reassortment affected by the real one i
    pair = {i:{j:0 for j in range(num_reassort)} for i in range(num_reassort)}
    for i in range(num_reassort-1):
        for j in range(i+1,num_reassort):
            mi = match_info[j][i] # which seg of i have j
            mj = match_info[i][j] # which seg of j have i
            if mi != 0 and mj != 0:
                pair[i][j] = 1
                pair[j][i] = 1
                vali = compare_length(i,mi,s_vec)
                valj = compare_length(j,mj,s_vec)
                ratiol = np.count_nonzero((vali-valj)[:20]>0)/20
                flagl = 1 if ratiol > 0.5 else 0
                if flagl==1:
                    mask_clades.append(i)
                    fake_clades[i].append(j)
                else:
                    mask_clades.append(j)
                    fake_clades[j].append(i)
    frequency = Counter(mask_clades)
    sorted_unique = sorted(frequency.keys(), key=lambda x: frequency[x],reverse=True)
    final_clades = []
    fp_clades = []
    for i in sorted_unique:
        if i in fp_clades:
            continue
        else:
            final_clades.append(i)
            fp_clades = fp_clades + fake_clades[i]
    print(final_clades)
    for i in range(num_reassort):
        if sum(pair[i].values()) ==0:
            final_clades.append(i)
    print(final_clades)
    return final_clades,fake_clades


def subclade(all_clade,clade2):
    res = []
    for i in range(len(all_clade)):
        flag = 0
        for j in range(len(clade2)):
            c1 = all_clade[i]
            c2 = clade2[j]
            if len(c1) > len(c2): # not include the reassorted one
                continue
            share = set(c1).intersection(set(c2))
            if len(share) == len(c1):
                flag = 1
        if flag == 0:
            res.append(i)
    return res



def identify(config):
    """ loading config """
    thres_mask =  0.5
    leaf_for_clade_mask = 0.99
    ckpt_path = config.ckpt
    leaf_file = config.leaf_file
    clade_file = config.clade_file
    out_leaf = os.path.join(config.out, f"{config.prefix}_leaf.csv")
    out_clade = os.path.join(config.out, f"{config.prefix}_clade.csv")
    """ loading the trained model """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    info = torch.load(ckpt_path, map_location=device,weights_only=False)
    print("Running on:", device)
    network = info['network']
    model_parameter = info['model']
    siamese = Siamese_cnn(network)
    siamese.load_state_dict(model_parameter)
    siamese.to(device)           # <— move model to the same device as data
    siamese.eval()               # <— set to evaluation mode for inference
    n_thres = config.n_thres
    """ load the leaf feature """
    leaf_feature = pd.read_csv(leaf_file)
    ref_leaves = leaf_feature['name'].values
    """ leaves : segment 1"""
    mode = 1
    inputs = feature_processing(leaf_feature,mode)
    # first prediction
    mask_strains_l1 = []
    leaf11_result = prediction(inputs,siamese,network.input_size, config, mask_strains_l1)
    # second prediction
    mask_strains_l1 = list(ref_leaves[leaf11_result['score'].values>thres_mask]) # leaf masking 1
    leaf12_result = prediction(inputs,siamese,network.input_size, config, mask_strains_l1)
    """ leaves segment 2 """
    mode = 2
    inputs = feature_processing(leaf_feature,mode)
    # first prediction
    mask_strains_l2 = []
    leaf21_result = prediction(inputs,siamese,network.input_size, config, mask_strains_l2)
    # second prediction
    mask_strains_l2 = list(ref_leaves[leaf21_result['score'].values>thres_mask]) # leaf masking 2
    leaf22_result = prediction(inputs,siamese,network.input_size, config, mask_strains_l2)
    """ leaves: two pro """
    two_pro_leaf = two_pro(leaf12_result,leaf22_result)
    """ *************** clade prediction************** """
    """ loda the clades """
    clade_feature = pd.read_csv(clade_file)
    if clade_feature['na1'][0]== -1:
        print("No matched clades")
        exit()
    all_clades = [i.split(',') for i in clade_feature['share_leaves'].values]
    mask_strains_l = list(ref_leaves[two_pro_leaf['score'].values>leaf_for_clade_mask])
    """ clade: segment 1 """
    mode = 1
    inputs = feature_processing(clade_feature,mode)
    mask_strains_c1 = mask_strains_l
    clade11_result = prediction(inputs,siamese,network.input_size, config, mask_strains_c1)
    """ clade: segment 2 """
    mode = 2
    mask_strains_c2 = mask_strains_l
    inputs = feature_processing(clade_feature,mode)
    clade21_result = prediction(inputs,siamese,network.input_size, config, mask_strains_c2)
    """ clade: two pro """
    two_pro_clade = two_pro(clade11_result,clade21_result)
    info_key = list(clade_feature.keys())[14:]
    for key in info_key:
        two_pro_clade[key] = clade_feature[key] #first clade
    """ filtering clades """
    pro = two_pro_clade['score'].values/2
    clade_mask_thres = 0.5
    reassort_index = np.array([i for i in range(len(pro))])[pro>thres_mask]
    if len(reassort_index) == 0:
        two_pro_leaf['score'] = two_pro_leaf['score']/2
        two_pro_leaf.to_csv(out_leaf,index=False)
        two_pro_clade['score'] = two_pro_clade['score']/2
        two_pro_clade.to_csv(out_clade,index=False)
        print("**** Finish ****")
        exit()
    candidate_clades, match_info, s_vec, r_vec = clade_neighbors_match(clade_feature,reassort_index) # association between clades
    final_clades,fake_clades = real_mask_clades(match_info,s_vec,r_vec)
    mask_strains_c = []
    for i in final_clades:
        mask_strains_c = mask_strains_c + candidate_clades[i]
    """ re-predict leaves and clades"""
    """ leaves """
    replace_index = [i for i in range(len(ref_leaves)) if ref_leaves[i] not in mask_strains_c]
    # segment 1
    mode = 1
    inputs = feature_processing(leaf_feature,mode)
    mask_strains_l1 = mask_strains_c + mask_strains_l1
    leaf13_result = prediction(inputs,siamese,network.input_size, config, mask_strains_l1)
    leaf12_result.loc[replace_index] = leaf13_result.loc[replace_index] # keep the original results of leaves in clades 
    # segment 2
    mode = 2
    inputs = feature_processing(leaf_feature,mode)
    mask_strains_l2 = mask_strains_c + mask_strains_l2
    leaf23_result = prediction(inputs,siamese,network.input_size, config, mask_strains_l2)
    leaf22_result.loc[replace_index] = leaf23_result.loc[replace_index]
    """ leaves: two pro """
    two_pro_leaf = two_pro(leaf12_result,leaf22_result)
    """ clades """
    mask_strains_l = list(ref_leaves[two_pro_leaf['score'].values>leaf_for_clade_mask])
    reassort_clades = [candidate_clades[i] for i in final_clades]
    replace_index = subclade(all_clades,reassort_clades)
    # segment 1
    mode = 1
    inputs = feature_processing(clade_feature,mode)
    mask_strains_c1 = mask_strains_l + mask_strains_c
    clade12_result = prediction(inputs,siamese,network.input_size, config, mask_strains_c1)
    clade11_result.loc[replace_index] = clade12_result.loc[replace_index] 
    # segment 2
    mode = 2
    inputs = feature_processing(clade_feature,mode)
    mask_strains_c1 = mask_strains_l + mask_strains_c
    clade22_result = prediction(inputs,siamese,network.input_size, config, mask_strains_c2)
    clade21_result.loc[replace_index] = clade22_result.loc[replace_index]
    """ clade: two pro """
    two_pro_clade = two_pro(clade11_result,clade21_result)
    info_key = list(clade_feature.keys())[14:]
    for key in info_key:
        two_pro_clade[key] = clade_feature[key]
    """ save the result """
    two_pro_leaf['score'] = two_pro_leaf['score']/2
    two_pro_leaf.to_csv(out_leaf,index=False)
    two_pro_clade['score'] = two_pro_clade['score']/2
    two_pro_clade.to_csv(out_clade,index=False)
    
    

    
    
