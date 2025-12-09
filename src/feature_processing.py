import os
from ete3 import Tree
import dendropy
import multiprocessing as mp
import concurrent.futures
from tqdm import tqdm
import numpy as np
import copy
import pandas as pd
import pickle
from scipy.stats import rankdata
from types import SimpleNamespace

import argparse



""" main: calling single feature """
def call_single(args,seg_name,tree_path):
    feature_path = os.path.join(args.out,'feature') 
    os.system(f"mkdir -p {feature_path}")
    ref_tree = f"{tree_path}/{seg_name[0]}.treefile"
    query_list = []
    for i in range(len(seg_name)):
        tree_file = f"{tree_path}/{seg_name[i]}.treefile"
        out_pic = f"{feature_path}/{seg_name[i]}.pickle"
        b_thres1 = args.b_thres1
        b_thres2 = args.b_thres2
        para_single = SimpleNamespace(tree1=ref_tree,tree2=tree_file,output=out_pic,b_thres1=b_thres1,b_thres2=b_thres2)
        query_list.append(para_single)
    with concurrent.futures.ProcessPoolExecutor(min(len(query_list),args.thread)) as executor:
        list(tqdm(executor.map(feature_single, query_list), total=len(query_list)))
    return feature_path


""" main: calling pair feature """
def call_pair(args,seg_name,tree_path,feature_path):
    query_list = []
    for i in range(len(seg_name)-1):
        for j in range(i+1,len(seg_name)):
            tree1_file = f"{tree_path}/{seg_name[i]}.treefile"
            tree2_file = f"{tree_path}/{seg_name[j]}.treefile"
            leaf_file = f"{feature_path}/{seg_name[i]}_{seg_name[j]}_leaf.csv"
            clade_file = f"{feature_path}/{seg_name[i]}_{seg_name[j]}_clade.csv"
            para_pair = SimpleNamespace(tree1=tree1_file,tree2=tree2_file,vec_size=1000,out_leaf=leaf_file,out_clade=clade_file)
            para_pair.min_clade_size = args.min_clade_size
            para_pair.max_clade_size = args.max_clade_size
            para_pair.max_clade_diff = args.max_clade_diff
            para_pair.jd = args.jd
            para_pair.pickle1 = f"{feature_path}/{seg_name[i]}.pickle"
            para_pair.pickle2 = f"{feature_path}/{seg_name[j]}.pickle"
            query_list.append(para_pair)
    with concurrent.futures.ProcessPoolExecutor(min(len(query_list),args.thread)) as executor:
        list(tqdm(executor.map(save_feature_pair, query_list), total=len(query_list)))


##################### single #####################################

""" basic: calling splits,clades """

def call_split(tree_file,ref_leaves,b_thres1=0.0003,b_thres2=0.001):
    # max_size=5000
    """ min_size: include it """
    t = dendropy.Tree.get(path=tree_file,schema='newick')
    fake_root = t.edges()[0]
    fake_root.length = 0
    leaves = np.array([i for i in t.leaf_nodes() if i.is_leaf])
    leaves_name = np.array([i.taxon.label for i in leaves])
    num = len(leaves_name)
    """ re-index """
    mapping = {leaves_name[i]:i for i in range(num)}
    mapping_index = [mapping[leaf] for leaf in ref_leaves]
    leaves = leaves[mapping_index]
    leaves_name = np.array([i.label for i in leaves])
    """ call internal nodes and edges """
    edges_map = t.bipartition_edge_map
    bipartitions = list(edges_map.keys())
    final_split = []
    final_length = []
    final_size = []
    all_split = []
    all_size = []
    all_string = []
    all_length = []
    internal_nodes = t.internal_nodes()
    """ analyze bipartitions """
    for node in internal_nodes:
        ed = node.edge
        split_string = node.bipartition.split_as_bitstring() # not reverse
        all_string.append(split_string)
        split = np.array([bool(int(i)) for i in split_string[::-1]])
        # split = split[mapping_index]
        n1 = np.count_nonzero(split)
        n2 = num - n1
        """ get the small one """
        if n1 <= n2:
            current_split = split[mapping_index]
            current_size = n1
        else:
            split = np.invert(split)
            current_split = split[mapping_index]
            current_size = n2
        all_split.append(current_split)
        all_size.append(current_size)
        all_length.append(ed.length)
        """ ignore the natural splits """
        if current_size < 2:
            continue
        else:
            """ edge's length >= b_thres2:retain """
            ed_length = np.round(ed.length,5)
            if np.round(ed_length,4) >= b_thres2:
                final_split.append(current_split)
                final_length.append(ed_length)
                final_size.append(current_size)
                continue
            elif ed_length < b_thres1: # b_thres1 < edge's length:ignore
                continue
            """ b_thres1 < edge's length < b_thres2: ignore some"""
            adj_eds = ed.adjacent_edges
            nn = []
            for a_ed in adj_eds:
                n_size = np.count_nonzero([int(i) for i in a_ed.leafset_as_bitstring()])
                nn.append(min(n_size,num-n_size))
            nn = np.array(nn)
            if sum(nn[:2]) < sum(nn[2:]):
                adj_eds = adj_eds[:2]
            else:
                adj_eds = adj_eds[2:]
            adj_eds_length = [np.round(a_ed.length,5) for a_ed in adj_eds if np.round(a_ed.length,5)>=0]
            flag = 0 # 0 ==> the edge's length < all its adjacent edges
            for a_ed_length in adj_eds_length:
                if ed_length > a_ed_length:
                    flag = 1
                    break
            if flag == 0:
                continue
            final_split.append(current_split)
            final_length.append(ed_length)
            final_size.append(current_size)
    all_info = [np.array(all_split),np.array(all_size),np.array(all_length),np.array(all_string)]
    return np.array(final_split),np.array(final_length),np.array(final_size),leaves_name,all_info

""" leaf: rank for strain count (ignore some edges with extremely small values) """

def split_to_rank(split,split_size):
    num = len(split[0])
    mat_num1 = np.zeros((num,num))
    mat_num2 = np.zeros((num,num))
    for i in range(len(split)):
        if split_size[i]<2:
            print(111)
            continue
        sp1 = split[i]
        temp1 = copy.copy(mat_num1[sp1])
        temp1[:,sp1] = temp1[:,sp1] + 1
        mat_num1[sp1,:] = temp1
        ### larger one
        sp2 = np.invert(split[i])
        temp2 = copy.copy(mat_num2[sp2])
        temp2[:,sp2] = temp2[:,sp2] + 1
        mat_num2[sp2,:] = temp2
    mat_num2 = mat_num2*0.9 / (mat_num2.max(axis=1, keepdims=True)+1e-6) # normalized, 0.9 is to make sure the rank prioritize small size
    mat_num2 = np.round(mat_num2,4)
    mat_num2 = np.where(mat_num1 != 0, mat_num1, mat_num2)
    # mat_num = np.where(mat_num1.max(axis=1, keepdims=True) < 2, mat_num2, mat_num1)
    mat_num = mat_num2
    mat_rank = []
    for vec in mat_num:
        mat_rank.append(rankdata(-vec,method='dense'))
    mat_rank = np.array(mat_rank)
    return mat_num,mat_rank

""" leaf: branch distance """

def call_branch_dist(file,ref_leaves):
    """ load the tree """
    t = dendropy.Tree.get(path=file,schema='newick')
    """ leaves name """
    leaves = t.taxon_namespace
    leaves_name = np.array([i.label for i in leaves])
    num_leaves = len(leaves_name)
    """ re-index """
    mapping = {leaves_name[i]:i for i in range(num_leaves)}
    temp = [leaves[mapping[leaf]] for leaf in ref_leaves]
    leaves = temp.copy()
    leaves_name = np.array([i.label for i in leaves])
    """ initialize the distance matrix """
    b_dist = np.zeros((num_leaves,num_leaves))
    pdm = t.phylogenetic_distance_matrix()
    """ find the edge length of leaves to their LA """
    edge_length = []
    for i in range(num_leaves):
        e_l = t.find_node_for_taxon(leaves[i])
        edge_length.append(e_l.edge_length)
    """ find the distance of leaves """
    for i in range(num_leaves-1):
        for j in range(i+1,num_leaves):
            d = pdm.patristic_distance(leaves[i],leaves[j])
            b_dist[i][j] = d-edge_length[i]-edge_length[j]
            b_dist[j][i] = d-edge_length[i]-edge_length[j]
    return np.round(b_dist,10)


""" leaf: strain counts """

def call_diff(mat_rank):
    """calling number of strains between two groups """
    posi = (mat_rank - mat_rank.T)>0 # strains above the target
    same_path = mat_rank.T<=1 # strains directly connected to the path
    mat_diff = np.zeros(mat_rank.shape, dtype=int)
    for i in range(mat_rank.shape[0]):
        for j in range(mat_rank.shape[1]):
            a1 = (mat_rank[i] < mat_rank[i, j]) & (posi[i]) & same_path[i] # not including the LCA one
            a2 = (mat_rank[j] < mat_rank[j, i]) & (posi[j]) & same_path[j] # not including the LCA one
            a3 = a1*a2
            mat_diff[i, j] = np.count_nonzero(a1) + np.count_nonzero(a2) - np.count_nonzero(a3)
    return mat_diff


""" leaf: order distance """
def split_to_dist(split,split_size,split_length):
    num = len(split[0])
    mat_dist = np.zeros((num,num))
    index = np.argsort(-split_size)
    for idx in index:
        n = split_size[idx]
        edge_len = split_length[idx]
        sp1 = split[idx]
        sp2 = np.invert(sp1)
        ratio1 = n*2/num
        ratio = ratio1**1 # double
        dist_offset = edge_len*ratio
        if edge_len > 0.01: #directly use the branch length
            dist_offset = edge_len
        dist = copy.copy(mat_dist[sp1])
        dist[:,sp2] = dist[:,sp2] + dist_offset
        mat_dist[sp1] = dist
    mat_dist = mat_dist + mat_dist.T
    mat_dist = np.log10(mat_dist+1e-100)
    for i in range(len(mat_dist)):# mask the strain itself
        mat_dist[i][i] = 100
    return mat_dist


""" basic: split sorting largest -> smallest """
def sort_split(all_info):
    split = all_info[0]
    split_size = all_info[1]
    split_length = all_info[2]
    split_string = all_info[3]
    idx = np.argsort(-split_size)
    split = split[idx]
    split_size = split_size[idx]
    split_length = split_length[idx]
    split_string = split_string[idx]
    sorted_all_info = [split,split_size,split_length,split_string]
    return sorted_all_info

""" basic: save the basic feature """

def feature_single(config):
    """ general parameters """
    tree1_file = config.tree1 # as the refernce file
    tree2_file = config.tree2
    out_file = config.output
    b_thres1 = config.b_thres1
    b_thres2 = config.b_thres2
    """ load data """
    ref_tree = dendropy.Tree.get(path= tree1_file, schema='newick')
    ref_leaves = np.array([i.label for i in ref_tree.taxon_namespace])
    num_leaves = len(ref_leaves)
    """ 1. basic: call split """
    split,length,split_size,name,all_info = call_split(tree2_file,ref_leaves,b_thres1,b_thres2)
    # all info ==> [split, size, length, string]
    all_info = sort_split(all_info)
    """ 2. basic: split rank """
    _, mat_rank = split_to_rank(split,split_size)
    """ 3. leaf: branch length distance """
    mat_length = call_branch_dist(tree2_file,ref_leaves)
    """ 4. mat: diff """
    mat_diff = call_diff(mat_rank)+1
    """ 5. leaf: order """
    mat_order = split_to_dist(all_info[0],all_info[1],all_info[2])
    """ save the variables """
    data = {'all_info':all_info,'mat_rank':mat_rank,'mat_length':mat_length,'ref_leaves':ref_leaves,'mat_order':mat_order,'mat_diff':mat_diff}
    f_pic = open(f"{out_file}",'wb')
    pickle.dump(data,f_pic,protocol=pickle.HIGHEST_PROTOCOL)


################### pair ##########################

""" basic: branch length coding """

def length_to_code(mat_length,config):
    thres = config.theta
    base = config.B
    count = mat_length/thres
    # mat_lcode = np.round(count*base)
    flag = count>1
    mat_lcode = np.round(count*flag*base)
    return mat_lcode

""" basic: generating the feature """

def extract_neighbors(config,ref_leaves,mat_scode,mat_lcode,mat_rank,mat_order1,mat_order2,idx,vec_size):
    index1 = np.argsort(mat_order1[idx])
    index2 = np.argsort(mat_order2[idx])
    mat_rank_sym = mat_rank
    # mat_rank_sym = np.maximum(mat_rank,mat_rank.T)
    r1 = mat_rank_sym[idx][index1][:vec_size]
    r2 = mat_rank_sym[idx][index2][:vec_size]
    leaf1 = ref_leaves[index1][:vec_size]
    leaf2 = ref_leaves[index2][:vec_size]
    # effi = rank_masking(r1,r2,min_rank,max_rank)
    vec1 = (mat_scode[idx][index1][:vec_size])
    vec2 = (mat_scode[idx][index2][:vec_size])
    vec1 = list(vec1) + list(mat_lcode[idx][index1][:vec_size])
    vec2 = list(vec2) + list(mat_lcode[idx][index2][:vec_size])
    return vec1,vec2,leaf1,leaf2,r1,r2


""" clade: match clade in two trees"""

def compare_clade(all_info1,all_info2,config):
    min_clade_size=config.min_clade_size
    max_clade_size=config.max_clade_size
    max_clade_diff=config.max_clade_diff
    jd_thres=config.jd
    """ load clade """
    split1 = all_info1[0]
    size1 = all_info1[1]
    string1 = all_info1[-1]
    split2 = all_info2[0]
    size2 = all_info2[1]
    string2 = all_info2[-1]
    """ filter small clade and large clade"""
    max_clade = min(max_clade_size,len(split1[0])/3)
    flag1 = (size1>=min_clade_size)*(size1<=max_clade)
    split1 = split1[flag1]
    size1 = size1[flag1]
    string1 = string1[flag1]
    flag2 = (size2>=min_clade_size)*(size2<=max_clade)
    split2 = split2[flag2]
    size2 = size2[flag2]
    string2 = string2[flag2]
    common_clade = []
    info_clade = []
    for i in range(len(size1)):
        sp1 = split1[i]
        si1 = size1[i]
        found_list = []
        found_jd = 0
        for j in range(len(size2)):
            sp2 = split2[j]
            si2 = size2[j]
            diff = si1-si2
            if diff>max_clade_diff:
                break
            elif abs(diff)>max_clade_diff:
                continue
            share = sp1*sp2
            union = sp1+sp2
            share_num = np.count_nonzero(share)
            union_num = np.count_nonzero(union)
            jd = np.round(share_num/union_num,3)
            if jd>jd_thres and jd>found_jd:#found_jd the best hit among j for i 
                found_jd = jd
                found_list.append([[share,union,sp1,sp2],[si1,si2,share_num,i,j,jd]])
        if found_jd > 0:
            common_clade.append(found_list[-1][0])
            info_clade.append(found_list[-1][1])
    info_clade = np.array(info_clade)
    """ repeat internal node filtering """
    selected = []
    if len(info_clade) == 0:
        return [],[],[]
    index = np.argsort(-info_clade[:,-1])
    share_all0 = []
    union_all0 = []
    sp1_all0 = []
    sp2_all0 = []
    final_info0 = []
    for idx in index:
        jth = info_clade[idx][-2]
        if jth not in selected:
            share_all0.append(common_clade[idx][0])
            union_all0.append(common_clade[idx][1])
            sp1_all0.append(common_clade[idx][2])
            sp2_all0.append(common_clade[idx][3])
            final_info0.append(info_clade[idx])
            selected.append(jth)
    final_info0 = np.array(final_info0)
    """ repeat clade filtering """
    share_all = []
    union_all = []
    sp1_all = []
    sp2_all = []
    final_info = []
    retain_clade = unique_clade(final_info0,share_all0,jd_thres)
    for i in retain_clade:
        final_info.append(final_info0[i])
        share_all.append(share_all0[i])
        union_all.append(union_all0[i])
        sp1_all.append(sp1_all0[i])
        sp2_all.append(sp2_all0[i])
    final_info = np.array(final_info)
    ix1 = [int(i) for i in final_info[:,-3]]
    ix2 = [int(i) for i in final_info[:,-2]]
    final_clade = [share_all,union_all,sp1_all,sp2_all]
    final_string = [string1[ix1],string2[ix2]]
    return final_clade,final_info,final_string



""" repeat clade filtering """
def unique_clade(final_info0,share_all0,thres=0.8):
    clade_size = final_info0[:,-4]
    size_idx = np.argsort(-clade_size) # size from largest to smallest
    retain_clade = [size_idx[0]] # the first one
    for i in range(1,len(clade_size)):
        site = size_idx[i]
        check_range = (clade_size[site]/clade_size[retain_clade])>thres
        flag = 0
        for j in np.array(retain_clade)[check_range]:
            if np.count_nonzero(share_all0[j]*share_all0[site])/np.count_nonzero(share_all0[j])>thres:
                flag = 1
                break
        if flag == 0:
            retain_clade.append(site)
    return retain_clade


""" clade: order distance"""

def split_to_dist_clade(all_info,clade,clade_size,union_clade):
    """ parameters
    split: all splits ==> all_info[0]
    split_size: the size ==> all_info[1]
    split_length: the edge length ==> all_info[2]
    clade: the bool array for clades
    clade_size: the size array
    union_clade: the union clade between two trees
    """
    split = all_info[0]
    split_size = all_info[1]
    split_length = all_info[2]
    num = len(split[0])
    mat_dist = np.zeros((len(clade),num))
    index = np.argsort(-split_size)
    split_length = split_length
    for idx in index:
        n = split_size[idx]
        edge_len = split_length[idx]
        sp1 = split[idx]
        sp2 = np.invert(sp1)
        """ clade index """
        large_idx = np.invert(np.array([i.any() for i in clade*sp1])) #clade in large group
        small_idx = np.invert(np.array([i.any() for i in clade*sp2])) # clade in small group
        ratio1 = (n-clade_size+1)*2/(num-clade_size+1) #clade in small group
        ratio2 = n*2/num # clade in large group
        dist_offset1 = edge_len*ratio1[small_idx] #clade in small group
        dist_offset2 = edge_len*ratio2 # clade in large group
        if edge_len > 0.01:
            dist_offset1 = np.array([edge_len]*len(dist_offset1))
            dist_offset2 = edge_len
        """ clade in small group """
        dist = copy.copy(mat_dist[small_idx])
        dist[:,sp2] = dist[:,sp2] + dist_offset1[:, np.newaxis]
        mat_dist[small_idx] = dist
        """ clade in large group """
        dist = copy.copy(mat_dist[large_idx])
        dist[:,sp1] = dist[:,sp1] + dist_offset2
        mat_dist[large_idx] = dist
    mat_dist = np.maximum(mat_dist,0)
    mat_dist = np.log10(mat_dist+1e-100)
    for i in range(len(mat_dist)):# mask the strain itself
        # print(mat_dist[i])
        mat_dist[i][union_clade[i]] = 100
    return mat_dist

""" clade: clade to leaf mapping """

def call_clade_mapping(file,clade_string):
    """ load the tree """
    t = dendropy.Tree.get(path=file,schema='newick')
    """ leaves name """
    leaves = [i for i in t.leaf_nodes() if i.is_leaf]
    leaves = np.array(leaves)
    leaves_name = np.array([i.taxon.label for i in leaves])
    num_leaves = len(leaves_name)
    """ clade name """
    _ = t.bipartition_edge_map
    internal_nodes = t.internal_nodes()
    mapping = {n.bipartition.split_as_bitstring():n for n in internal_nodes}
    clade = [mapping[s] for s in clade_string]
    num_clade = len(clade)
    """ clade mapping """
    clade_mapping = []
    clade_edge_length = []
    for i in range(num_clade):
        node = clade[i]
        clade_size = len(node.leaf_nodes())
        bit = [bool(int(b)) for b in clade_string[i][::-1]]
        bit_size = np.count_nonzero(bit)
        comp_size = num_leaves - bit_size
        final_size = bit_size
        if bit_size > comp_size:
            bit = np.invert(bit)
            final_size = comp_size
        cleaf = node.leaf_nodes()[0]
        if clade_size == final_size:
            minus_dist = cleaf.distance_from_root() - node.distance_from_root() + node.edge_length - cleaf.edge_length
        else:
            parent = node.parent_node
            bit = [bool(int(b)) for b in parent.bipartition.split_as_bitstring()[::-1]]
            cleaf = node.sister_nodes()[0].leaf_nodes()[0]
            minus_dist = cleaf.distance_from_root() - parent.distance_from_root() + node.edge_length - cleaf.edge_length
        clade_mapping.append([cleaf,minus_dist])
        clade_edge_length.append(node.edge_length)
    return clade_edge_length,clade_mapping

"""clade: branch distance """ 

def call_clade_branch_dist(leaf_dist,clade_bit,cleaf_index,clade_mapping):
    clade_dist = copy.copy(leaf_dist[cleaf_index])
    minus_value = np.array([i[1] for i in clade_mapping])
    clade_dist = clade_dist - minus_value[:, np.newaxis]
    for i in range(len(clade_bit)):
        clade_dist[i][clade_bit[i]] = 0
    return np.round(clade_dist,10)


"""clade: strain count """ 

def call_diff_clade(mat_rank,clade_bit,cleaf_index):
    num_clade = len(clade_bit)
    num_leaf = len(mat_rank)
    """calling number of strains between two groups """
    posi = (mat_rank - mat_rank.T)>0 # strains above the target (include those have deep clade)
    same_path = mat_rank.T<=1 #strains have the rank 1 with the target (directed connected to the path)
    mat_diff = np.zeros((num_clade,num_leaf), dtype=int)
    chara = np.array([i for i in range(num_leaf)])
    clade_rank = []
    for k in range(num_clade):
        i = cleaf_index[k]
        effective = np.invert(clade_bit[k])
        c_rank = (mat_rank[i] + mat_rank.T[i]-1)
        c_rank[clade_bit[k]] = 0
        # clade_rank.append(mat_rank)
        clade_rank.append(rankdata(c_rank,method='dense')-1)
        for j in chara[effective]:
            a1 = (mat_rank[i] < mat_rank[i, j]) & (posi[i]) & same_path[i] & effective
            a2 = (mat_rank[j] < mat_rank[j, i]) & (posi[j]) & same_path[j] & effective
            a3 = a1*a2
            mat_diff[k, j] = np.count_nonzero(a1) + np.count_nonzero(a2) - np.count_nonzero(a3)
    clade_rank = np.array(clade_rank)
    return mat_diff,clade_rank


""" clade: all variables """

def clade_all(config,ref_leaves,tree_file1,tree_file2,all_info1,all_info2,mat_length1,mat_length2,mat_rank1,mat_rank2):
    """ 1. call matched clades """
    common_clade,info_clade,clade_string = compare_clade(all_info1,all_info2,config) # filtering clades by size and similarity
    if len(info_clade) ==0:
        return [],[],[],[],[],[]
    """
    common_clade = [share,union,clade1,clade2]
    info: size1,size2,share_num, i-th,j-th,jarcard_similarity
    """
    union_clade = common_clade[1]
    match_clade1 = common_clade[2]
    match_clade_size1 = info_clade[:,1]
    clade_string1 = clade_string[0]
    match_clade2 = common_clade[3]
    match_clade_size2 = info_clade[:,2]
    clade_string2 = clade_string[1]
    """ 2. call the order distance """
    clade_dist1 = split_to_dist_clade(all_info1,match_clade1,match_clade_size1,union_clade)
    clade_dist2 = split_to_dist_clade(all_info2,match_clade2,match_clade_size2,union_clade)
    """ 3. call the clade to leaf mapping """
    clade_edge1,clade_mapping1 = call_clade_mapping(tree_file1,clade_string1)
    clade_edge2,clade_mapping2 = call_clade_mapping(tree_file2,clade_string2)
    leaf_mapping = {ref_leaves[i]:i for i in range(len(ref_leaves))}
    cleaf_index1 = [leaf_mapping[i[0].taxon.label] for i in clade_mapping1]
    cleaf_index2 = [leaf_mapping[i[0].taxon.label] for i in clade_mapping2]
    """ 4. call the branch length distance and coding """
    clade_length1 = call_clade_branch_dist(mat_length1,match_clade1,cleaf_index1,clade_mapping1)
    clade_length2 = call_clade_branch_dist(mat_length2,match_clade2,cleaf_index2,clade_mapping2)
    """ 5. call the strain count """
    clade_diff1,clade_rank1 = call_diff_clade(mat_rank1,match_clade1,cleaf_index1)
    clade_diff2,clade_rank2 = call_diff_clade(mat_rank2,match_clade2,cleaf_index2)
    """ output """
    clade_meta = [common_clade,info_clade,ref_leaves,clade_edge1,clade_edge2]
    clade_order = [clade_dist1,clade_dist2]
    clade_length = [clade_length1,clade_length2]
    clade_diff = [clade_diff1+1,clade_diff2+1]
    clade_rank = [clade_rank1,clade_rank2]
    return clade_meta,clade_order,clade_diff,clade_rank,clade_length

""" clade: info saveing """
def clade_info_saving(clade_meta):
    common_clade = clade_meta[0]
    info_clade = clade_meta[1]
    ref_leaves = clade_meta[2]
    num_clades = len(info_clade)
    name_code = [f"clade{i}" for i in range(num_clades)]
    meta_info = pd.DataFrame()
    share = []
    union = []
    c1 = []
    c2 = []
    for i in range(num_clades):
        share.append(",".join([f"{s}" for s in ref_leaves[common_clade[0][i]]]))
        # union.append(",".join([f"{s}" for s in ref_leaves[common_clade[1][i]]]))
        c1.append(",".join([f"{s}" for s in ref_leaves[common_clade[2][i]]]))
        c2.append(",".join([f"{s}" for s in ref_leaves[common_clade[3][i]]]))
    meta_info['name'] = name_code
    meta_info['share_num'] = info_clade[:,2]
    meta_info['tree1_num'] = info_clade[:,0]
    meta_info['tree2_num'] = info_clade[:,1]
    meta_info['Jaccard_index'] = info_clade[:,-1]
    meta_info['tree1_edge'] = clade_meta[3]
    meta_info['tree2_edge'] = clade_meta[4]
    meta_info['share_leaves'] = share
    meta_info['tree1_leaves'] = c1
    meta_info['tree2_leaves'] = c2
    return meta_info


""" basic: feature generation """

def feature_output(config,ref_leaves,mat_diff,mat_lcode,mat_rank,mat_order,target_idx):
    vec_size = config.vec_size
    # target_idx = [] for training
    s1 = []
    s2 = []
    s3 = []
    s4 = []
    na1 = []
    na2 = []
    na3 = []
    na4 = []
    r1 = []
    r2 = []
    r3 = []
    r4 = []
    for idx in target_idx:
        v1,v2,leaf1,leaf2,ra1,ra2 = extract_neighbors(config,ref_leaves,mat_diff[0],mat_lcode[0],mat_rank[0],mat_order[0],mat_order[1],idx,vec_size)
        v3,v4,leaf3,leaf4,ra3,ra4 = extract_neighbors(config,ref_leaves,mat_diff[1],mat_lcode[1],mat_rank[1],mat_order[0],mat_order[1],idx,vec_size)
        s1.append(",".join([str(i) for i in v1]))
        s2.append(",".join([str(i) for i in v2]))
        s3.append(",".join([str(i) for i in v3]))
        s4.append(",".join([str(i) for i in v4]))
        na1.append(",".join([str(i) for i in leaf1]))
        na2.append(",".join([str(i) for i in leaf2]))
        na3.append(",".join([str(i) for i in leaf3]))
        na4.append(",".join([str(i) for i in leaf4]))
        r1.append(",".join([str(i) for i in ra1]))
        r2.append(",".join([str(i) for i in ra2]))
        r3.append(",".join([str(i) for i in ra3]))
        r4.append(",".join([str(i) for i in ra4]))
    return s1,s2,s3,s4,na1,na2,na3,na4,r1,r2,r3,r4



def feature_pair(config):
    """ general parameters """
    #distant_thres = config.distant_thres # for distinguishing inter-subtype 0.01
    #base = config.base # base value for long length coding 5
    vec_size = config.vec_size # feature size
    min_clade_size = config.min_clade_size
    max_clade_size = config.max_clade_size
    max_clade_diff = config.max_clade_diff
    jd_thres = config.jd
    """ load tree data """
    tree1_file = config.tree1
    tree2_file = config.tree2
    """ load the pickle data """
    pickle_file1 = f"{config.pickle1}"
    pickle_file2 = f"{config.pickle2}"
    with open(pickle_file1, 'rb') as f1:
        data_pickle1 = pickle.load(f1)
    with open(pickle_file2, 'rb') as f2:
        data_pickle2 = pickle.load(f2)
    ref_leaves = data_pickle1['ref_leaves']
    num_leaves = len(ref_leaves)
    """ extract feature """
    """ 1. basic: call split """
    all_info1 = data_pickle1['all_info']
    all_info2 = data_pickle2['all_info']
    """ 2. basic: split rank """
    mat_rank1 = data_pickle1['mat_rank']
    mat_rank2 = data_pickle2['mat_rank']
    """ 3. leaf: branch length distance """
    mat_length1 = data_pickle1['mat_length']
    mat_length2 = data_pickle2['mat_length']
    """ 4. mat: diff """
    mat_diff1 = data_pickle1['mat_diff']
    mat_diff2 = data_pickle2['mat_diff']
    """ 5. leaf: order """
    mat_order1 = data_pickle1['mat_order']
    mat_order2 = data_pickle2['mat_order']
    """ 6. leaf: generate the feature vector """
    mat_length = [mat_length1,mat_length2]
    mat_diff = [mat_diff1,mat_diff2]
    mat_order = [mat_order1,mat_order2]
    mat_rank = [(mat_rank1+mat_rank1.T-1),(mat_rank2+mat_rank2.T-1)]
    target_idx = np.array([i for i in range(num_leaves)])
    """ leaf_feature = [s1,s2,s3,s4,na1,na2,na3,na4,r1,r2,r3,r4] """
    leaf_feature = feature_output(config,ref_leaves,mat_diff,mat_length,mat_rank,mat_order,target_idx)
    name = ref_leaves.copy()
    """ 7. clade: generate the feature vector """
    clade_meta,clade_order,clade_diff,clade_rank,clade_length = \
    clade_all(config,ref_leaves,tree1_file,tree2_file,all_info1,all_info2,mat_length1,mat_length2,mat_rank1,mat_rank2)
    if len(clade_meta) ==0:
        clade_meta_info = []
        clade_feature = []
    else:
        num_clades = len(clade_order[0])
        target_idx = np.array([i for i in range(num_clades)])
        clade_meta_info = clade_info_saving(clade_meta)
        # clade_feature = feature_output(config,ref_leaves,clade_diff,clade_lcode,clade_rank,clade_order,target_idx)
        clade_feature = feature_output(config,ref_leaves,clade_diff,clade_length,clade_rank,clade_order,target_idx)
    return leaf_feature,name,clade_feature,clade_meta_info



def save_feature_pair(config):
    """ loading parameters """
    file_leaf = config.out_leaf
    file_clade = config.out_clade
    """ data loading """
    leaf_feature,name,clade_feature,clade_meta_info = feature_pair(config)
    """ save the leaves feature """ 
    s1,s2,s3,s4,na1,na2,na3,na4,r1,r2,r3,r4 = leaf_feature
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
    data['name'] = name
    data.to_csv(file_leaf,index=False)
    """ save the clade feature """
    data_clade = pd.DataFrame()
    if len(clade_feature) ==0:
        s1,s2,s3,s4,na1,na2,na3,na4,r1,r2,r3,r4 = [[-1,-1]]*12
    else:
        s1,s2,s3,s4,na1,na2,na3,na4,r1,r2,r3,r4 = clade_feature
    data_clade['na1'] = na1
    data_clade['na2'] = na2
    data_clade['na3'] = na3
    data_clade['na4'] = na4
    data_clade['s1'] = s1
    data_clade['s2'] = s2
    data_clade['s3'] = s3
    data_clade['s4'] = s4
    data_clade['r1'] = r1
    data_clade['r2'] = r2
    data_clade['r3'] = r3
    data_clade['r4'] = r4
    for key in clade_meta_info:
        data_clade[key] = clade_meta_info[key]
    data_clade.to_csv(file_clade,index=False)


################### training function #######################

def extract_feature_train(para):
    """ parameters """
    path = para[0]
    config = para[1]
    times=config.times
    vec_size = config.vec_size
    b_thres1 = config.b_thres1
    b_thres2 = config.b_thres2
    """ load data """
    tree1_file = f"{path}/seg1.treefile"
    tree2_file = f"{path}/seg2.treefile"
    label_file = f"{path}/reassort_leaves.txt"
    ref_tree = dendropy.Tree.get(path= tree1_file, schema='newick')
    ref_leaves = np.array([i.label for i in ref_tree.taxon_namespace])
    num_leaves = len(ref_leaves)
    """ sampling target strains """
    f = open(label_file,'r')
    reassorted = [j for j in f.readline().split()]
    normal = np.array([k for k in ref_leaves if k not in reassorted])
    rng = np.random.default_rng()
    flag = rng.permutation(len(normal))
    target_leaves = reassorted + [str(s) for s in normal[flag][:times*len(reassorted)]]
    mapping = {ref_leaves[k]:k for k in range(num_leaves)}
    target_idx = [mapping[k] for k in target_leaves]
    y = []
    name = []
    for leaf in target_leaves:
        name.append(f"{path}_{leaf}")
        if leaf in reassorted:
            y.append(1)
        else:
            y.append(0)
    """ extract feature """
    """ 1. call split """
    split1,length1,split_size1,name1,all_info1 = call_split(tree1_file,ref_leaves,b_thres1,b_thres2)
    split2,length2,split_size2,name2,all_info2 = call_split(tree2_file,ref_leaves,b_thres1,b_thres2)
    # all info ==> [split, size, length]
    """ 2. basic: split rank """
    mat_num1, mat_rank1 = split_to_rank(split1,split_size1)
    mat_num2, mat_rank2 = split_to_rank(split2,split_size2)
    """ 3. basic: branch length distance """
    mat_length1 = call_branch_dist(tree1_file,ref_leaves)
    mat_length2 = call_branch_dist(tree2_file,ref_leaves)
    """ 4. leaf: strain counts """
    mat_diff1 = call_diff(mat_rank1)+1
    mat_diff2 = call_diff(mat_rank2)+1
    """ 5. leaf: order """
    mat_order1 = split_to_dist(all_info1[0],all_info1[1],all_info1[2])
    mat_order2 = split_to_dist(all_info2[0],all_info2[1],all_info2[2])
    """ 6 leaf: rank to path nodes """
    mat_rank1 = mat_rank1 + mat_rank1.T -1
    mat_rank2 = mat_rank2 + mat_rank2.T -1
    """ 7. generate the feature vector """
    s1 = []
    s2 = []
    s3 = []
    s4 = []
    na1 = []
    na2 = []
    na3 = []
    na4 = []
    r1 = []
    r2 = []
    r3 = []
    r4 = []
    for idx in target_idx:
        v1,v2,leaf1,leaf2,ra1,ra2 = extract_neighbors(config,ref_leaves,mat_diff1,mat_length1,mat_rank1,mat_order1,mat_order2,idx,vec_size)
        v3,v4,leaf3,leaf4,ra3,ra4 = extract_neighbors(config,ref_leaves,mat_diff2,mat_length2,mat_rank2,mat_order1,mat_order2,idx,vec_size)
        s1.append(",".join([str(i) for i in v1]))
        s2.append(",".join([str(i) for i in v2]))
        s3.append(",".join([str(i) for i in v3]))
        s4.append(",".join([str(i) for i in v4]))
        na1.append(",".join([str(i) for i in leaf1]))
        na2.append(",".join([str(i) for i in leaf2]))
        na3.append(",".join([str(i) for i in leaf3]))
        na4.append(",".join([str(i) for i in leaf4]))
        r1.append(",".join([str(i) for i in ra1]))
        r2.append(",".join([str(i) for i in ra2]))
        r3.append(",".join([str(i) for i in ra3]))
        r4.append(",".join([str(i) for i in ra4]))
    return s1,s2,s3,s4,na1,na2,na3,na4,r1,r2,r3,r4,y,name


def run_train(config):
    """ general parameters """
    main_folder = config.main_folder
    b_thres1 = config.b_thres1
    b_thres2 = config.b_thres2
    vec_size = config.vec_size
    times = config.times
    thread = config.thread

    # --- Get all subfolders ---
    subfolders = [
        os.path.join(main_folder, d)
        for d in os.listdir(main_folder)
        if os.path.isdir(os.path.join(main_folder, d))
    ]

    # --- Get all target folders ---
    target_paths = []
    for subfolder in subfolders:
        targets = [
            os.path.join(subfolder, t)
            for t in os.listdir(subfolder)
            if os.path.isdir(os.path.join(subfolder, t))
        ]
        target_paths.extend(targets)
    print(f"total samples:{len(target_paths)}")
    # --- Build query list ---
    query_list = [
        (path, config)
        for path in target_paths
    ]
    # print(query_list)
    # --- Extract the feature ---
    with concurrent.futures.ProcessPoolExecutor(thread) as executor:
            result = list(tqdm(executor.map(extract_feature_train, query_list), total=len(query_list)))
    s1 = []
    s2 = []
    s3 = []
    s4 = []
    na1 = []
    na2 = []
    na3 = []
    na4 = []
    label = []
    name = []
    r1 = []
    r2 = []
    r3 = []
    r4 = []
    for res in result:
        s1 = s1 + res[0]
        s2 = s2 + res[1]
        s3 = s3 + res[2]
        s4 = s4 + res[3]
        na1 = na1 + res[4]
        na2 = na2 + res[5]
        na3 = na3 + res[6]
        na4 = na4 + res[7]
        label = label + res[-2]
        name = name + res[-1]
        r1 = r1 + res[8]
        r2 = r2 + res[9]
        r3 = r3 + res[10]
        r4 = r4 + res[11]
    return s1,s2,s3,s4,label,name,na1,na2,na3,na4,r1,r2,r3,r4

        