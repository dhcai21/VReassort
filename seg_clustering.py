import os
import sys
import random
import argparse
from collections import Counter
import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities

def non_overlapping_cliques(graph,seg_list):
    # Step 1: Find all maximal cliques
    all_cliques = list(nx.find_cliques(graph))
    # Step 2: Sort by size (descending)
    clique_scores = [
        (sum(dict(graph.degree(c)).values()), c)
        for c in all_cliques
    ]
    ranked_cliques = sorted(clique_scores, reverse=True)
    sorted_cliques = [i[1] for i in ranked_cliques]
    # Step 3: Greedy selection of non-overlapping cliques
    selected_nodes = set(sorted_cliques[0]) # initialize the first node
    selected_cliques = [sorted_cliques[0]] # initialize the first selected_cliques
    selected_size = [len(sorted_cliques[0])] # recode the size of the first selected cliques
    while len(sorted_cliques)>0:
        # for j in range(len(sorted_cliques)):
        clique = sorted_cliques[0]
        # if all(node in selected_nodes for node in clique):
        #     continue
        num = len(selected_cliques)
        flag = 0
        for i in range(len(selected_cliques)): # compare with all cliques in the selected list
            ref = selected_cliques[i]
            share = set(clique).intersection(set(ref))
            ### ---------------- no share ------------------ ###
            if len(share) == 0:
                continue
            ### ---------------- subset ------------------ ###
            elif len(share) == min(len(ref),len(clique)):
                flag = 1
                ref = list(set(ref).union(clique))
                selected_cliques[i] = ref
                selected_nodes.update(clique)
                selected_size[i] = max(selected_size[i],len(clique))
                sorted_cliques.remove(clique)
                break
            ### ---------------- 50% overlap ------------------ ###
            elif len(share) >= 0.5 * max(selected_size[i],len(clique)):
            # elif len(share) >= 0.5 * max(len(ref),len(clique)):
                flag = 1
                ref = list(set(ref).union(clique))
                selected_cliques[i] = ref
                selected_nodes.update(clique)
                selected_size[i] = max(selected_size[i],len(clique))
                sorted_cliques.remove(clique)
                break
            else:
            ### -------------- cutoff -----------------------###
                flag = 1
                diff = list(set(clique).difference(set(ref)))
                sorted_cliques.remove(clique)
                sorted_cliques.append(diff)
                clique_scores = [
                    (sum(dict(graph.degree(c)).values()), c)
                    for c in sorted_cliques
                ]
                ranked_cliques = sorted(clique_scores, reverse=True)
                sorted_cliques = [i[1] for i in ranked_cliques]
                break
        if flag == 0:
            selected_cliques.append(clique)
            selected_nodes.update(clique)
            selected_size.append(len(clique))
            sorted_cliques.remove(clique)
    orignal_edge = len(graph.edges())
    optimal_edge = 0
    for i in range(len(seg_list)-1):
        for j in range(i+1,len(seg_list)):
            optimal_edge = optimal_edge + 1
    # optimal_edge = len(seg_list)
    for i in range(len(selected_cliques)-1):
        for j in range(i+1,len(selected_cliques)):
            optimal_edge = optimal_edge - len(selected_cliques[i])*len(selected_cliques[j])
    return selected_cliques,orignal_edge,optimal_edge



def build_graph(seg_label,seg_list):
    G = nx.complete_graph(seg_list)
    G.remove_edges_from(seg_label)
    edge_num = len(G.edges())
    return G,edge_num


###################  pattern from the graph ###########################


def summarize_leaf_predictions(df_input, args):
    all_results = {}  # Collect predictions per strain name
    for idx, row in df_input.iterrows():
        seg1 = row['segi']
        seg2 = row['segj']
        path = row['leaf']
        if not os.path.exists(path):
            print(f"⚠️ File not found: {path}")
            continue
        try:
            df_leaf = pd.read_csv(path)
        except Exception as e:
            print(f"❌ Error reading {path}: {e}")
            continue
        # Ensure expected columns
        required_cols = {'name', 'score'}
        if not required_cols.issubset(df_leaf.columns):
            print(f"⚠️ Skipping file {path}, missing required columns {required_cols}")
            continue
        # Filter by score
        df_selected = df_leaf[df_leaf['score'] >= args.thres_leaf]
        # Label each selected strain
        label = f"({seg1},{seg2})"
        for strain in df_selected['name']:
            all_results.setdefault(strain.split()[0], []).append(label)
            # all_results.setdefault(strain, []).append(label)
    # Combine results
    summary_records = [
        {"name": strain, "labels": ";".join(labels)}
        for strain, labels in all_results.items()
    ]
    df_summary = pd.DataFrame(summary_records)
    # Collect unique segments from both segi and segj columns
    segi_unique = set(df_input["segi"].dropna().astype(str))
    segj_unique = set(df_input["segj"].dropna().astype(str))
    unique_segments = sorted(segi_unique.union(segj_unique), key=lambda x: int(x) if x.isdigit() else x)
    return df_summary, unique_segments

def summarize_leaf_predictions(df_input, args):
    all_results = {}  # Collect predictions per strain name
    for idx, row in df_input.iterrows():
        seg1 = row['segi']
        seg2 = row['segj']
        path = row['leaf']
        if not os.path.exists(path):
            print(f"⚠️ File not found: {path}")
            continue
        try:
            df_leaf = pd.read_csv(path)
        except Exception as e:
            print(f"❌ Error reading {path}: {e}")
            continue
        # Ensure expected columns
        required_cols = {'name', 'score'}
        if not required_cols.issubset(df_leaf.columns):
            print(f"⚠️ Skipping file {path}, missing required columns {required_cols}")
            continue
        # Filter by score
        df_selected = df_leaf[df_leaf['score'] >= args.thres_leaf]
        # Label each selected strain
        label = f"({seg1},{seg2})"
        for strain in df_selected['name']:
            all_results.setdefault(strain.split()[0], []).append(label)
            # all_results.setdefault(strain, []).append(label)
    # Combine results
    summary_records = [
        {"name": strain, "labels": ";".join(labels)}
        for strain, labels in all_results.items()
    ]
    df_summary = pd.DataFrame(summary_records)
    # Collect unique segments from both segi and segj columns
    segi_unique = set(df_input["segi"].dropna().astype(str))
    segj_unique = set(df_input["segj"].dropna().astype(str))
    unique_segments = sorted(segi_unique.union(segj_unique), key=lambda x: int(x) if x.isdigit() else x)
    return df_summary, unique_segments




def jaccard_index(set1, set2):
    intersection = set1 & set2
    union = set1 | set2
    return len(intersection) / len(union) if union else 0

def summarize_clade_predictions(df_input, args):
    all_clades = []
    for idx, row in df_input.iterrows():
        seg1 = row['segi']
        seg2 = row['segj']
        path = row['clade']
        if not os.path.exists(path):
            print(f"⚠️ File not found: {path}")
            continue
        try:
            df_clade = pd.read_csv(path)
        except Exception as e:
            print(f"❌ Error reading {path}: {e}")
            continue
        required_cols = {'share_leaves', 'score'}
        if not required_cols.issubset(df_clade.columns):
            print(f"⚠️ Skipping file {path}, missing required columns {required_cols}")
            continue
        df_selected = df_clade[df_clade['score'] >= args.thres_clade]
        label = f"({seg1},{seg2})"
        for leaves_str in df_selected['share_leaves']:
            current_leaves = set(x.strip() for x in str(leaves_str).split(',') if x.strip())
            if not current_leaves:
                continue
            matched = False
            for existing_clade in all_clades:
                if not existing_clade['all_sets']:
                    continue
                jaccs = [jaccard_index(current_leaves, s) for s in existing_clade['all_sets']]
                mean_jacc = sum(jaccs) / len(jaccs)
                if mean_jacc >= args.jaccard:
                    existing_clade['labels'].append(label)
                    existing_clade['all_sets'].append(current_leaves)
                    existing_clade['jacc_history'].append(mean_jacc)
                    matched = True
                    break
            if not matched:
                all_clades.append({'all_sets': [current_leaves], 'labels': [label], 'jacc_history': []})
    summary_records = []
    for clade in all_clades:
        intersection_leaves = set.intersection(*clade['all_sets']) if clade['all_sets'] else set()
        name = ",".join(sorted(intersection_leaves))
        labels = ";".join(clade['labels'])
        intersection_size = len(intersection_leaves)
        mean_jaccard = sum(clade['jacc_history']) / len(clade['jacc_history']) if clade['jacc_history'] else 1.0
        summary_records.append({"name": name, "labels": labels, "intersection_size": intersection_size, "mean_jaccard": mean_jaccard})
    df_summary = pd.DataFrame(summary_records)
    segi_unique = set(df_input["segi"].dropna().astype(str))
    segj_unique = set(df_input["segj"].dropna().astype(str))
    unique_segments = sorted(segi_unique.union(segj_unique), key=lambda x: int(x) if x.isdigit() else x)
    return df_summary, unique_segments



def parse_label_string(label_str):
    """
    Convert a label string like:
      '(seg1,seg2);(seg1,seg6);(seg6,seg8)'
    into a list of tuples:
      [('seg1', 'seg2'), ('seg1', 'seg6'), ('seg6', 'seg8')]
    """
    if not isinstance(label_str, str) or not label_str.strip():
        return []
    # Split by semicolon, strip parentheses/spaces, then split by comma
    pairs = []
    for item in label_str.split(';'):
        item = item.strip().strip('()')
        if not item:
            continue
        segs = [s.strip() for s in item.split(',')]
        if len(segs) == 2:
            pairs.append((segs[0], segs[1]))
    return pairs

def call_para():
    # --- Define Argument Parser ---
    parser = argparse.ArgumentParser(
        description="Segment clustering based on the predicted reassortments."
    )
    # --- Required arguments ---
    parser.add_argument(
        "-i", "--input",
        required=True,
        # default='result_summary.csv',
        help="Input CSV file containing the result from VReassort (example:'result_summary.csv'). "
    )
    parser.add_argument(
        "-o", "--out",
        # required=True,
        default='cluster',
        help="Output prefix or path and prefix (default: cluster)."
    )
    # --- Optional arguments ---
    parser.add_argument(
        "-tl", "--thres_leaf",
        type=float,
        default=0.6,
        help="Threshold value for leaf analysis (default: 0.6)."
    )
    parser.add_argument(
        "-tc", "--thres_clade",
        type=float,
        default=0.6,
        help="Threshold value for clade analysis (default: 0.6)."
    )
    parser.add_argument(
        "-j", "--jaccard",
        type=float,
        default=0.8,
        help="The least Jaccard Index between clades (default: 0.8)."
    )
    # --- Analysis flags ---
    parser.add_argument(
        "-dl", "--do_leaf",
        action="store_true",
        default=True,
        help="Perform leaf-level clustering if this flag is set."
    )
    parser.add_argument(
        "-dc", "--do_clade",
        action="store_true",
        help="Perform clade-level clustering if this flag is set."
    )
    # --- Parse arguments ---
    args = parser.parse_args()
    return args  # ✅ return the parsed arguments object


if __name__ == "__main__":
    # --- Get command-line parameters and the results ---
    args = call_para()
    files = args.input
    df_input = pd.read_csv(files)


    if args.do_leaf:
        df_leaf,seg_list = summarize_leaf_predictions(df_input, args)
        name = df_leaf['name'].values
        label = df_leaf['labels'].values
        num = len(name)
        clusters = []
        clu_num = []
        for i in range(num):
            seg_label = parse_label_string(label[i])
            graph,num_edge = build_graph(seg_label,seg_list)
            clique,orignal_num,revised_num = non_overlapping_cliques(graph,seg_list)
            clusters.append("\t||\t".join([",".join(c) for c in clique]))
            clu_num.append(len(clique))
        #--------------- save the clustering results ----------------------#
        df_out = pd.DataFrame()
        df_out['name'] = name
        df_out['num'] = clu_num
        df_out['cluster'] = clusters
        df_out.to_csv(f"{args.out}_leaf.csv",index=False)


    if args.do_clade:
        df_clade, seg_list = summarize_clade_predictions(df_input, args)
        if len(df_clade) == 0:
            print("No found reassortant clade pattern")
            exit()
        name = df_clade['name'].values
        label = df_clade['labels'].values
        inter_size = df_clade['intersection_size'].values
        mean_jacc = df_clade['mean_jaccard'].values
        num = len(name)
        clusters = []
        clu_num = []
        for i in range(num):
            seg_label = parse_label_string(label[i])
            graph, num_edge = build_graph(seg_label, seg_list)
            clique, original_num, revised_num = non_overlapping_cliques(graph, seg_list)
            clusters.append("\t||\t".join([",".join(c) for c in clique]))
            clu_num.append(len(clique))
        df_out = pd.DataFrame()
        df_out['name'] = name
        df_out['intersection_size'] = inter_size
        df_out['mean_jaccard'] = mean_jacc
        df_out['num'] = clu_num
        df_out['cluster'] = clusters
        df_out.to_csv(f"{args.out}_clade.csv", index=False)


    if not args.do_leaf and not args.do_clade:
        print("⚠️ No analysis selected. Use --do_leaf and/or --do_clade to perform analyses.")





