#!/usr/bin/env python3
import os
import sys
import numpy as np
from ete3 import Tree
import argparse
import csv
import platform
import subprocess
import argparse
from pathlib import Path

####### tree function ##########

def new_tree(file_in, num_reassort, min_thres, max_thres, min_size, max_size):
    """Generate a reassorted tree from the input tree."""
    ref_tree = Tree(file_in)
    ref_tree.standardize()

    # Obtain candidate clades to be relocated
    internals = [i for i in ref_tree.get_descendants() if not i.is_leaf()]
    leaves = ref_tree.get_leaves()
    final_leaves = [i.name for i in leaves]

    # Candidate selection based on size thresholds
    if min_size > 1:
        candidate = internals
    elif max_size == 1:
        candidate = leaves
    else:
        candidate = internals + leaves

    # Filter clades based on size
    node_num = np.array([len(i.get_leaves()) for i in candidate])
    candidate = np.array(candidate)
    valid = (node_num >= min_size) & (node_num <= max_size)
    final_candidate = candidate[valid]
    num_candidate = len(final_candidate)
    repeats = 0
    while True:
        repeats = repeats + 1
        # Randomly select clades for reassortment
        candidate_index = np.arange(num_candidate)
        num_finished = 0
        while num_finished != num_reassort:
            np.random.shuffle(candidate_index)
            included_strains = []
            selected_index = []
            num_finished = 0
            for index in candidate_index:
                temp = [s.name for s in final_candidate[index].get_leaves()]
                if any(name in included_strains for name in temp):
                    continue
                selected_index.append(index)
                included_strains.extend(temp)
                num_finished += 1
                if num_finished == num_reassort:
                    break

        # Mask the selected clades and ancestors
        mask = []
        for index in selected_index:
            selected_strains = final_candidate[index]
            mask += [selected_strains] + selected_strains.get_descendants() + [selected_strains.get_ancestors()[0]]

        # Select target locations
        info = []
        selected_target = []
        for index in selected_index:
            selected_strains = final_candidate[index]
            target_candidate = np.array([i for i in ref_tree.get_descendants() if i not in mask])
            dist_b = np.array([ref_tree.get_distance(i, selected_strains) for i in target_candidate]) - selected_strains.dist
            dist = dist_b / max(dist_b)
            flag = (dist > min_thres) & (dist <= max_thres)
            target_candidate = target_candidate[flag]
            if len(target_candidate) == 0:
                # print("no_candidate")
                continue
            target_index = np.random.randint(0, len(target_candidate))
            target = target_candidate[target_index]
            temp_b = ref_tree.get_distance(selected_strains, target) - selected_strains.dist
            reassorted = [i.name for i in selected_strains.get_leaves()]
            reassorted_distance = temp_b / max(dist_b)
            info.append([
                reassorted,
                np.round(temp_b, 6),
                np.round(reassorted_distance, 6),
                len(reassorted)
            ])
            selected_target.append(target)
            mask.append(target)
        if len(selected_target) == len(selected_index):
            break
        if repeats == 100:
            print("no suitable cases with the setting after 100 runs:")
            print(f"num_reassort: {num_reassort}")
            print(f"min_size: {min_size}")
            print(f"max_size: {max_size}")
            print("Highly possible: distance range is too small, extend it")
            print(f"min_dist: {min_thres}")
            print(f"max_dist: {max_thres}")
            return [],[]
    # Relocate selected clades
    for i, index in enumerate(selected_index):
        selected_strains = final_candidate[index]
        target = selected_target[i]
        subtree = selected_strains.detach()
        parent = target.get_ancestors()[0]
        target.detach()
        length = target.dist
        anchor = parent.add_child(dist=length * 0.9)
        target.dist = length * 0.1
        subtree.dist = length * 0.1
        anchor.add_child(target)
        anchor.add_child(subtree)

    ref_tree.prune(final_leaves, preserve_branch_length=True)
    ref_tree.standardize()

    return ref_tree, info


def generation(file_in, file_out, num_reassort, min_thres, max_thres, min_size=1, max_size=5):
    """Generate and save reassorted trees."""
    res = new_tree(file_in, num_reassort, min_thres, max_thres, min_size, max_size)

    if len(res[1]) < num_reassort:
        print("Regenerating tree (insufficient reassortants)")
        return
        # res = new_tree(file_in, num_reassort, min_thres, max_thres, min_size, max_size)

    reassorted_tree, info = res

    # Save the reassorted tree
    reassorted_tree.write(format=0, outfile=f"{file_out}/reassort_tree.nwk")

    # Save labels and distances
    with open(f"{file_out}/reassort_leaves.txt", "w") as fl, \
         open(f"{file_out}/reassort_distance.txt", "w") as fd, \
         open(f"{file_out}/reassort_detail.txt", "w") as fi:
        for strains, branch_distance, norm_distance, num_strains in info:
            for j in strains:
                fl.write(f"{j}\t")
            fd.write(f"{branch_distance}\t{norm_distance}\n")
            fi.write(f"{strains}\t{branch_distance}\t{norm_distance}\t{num_strains}\n")


def run_generation(tree, num_reassort, repeats, start, min_size, max_size, min_dist, max_dist, output):
    """Execute generation for one configuration."""
    start_idx = start
    end_idx = start_idx + repeats - 1

    os.makedirs(output, exist_ok=True)

    for i in range(start_idx, end_idx + 1):
        sample_dir = f"{output}/sample{i}"
        os.makedirs(sample_dir, exist_ok=True)

        os.system(f"cp {tree} {sample_dir}/ref_tree.nwk")

        print(f"[INFO] Generating sample{i} ({output})")
        generation(
            file_in=f"{sample_dir}/ref_tree.nwk",
            file_out=sample_dir,
            num_reassort=num_reassort,
            min_thres=min_dist,
            max_thres=max_dist,
            min_size=min_size,
            max_size=max_size,
        )

        print(f"[DONE] Saved sample{i} → {sample_dir}")


####### Sequence simulation functions #########


def run_cmd(cmd, cwd=None):
    """Run a shell command and print it."""
    print(f"[CMD] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def get_fasttree_path():
    """Detect platform and find the correct FastTree binary."""
    # Path of the running Python file
    script_dir = Path(__file__).resolve().parent

    system = platform.system().lower()
    if "linux" in system:
        fasttree_file = script_dir / "FastTreeMP_linux"
    elif "darwin" in system:  # macOS appears as 'Darwin'
        fasttree_file = script_dir / "FastTreeMP_mac"
    else:
        raise OSError(f"Unsupported platform: {system}")

    if not fasttree_file.exists():
        raise FileNotFoundError(f"FastTree binary not found: {fasttree_file}")

    # Make sure FastTree is executable
    run_cmd(["chmod", "+x", str(fasttree_file)])
    return fasttree_file


def process_sample(sample_dir, model, threads, length, do_msa, do_tree, fasttree_path):
    """Simulate, optionally align, and optionally rebuild trees for one sample."""
    sample_dir = Path(sample_dir)
    ref_tree = sample_dir / "ref_tree.nwk"
    reassort_tree = sample_dir / "reassort_tree.nwk"

    if not ref_tree.exists() or not reassort_tree.exists():
        print(f"[SKIP] {sample_dir.name}: missing tree files.")
        return

    print(f"\n[INFO] Processing {sample_dir.name}")

    # Step 1: Simulation with iqtree2 --alisim
    ref_prefix = sample_dir / "seg1"
    reassort_prefix = sample_dir / "seg2"

    run_cmd([
        "iqtree2", "--alisim", str(ref_prefix),
        "-m", model,
        "-t", str(ref_tree),
        "-af", "fasta",
        "-nt", str(threads),
        "--length", str(length)
    ])

    run_cmd([
        "iqtree2", "--alisim", str(reassort_prefix),
        "-m", model,
        "-t", str(reassort_tree),
        "-af", "fasta",
        "-nt", str(threads),
        "--length", str(length)
    ])

    ref_fa = ref_prefix.with_suffix(".fa")
    reassort_fa = reassort_prefix.with_suffix(".fa")

    if do_msa:
        # Step 2: Alignment with MAFFT
        ref_msa = sample_dir / "seg1_msa.fa"
        reassort_msa = sample_dir / "seg2_msa.fa"

        print(f"[INFO] Running MAFFT for {sample_dir.name}")
        with open(ref_msa, "w") as fout:
            subprocess.run(
                ["mafft", "--quiet", "--thread", str(threads), str(ref_fa)],
                stdout=fout,
                check=True
            )
        with open(reassort_msa, "w") as fout:
            subprocess.run(
                ["mafft", "--quiet", "--thread", str(threads), str(reassort_fa)],
                stdout=fout,
                check=True
            )
    else:
        ref_msa = ref_fa
        reassort_msa = reassort_fa

    if do_tree:
        # Step 3: Tree building with FastTreeMP
        ref_treefile = sample_dir / "seg1.treefile"
        reassort_treefile = sample_dir / "seg2.treefile"
        os.system(f"export OMP_NUM_THREADS={threads}")

        print(f"[INFO] Building trees with {fasttree_path.name} for {sample_dir.name}")
        with open(ref_treefile, "w") as fout:
            subprocess.run(
                [str(fasttree_path), "-cat", "4", "-nt", "-gtr", "-gamma", str(ref_msa)],
                stdout=fout,
                check=True
            )
        with open(reassort_treefile, "w") as fout:
            subprocess.run(
                [str(fasttree_path), "-cat", "4", "-nt", "-gtr", "-gamma", str(reassort_msa)],
                stdout=fout,
                check=True
            )

    print(f"[DONE] {sample_dir.name}")



def parse_arguments():
    """Parse command-line arguments or CSV input, distinguishing tree vs. sequence parameter sets."""
    parser = argparse.ArgumentParser(
        description="Generate reassorted trees or simulated sequences based on provided parameters or CSV input."
    )

    # Shared/general arguments
    parser.add_argument(
        "-csv", "--csv",
        default='',
        help="Path to a CSV file listing multiple parameter sets (applies to both tree and sequence modes)."
    )
    parser.add_argument(
        "-threads", "--threads",
        type=int,
        default=4,
        help="Number of threads for sequence tools (default: 4)"
    )

    # -------------------------------
    # Tree-related arguments
    # -------------------------------
    tree_group = parser.add_argument_group(
        "Tree parameters", 
        "Parameters controlling reassortment and tree generation"
    )

    tree_group.add_argument(
        "-tree", "--tree",
        help="Path to the input reference tree (.nwk file)."
    )
    tree_group.add_argument(
        "-num_reassort", "--num_reassort",
        type=int,
        default=3,
        help="(default: 3) Number of reassortments per tree."
    )
    tree_group.add_argument(
        "-repeats", "--repeats",
        type=int,
        default=10,
        help="(default: 3) Number of random replicate trees to generate."
    )
    tree_group.add_argument(
        "-min_size", "--min_size",
        type=int,
        default=1,
        help="(default: 10) Minimum clade size to be relocated."
    )
    tree_group.add_argument(
        "-max_size", "--max_size",
        type=int,
        default=10,
        help="(default: 10) Maximum clade size to be relocated."
    )
    tree_group.add_argument(
        "-min_dist", "--min_dist",
        type=float,
        default=0.0,
        help="(0-1, default: 0) Minimum distance threshold for relocation."
    )
    tree_group.add_argument(
        "-max_dist", "--max_dist",
        type=float,
        default=1.0,
        help="(0-1, default 1) Maximum distance threshold for relocation."
    )
    tree_group.add_argument(
        "-output", "--output",
        help="Main output folder for generated trees."
    )
    tree_group.add_argument(
        "-start", "--start",
        type=int,
        default=1,
        help="(default: 1) Starting index for sample numbering."
    )

    # -------------------------------
    # Sequence-related arguments
    # -------------------------------
    seq_group = parser.add_argument_group(
        "Sequence parameters", 
        "Parameters controlling sequence simulation and processing"
    )

    seq_group.add_argument(
        "-model", "--model",
        default="GTR+F+I+G4",
        help="Substitution model (default: GTR+F+I+G4)"
    )
    seq_group.add_argument(
        "-length", "--length",
        type=int,
        default=1500,
        help="Sequence length to simulate (default: 1500)"
    )
    seq_group.add_argument(
        "-do_msa","--do_msa",
        action="store_true",
        help="Perform MAFFT alignment step"
    )
    seq_group.add_argument(
        "-do_tree","--do_tree",
        action="store_true",
        help="Perform FastTreeMP tree building step from simulated sequences"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    
    print(args)

    if args.csv:
        # ====== Multiple runs from CSV ======
        with open(args.csv, newline='') as f:
            # Changed: now standard CSV (comma-delimited)
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                output_path = Path(row["output"]).resolve()

                print(f"\n[INFO] Processing row for base output: {output_path}")
                run_generation(
                    tree=row["tree"],
                    num_reassort=int(row["num_reassort"]),
                    repeats=int(row["repeats"]),
                    start=int(row["start"]),
                    min_size=int(row["min_size"]),
                    max_size=int(row["max_size"]),
                    min_dist=float(row["min_dist"]),
                    max_dist=float(row["max_dist"]),
                    output=output_path
                )

                print(f"[INFO] Generating sequences for {output_path}")

                # Detect FastTree binary once per run
                fasttree_path = get_fasttree_path()

                # Determine sample folders for this row
                start_idx = int(row["start"])
                end_idx = start_idx + int(row["repeats"]) - 1
                for i in range(start_idx, end_idx + 1):
                    sample_dir = output_path / f"sample{i}"
                    if not sample_dir.exists():
                        print(f"[WARN] Missing {sample_dir}, skipping sequence processing.")
                        continue

                    process_sample(
                        sample_dir,
                        args.model,
                        args.threads,
                        args.length,
                        args.do_msa,
                        args.do_tree,
                        fasttree_path
                    )

    else:
        # ====== Single command-line run ======
        required = [
            args.tree
        ]
        print(required)
        if not all(required):
            print("[ERROR] Missing arguments. Use either --csv <file> or -tree <ref.treefile>.")
            sys.exit(1)

        output_path = Path(args.output).resolve()
        run_generation(
            tree=args.tree,
            num_reassort=args.num_reassort,
            repeats=args.repeats,
            start=args.start,
            min_size=args.min_size,
            max_size=args.max_size,
            min_dist=args.min_dist,
            max_dist=args.max_dist,
            output=output_path
        )


        fasttree_path = get_fasttree_path()

        start_idx = args.start
        end_idx = start_idx + args.repeats - 1
        for i in range(start_idx, end_idx + 1):
            sample_dir = output_path / f"sample{i}"
            if not sample_dir.exists():
                print(f"[WARN] Missing {sample_dir}, skipping sequence processing.")
                continue

            process_sample(
                sample_dir,
                args.model,
                args.threads,
                args.length,
                args.do_msa,
                args.do_tree,
                fasttree_path
            )

    print("\n✅ Sequence generation and processing completed.")
    sys.exit(0)