# VReassort: A Tool to Identify Virus Reassortment

**E-mail:** dhcai2@cityu.edu.hk  
**Repository:** [https://github.com/dhcai21/VReassort](https://github.com/dhcai21/VReassort)

---

VReassort is a bioinformatics tool designed to **identify viral reassortment events** from genome sequences or constructed phylogenetic trees.  
It utilizes tree-based analysis to reveal potential reassortment events across multiple genomic segments by pair-wise comparison.  
The tool also supports data simulation and neural network model training workflows for advanced customization.

---

## Installation

### Main Dependencies

* Conda  
* Python ≥ 3.10.18  
* MAFFT == 7.525  
* TrimAl == 1.5.0  
* [FastTreeMP](https://morgannprice.github.io/fasttree/) == 2.1.11  
* DendroPy == 5.0.1  
* ETE3 == 3.1.3  
* PyTorch  

### Supporting Platforms

* **Linux**
* **macOS**

---

### Installation via Conda

[![Anaconda-Server Badge](https://anaconda.org/dhcai21/vreassort/badges/version.svg)](https://anaconda.org/dhcai21/vreassort)

```bash
conda create -n vreassort
conda activate vreassort
conda install dhcai21::vreassort
```

---

## Quick Start Example

### Example 1 – Sequence‑based Reassortment Detection

Set the number of threads:

```bash
export OMP_NUM_THREADS=8
```

Run VReassort using segment sequence FASTA files specified in a CSV file:

```bash
python vreassort.py --input input_file/input_seq.csv --out result --type seq --thread 8
```

If VReassort has been installed via Conda, you may also run:

```bash
vreassort --input input.csv --out result --type seq --thread 8
```

### Example 2 – Tree‑based Reassortment Detection

If you already have phylogenetic trees built with tools such as **FastTreeMP**, **IQ‑TREE 2**, or **RAxML**,  
you can use those trees directly instead of raw sequences:

```bash
python vreassort.py --input input_file/input_tree.csv --out result_tree --type tree --thread 8
```

---

## Input Format

Edit the CSV file `input_file/input_seq.csv` to indicate your segment names and corresponding file paths.

### Example: `input_seq.csv`

| segment | path                 |
|----------|----------------------|
| HA       | test_data/seg4.fasta |
| NA       | test_data/seg6.fasta |
| ...      | ...                  |

> ⚠️ **Important:** If the names of segment sequences differ across files,  
> VReassort will treat the sequences in the same line order as corresponding to the same strain.

### Example FASTA files

```
# seg4.fasta
>Strain1
ACGTATCTACGACGT
>Strain2
ACGTACGTTCTACTC
...

# seg6.fasta
>Strain1
ACGTTCTACACGT
>Strain2
ACGTACTCTACGT
...
```

---

## Output


All results are written to the specified output folder (e.g., `result_multi/`):

```
result_multi/
├── fasta/
├── feature/
├── identification/
├── name_mapping.csv
├── result_files.csv
└── tree/
```

**Description:**

- **`identification/`**  
  Contains the reassortment **scores for each strain and clade** derived from the model’s predictions.  
  Each file in this folder corresponds to a segment pair and lists the confidence or probability scores indicating potential reassortment events.

- **`name_mapping.csv`**  
  Records the mapping between the **original sequence names** and the **consistent internal names** used across all analysis outputs.  
  This ensures sequence correspondence among different segments and standardizes identifiers within the results.

- **`tree/`**
  Contains phylogenetic trees for all segments.


- **`result_files.csv`**  
  Aggregates the reassortment identification results for **all segment pairs**.  
  This file serves as a summary index for downstream analyses.  
  You can use it directly to cluster segments based on reassortment patterns:

  ```bash
  python src/seg_clustering.py -i result_files.csv
  ```

  Running this command will **cluster related segments** (see the “Segment Clustering” section below).

---

## Arguments

```
usage: vreassort.py [-h] --input INPUT --type TYPE [--out OUT]
                    [--model MODEL] [--thread THREAD] [--reshuffle]
                    [--n_thres N_THRES] [--b_thres1 B_THRES1]
                    [--b_thres2 B_THRES2] [--theta THETA] [--B B]
                    [--min_clade_size MIN_CLADE_SIZE]
                    [--max_clade_size MAX_CLADE_SIZE] [--jd JD]
                    [--max_clade_diff MAX_CLADE_DIFF]
```

| Option | Description |
|---------|-------------|
| `--input` | CSV file listing paths to segment sequences or trees. |
| `--type` | Input data type: `'seq'` (default) or `'tree'`. |
| `--out` | Output folder (default: `result`). |
| `--model` | Path to the trained model file (default: `src/trained.pt`). |
| `--thread` | Number of CPU threads. |
| `--reshuffle` | Randomly reshuffle sequence order to introduce topology variation. |
| `--n_thres` | Minimum normalization boundary (default: 20). |
| `--b_thres1` | Lower branch-length threshold (default: 0.0003). |
| `--b_thres2` | Upper branch-length threshold (default: 0.001). |
| `--theta` | Branch-length distance threshold for increment calculation (default: 0.01). |
| `--B` | Balance factor for the branch-length component (default: 5). |
| `--min_clade_size` | Minimum clade size (default: 5). |
| `--max_clade_size` | Maximum clade size (default: 150). |
| `--jd` | Minimum Jaccard index of matching clades between segments (default: 0.8). |
| `--max_clade_diff` | Maximum strain difference allowed between clades (default: 10). |

---
## Parameter Adjustment

VReassort provides several tunable parameters for fine‑grained control of the reassortment detection process.  
These parameters influence how branch lengths, normalization, and distance thresholds are interpreted during tree analysis.

| Parameter | Description | Default |
|------------|--------------|----------|
| `--n_thres N_THRES` | Minimum **normalization cap boundary**. Higher values increase robustness against noise in large, highly similar datasets by limiting the normalization range. | `20` |
| `--b_thres1 B_THRES1` | **Lower branch collapsing threshold.** Branches shorter than this value are collapsed. | `0.0003` |
| `--b_thres2 B_THRES2` | **Upper branch collapsing threshold.** Branches with lengths between `b_thres1` and `b_thres2` may be collapsed depending on balance and similarity configuration. | `0.001` |
| `--theta THETA` | **Branch length distance threshold** used when calculating incremental branch length changes between strains. | `0.01` |
| `--B B` | **Balance factor** that scales the contribution of branch length distance relative to node‑based components. Higher values make the model more sensitive to small changes in branch length. | `5` |

### Practical Guidelines

1. **Branch Collapsing (`b_thres1`, `b_thres2`)**  
   - Branches shorter than `0.0003` are always collapsed.  
   - Branches with lengths between `0.0003` and `0.001` may also be collapsed depending on the analysis mode.  
   - If you intend to detect reassortment among **highly similar strains** (e.g., with > 99.5 % sequence identity), consider **decreasing** these thresholds.  
   - Lower thresholds make the model more likely to preserve subtle short branches that could represent close parental lineages.  
   - ⚠️ Lower thresholds also have a **higher risk of producing false reassortment events** due to unstable short branches.

2. **Distance Increment Threshold (`theta`)**  
   - Controls how many “increments” in branch length are considered when comparing two strains.  
   - Smaller values make the model **more aggressive** in identifying potential reassortments among very similar strains (shorter evolutionary distances).  
   - However, excessively low values may amplify noise in unstable short branches, increasing false positives.

3. **Balance Factor (`B`)**  
   - Adjusts the weighting between **node‑based** structural signals and **branch length distance** contributions.  
   - Setting a **larger `B`** increases sensitivity to changes in branch length, useful when branch length distances are reliable.  
   - Because short branches are often less stable, combining a **small `theta`** value with a **high `B`** value can also **increase false positives** (spurious reassortments).

4. **Normalization Cap (`n_thres`)**  
   - Defines the lower bound (or cap) used during normalization.  
   - For large datasets (e.g., > 500 strains) with many highly similar sequences (> 99.5 % identity), increasing `n_thres` (e.g., 30 – 40) helps **reduce false reassortment calls** by dampening extreme normalization values.  
   - A **smaller `n_thres`** retains sensitivity to reassortments involving closely related parental strains, but may increase the number of false positives.

---

### Recommended Strategy

- **For typical datasets**  
  (moderate size, sequence divergence > 0.5 %):  
  Keep all default parameters.  
  The default configuration balances accuracy and computational efficiency for most use cases.

- **For large homogeneous datasets**  
  (hundreds of near‑identical strains):  
  Use the following adjusted setting to improve normalization and reduce false reassortments:

  ```bash
  --n_thres 40
  ```

- **For aggressively detecting reassortments with close parental strains**  
  Use the following settings to increase sensitivity when identifying reassortments among very closely related strains:

  ```bash
  --b_thres1 0.0001 --b_thres2 0.0005 --theta 0.005 --B 10
  ```

  ⚠️ **Caution:** This configuration greatly increases sensitivity but may produce **many false (spurious) reassortment events**, especially when branch lengths are extremely short or unstable.
  
  This configuration enhances the detection of reassortment among highly similar strains while minimizing noise introduced by small branch‑length fluctuations.

## Segment Clustering (for Multi‑Segment Comparison)

When analyzing **three or more segments**, cluster predicted reassortments by running:

```bash
python src/seg_clustering.py -i result_summary.csv -o cluster -tc 0.6 -tl 0.6 -j 0.8 -dc -dl
```

| Option | Description |
|---------|-------------|
| `-i` | Input CSV file containing reassortment results (`result_summary.csv`). |
| `-o` | Output prefix or path (default: `cluster`). |
| `-tl`, `-tc` | Score thresholds for leaf‑ and clade‑based clustering (default: 0.6). |
| `-j` | Minimum Jaccard index between clades (default: 0.8). |
| `-dl`, `-dc` | Flags to perform leaf‑ or clade‑level clustering. |

---

## Related Workflows

VReassort provides other components of the project:

- [**Data Simulation & Tree Generation**](https://github.com/dhcai21/VReassort/wiki/Simulation%E2%80%90and%E2%80%90Tree%E2%80%90Generation)
- [**Model Training Workflow**](https://github.com/dhcai21/VReassort/wiki/Model%E2%80%90Training)

These modules generate synthetic training data, extract structural features, and train the neural network model (`trained.pt`) used by VReassort for reassortment detection.

---

## License

VReassort is distributed under the **[GPL v3.0 License](https://github.com/dhcai21/VReassort?tab=GPL-3.0-1-ov-file)**.  
© 2025 dhcai21 — City University of Hong Kong

---
