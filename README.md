# VReassort: a tool to identify Virus Reassortment

### E-mail: dhcai2@cityu.edu.hk



## Installation
#### Main Dependencies:
* Conda
* Python >=3.10.18
* MAFFT==7.525
* TrimAl==1.5.0
* [FastTreeMP](https://morgannprice.github.io/fasttree/)==2.1.11
* Dendropy==5.0.1
* ETE3==3.1.3
* Pytorch

#### Supporting Platfroms
* Linux
* Mac 


#### Installation via Conda

[![Anaconda-Server Badge](https://anaconda.org/dhcai21/vreassort/badges/version.svg)](https://anaconda.org/dhcai21/vreassort)
#### 
```
conda create -n vreassort
conda activate vreassort
conda install dhcai21::vreassort
```


#### General installation
```
conda create -n vreassort python==3.10.18
conda activate vreassort
conda install -c bioconda mafft==7.525 trimal==1.5.0
conda install pytorch::pytorch
pip install numpy==1.26 pandas==2.3.1 ete3==3.1.3 dendropy==5.0.1 Bio==1.8.0 scipy==1.10.1 tqdm==4.67.1

git clone https://github.com/dhcai21/VReassort.git
cd VReassort
chmod +x src/*
```

## Usage
#### Example
set the threads
`export OMP_NUM_THREADS=8`

`python vreassort.py --input --input.csv --out result --thread 8`

if you install VReassort via Conda, you may use:

`vreassort --input --input.csv --out result --thread 8`

## Input
Please modify the input file `input_file/input_seq.csv` to indicate the segment names and the paths to the segment sequences. 

If you want to use built trees from other tools (e.g., IQ-Tree2) as input, please use the input file `input_file/input_tree.csv` and use `--type tree` in the command. 

`Example`

| segment | path                  |
|---------|-----------------------|
| HA      | test_data/seg4.fasta  |
| NA      | test_data/seg6.fasta  |
| ...     | ...                   |


> ⚠️ **Important:**  If the names of segment sequences across the segment FASTA files are inconsistent, VReassort will treat sequences in the same order across segments as belonging to the same strains.


```
**seg4.fasta**
>Strain1
ACGTATCTACGACGT 
>Strain2
ACGTACGTTCTACTC
...

**seg6.fasta**
>Strain1
ACGTTCTACACGT
>Strain2
ACGTACTCTACGT
...
```


## Output

```
result/
├── fasta/
├── feature/
├── identification/
├── name_mapping.csv
└── tree/
```


#### Name Mapping
To ensure consistency in sequence names across segments, we use pseudo-names (Code) to replace the original names of the segment sequences. The `name_mapping.csv` file provides the relevant information.
`Example: name_mapping.csv`

| Name | HA   | NA   | ...  |
|------|------|------|------|
| S1   | F0   | F0   | ...  |
| S2   | F1   | F1   | ...  |
| S3   | F2   | F2   | ...  |
| S4   | F3   | F3   | ...  |
| S5   | F4   | F4   | ...  |
| ...  | ...  | ...  | ...  |

#### Identification 
`Example: HA_NA_leaf.csv`
| Score     | S1         | S2         | Name  |
|-----------|------------|------------|-------|
| 0.000439  | 0.000524   | 0.000354   | S103  |
| 0.000613  | 0.000539   | 0.000686   | S112  |
| 0.000439  | 0.000524   | 0.000354   | S51   |
| 0.999999  | 0.999999   | 0.999999   | S82   |
| 0.999999  | 0.999999   | 0.999999   | S83   |
| 0.989963  | 0.979926   | 0.999999   | S98   |
| ...       | ...        | ...        | ...   |

`Example: HA_NA_clade.csv`
| Score      | S1         | S2         | Name   | Extra Info |
|------------|------------|------------|--------|------------|
| 0.00863296 | 0.00587875 | 0.01138717 | clade0 | ...        |
| 0.06610173 | 0.10411298 | 0.02809048 | clade1 | ...        |
| 0.02921185 | 0.05427545 | 0.00414825 | clade2 | ...        |
| 0.01904643 | 0.01755518 | 0.02053767 | clade3 | ...        |
| 0.00101161 | 0.00125027 | 0.00077295 | clade4 | ...        |
| 0.02334884 | 0.03315961 | 0.01353806 | clade5 | ...        |
| 0.00484115 | 0.00183183 | 0.00785047 | clade6 | ...        |
| 0.02662343 | 0.02405226 | 0.02919459 | clade7 | ...        |
| ...        | ...        | ...        | ...    | ...        |

Note: Only clades with sizes between 5 to 150 are reported.





## 📘&nbsp; License
VReassort is released under the terms of the [GPL v3.0 License](https://github.com/dhcai21/VReassort?tab=GPL-3.0-1-ov-file).
