# VirReassort: a tool to identify Virus Reassortment

### Citation
Dehan Cai, Yanni Sun, Reconstructing viral haplotypes using long reads, Bioinformatics, Volume 38, Issue 8, 15 April 2022, Pages 2127–2134, https://doi.org/10.1093/bioinformatics/btac089

### E-mail: dhcai2-c@my.cityu.edu.hk
### Version: V3 (2022-08-02 updated)
##### Use [Medaka](https://github.com/nanoporetech/medaka) to polish the final result
##### Support datasets with large sizes (>100k reads) by applying MCL to subgraphs of reads. (use **-sg #_of_subgraphs**)

#### Version: V2 (2022-05-22 updated)
##### Use a C package of MCL.
##### Support multiprocessing.

## Installation
### Dependencies:
* Conda
* Python >=3.8
* samtools >= 1.4.1
* pysam
* [Medaka](https://github.com/nanoporetech/medaka)
* MCL
* Required python package: pandas >= 1.1.3 tqdm, scipy

### An easiler way to install
After cloning this respository, you can use anaconda to install the **rvhaplo.yaml** (Linux). This will install all packages you need. The command is: `conda env create -f rvhaplo.yaml -n rvhaplo`

#### An optional way to install
```
conda create -n rvhaplo python==3.8

conda activate rvhaplo

conda install -c bioconda samtools

conda install mcl

pip install medaka scipy pandas tqdm pysam
```
