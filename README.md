# VReassort: a tool to identify Virus Reassortment

### E-mail: dhcai2@cityu.edu.hk

## Installation
### Main Dependencies:
* Conda
* Python >=3.10.18
* MAFFT==7.525
* TrimAl==1.5.0
* FastTree==2.1.11
* Dendropy==5.0.1
* ETE3==3.1.3
* Pytorch>=2.2.2


#### Installation via Conda

#### 
```
conda install 


```


#### General installation
```
conda create -n vreassort python==3.10.18
conda activate vreassort
conda install -c bioconda mafft==7.525 trimal==1.5.0
pip install torch==2.2.2 numpy==1.24.4 pandas ete3==3.1.3 dendropy==5.0.1 Bio scipy==1.10.1 tqdm

git clone https://github.com/dhcai21/VReassort.git
cd VReassort
chmod +x src/*

```



## 📘&nbsp; License
VReassort is released under the terms of the [GPL v3.0 License](https://github.com/dhcai21/VReassort?tab=GPL-3.0-1-ov-file).
