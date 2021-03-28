# ICDE21-HST
This repository stores the source code of the proposed algorithms to construct or update a Hierarchically Separated Tree (HST) in the following paper.
ICDE21-HST-Full.pdf is the full paper and ICDE21-HST-github.ppsx is our presentation slides.

[1] Yuxiang Zeng, Yongxin Tong, Lei Chen. "HST+: An Efficient Index for Embedding Arbitrary Metric Spaces." IEEE International Conference on Data Engineering (**ICDE**), 2021.

If you find this work helpful in your research, please consider citing our paper and the bibtex are listed below:
```  
@inproceedings{DBLP:conf/icde/ZengTC21,
  author    = {Yuxiang Zeng and Yongxin Tong and Lei Chen},
  title     = {HST+: An Efficient Index for Embedding Arbitrary Metric Spaces},
  booktitle   = {{ICDE}},
  year      = {2021},
}
```  

Usage of the algorithms
---------------

### Environment

gcc/g++ version: 7.4.0 

Python version: 2.7 (with Numpy)

OS: Linux

### Compile the algorithm

**cd algorithm/construct && make all**

BF: the algorithm proposed in STOC'03, ie, BASE in our experiments.

ODP: the algorithm proposed in our ICDE'21 paper, ie, HST+DP in our experiments.

OPT: the algorithm proposed in our ICDE'21 paper, ie, HST+DPO in our experiments.

**cd algorithm/update && make all**

reHSTO: the baseline used in our ICDE'21 paper, ie, reHST in our experiments.

HSF: the algorithm proposed in our ICDE'21 paper, ie, HST+HSF in our experiments.

**Notice**, the variable **MAX_DIM** in global.h should be modified based on the dimensions of datasets before compiling.

### Datasets

dataset/construct: real datasets and synthetic datasets used in our ICDE'21 paper. For example, checkinNYC_2 and checkinTKY_2 are Foursquare datasets, and gaiaChengdu_2 and gaiaHaikou_2 are Didi datasets.

dataset/update: real datasets and synthetic datasets used in our ICDE'21 paper. For example, checkinTKY_2 is the Foursquare dataset, gaiaHaikou_2 is the Didi dataset, and exp_2 is the synthetic dataset under the exponential distribution.

dataset/construct/&#42;.py: our data generators of construction experiments for synthetic datasets and scalability tests.

dataset/update/&#42;.py: our data generators of update experiments for synthetic datasets and real datasets.

### Run the algorithm by scripts

run-scripts/construct: For experiments of constructions, run&#42;.py is used to run the experiments; calc&#42;.py is used to calculate the experimental results; plot&#42;.py is used the plot the experimental figures.

run-scripts/update: For experiments of updates, run&#42;.py is used to run the experiments; calc&#42;.py is used to calculate the experimental results; plot&#42;.py is used the plot the experimental figures.

**Notice**, the variable **MAX_DIM** in global.h should be modified based on the dimensions of datasets before running any experiments.

Contact
------------
- Yuxiang Zeng: yzengal@connect.ust.hk


