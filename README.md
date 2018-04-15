# analysis of scRNAseq data of fetal coronary progenitors


Cluster cells using iterative rounds of RPCA to separate and define populations. 
Then determine whether each pair of subtypes is continuous or discrete using a statistical test on each pair of subtypes.

## Getting Started

Clone this repository into a new directory.
```
git clone https://github.com/gmstanle/coronary-progenitor-scRNAseq.git
cd coronary-progenitor-scRNAseq
```
Download two data objects from figshare into the data/ directory:
```
cd data
wget https://ndownloader.figshare.com/files/11086496
wget https://ndownloader.figshare.com/files/11089799
cd ..
```

### Installing


Open RStudio and install the requisite packages:


```
packages = read.table('packages.txt')
install.packages(packages)
```


## Running iRPCA

Open e12_iRPCA/e12_iRPCA.Rmd in RStudio

Run all chunks to group cells based on iterative rPCA.

To re-do the pipeline fully manually, set recalculate.pca <- TRUE. This may require changing some of the cutoffs in each step, but the file e12_iRPCA.pdf can be used as a guide for clustering.



## Authors

* **Geoffrey Stanley** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

