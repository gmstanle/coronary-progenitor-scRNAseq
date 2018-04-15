# analysis of scRNAseq data of fetal coronary progenitors


Cluster cells using iterative rounds of RPCA to separate and define populations. 
Then determine whether each pair of subtypes is continuous or discrete using a statistical test on each pair of subtypes.

## Getting Started

Clone this repository into a new directory.
```
git clone https://github.com/gmstanle/coronary-progenitor-scRNAseq.git
cd coronary-progenitor-scRNAseq
```
Download two data objects from figshare into the data/ directory


Linux:
```
cd data
wget https://ndownloader.figshare.com/files/11086496
wget https://ndownloader.figshare.com/files/11089799
cd ..
```

MacOS: use a web browser to open the above links, and they will begin downloading. Then move the data objects into the data folder:

```
cd data
mv ~/Downloads/fetal_venule_e12_seurat.Rdata .
mv ~/Downloads/fetal_venule_exp1_Nextseq.dt_goodCells.RData .
cd ..
```
Mac users who have not done so should download the XQuartz package via the XQuartz website.
### Installing


Start R and install the requisite packages:


```
pkg=as.character(read.table('packages.txt')[,1])
install.packages(pkg)
```


## Running iRPCA

Open e12_iRPCA/e12_iRPCA.Rmd in RStudio

Run all chunks to group cells based on iterative rPCA. Each iteration of rPCA has a folder with the results from rPCA on that group of cells. The first iteration is in the folder all_r1 and starts at the second chunk of the R notebook e12_iRPCA.Rmd. rPCA has been pre-computed for each step and the PC scores are stored in a file "PC_allscores.csv" in each iteration's folder (e.g. "all_r1/PC_allscores.csv"). These files contain the default PC scores (columns "PC1", "PC2", ...) as well as the sum-of-60 scores ("PC1.score", "PC2.score", ...). The pipeline uses only the modified sum-of-60 scores. The cutoffs used for iterating and defining cell types are hard-coded into each

To re-do the pipeline fully manually, set recalculate.pca <- TRUE. This may require changing some of the cutoffs in each step, and the file e12_iRPCA.pdf can be used as a guide for clustering. 

Cell cluster are added to the metadata file (data/e12_subtypeInfo.Rdata) in the second-to-last chunk. In the last chunk, the average expression of *in situ* validated marker genes is plotted across all clusters.


## Running PDT (pairwise discreteness test)

Open pdt/e12_pairwiseDiscretenessTest.Rmd.

Run all chunks - the first creates a Seurat object and adds the iRPCA-defined subtypes, the third chunk runs the pairwise discreteness test, and the fourth chunk plots the results of PDT.

## Authors

* **Geoffrey Stanley** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

