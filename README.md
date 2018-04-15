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

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

