# Genomic encoding of transcriptional burst kinetics
Repository for transcriptional burst kinetics inference and analysis from allelic single-cell RNA-seq described in the article **Genomic encoding of transcriptional burst kinetics** by Anton J. M. Larsson, Per Johnsson, Michael Hagemann-Jensen, Leonard Hartmanis, Omid R. Faridani, Björn Reinius, Åsa Segerstolpe, Chloe M. Rivera, Bing Ren and Rickard Sandberg Nature. 2019 Jan 2. doi: [10.1038/s41586-018-0836-1](http://dx.doi.org/10.1038/s41586-018-0836-1). [Epub ahead of print].

We are currently uploading Jupyter notebooks that will contain the complete analysis.

## txburst scripts
### System Requirements for inference and hypothesis testing scripts

_txburstML.py_, _txburstPL.py_ and _txburstTEST.py_  are all python3 scripts with dependencies:

```
pandas: 0.19.2
numpy: 1.9.0
scipy: 1.0.0
joblib: 0.11
```
No further installation is needed.

### txburstML.py

_txburstML.py_ estimates parameters of bursting kinetics given transcript counts for an allele of a gene.

If you have estimated the allelic transcript counts from the fraction of allelic reads, it is important to consider what is missing data in this case. Genes with expression (reads) but no allelic reads are different from genes without expression (and therefore no allelic reads). We handle the first case as missing data, since it is not possible to assign the expression. In the second case, we replace the NaN with a 0 since there is genuinely no detected expression. Omitting this step may have severe effects on the quality of the inference.

#### Usage

usage: txburstML.py [-h] [--njobs NJOBS] file

Maximum likelihood inference of bursting kinetics from scRNA-seq data

**positional arguments:**
  file           .csv file with allelic-resolution UMI counts

**optional arguments:**
  -h, --help     show this help message and exit
  --njobs NJOBS  Number of jobs for the parallelization, default 50

#### Output 

_txburstML.py_ outputs a pickled pandas dataframe which for each gene contains an array with three entries [k_on, k_off, k_syn] where the burst frequency = k_on and burst size = k_syn/k_off, and a boolean indicating whether gene passed a rudimentary filtering step based on the quality of the inference procedure. For example:

|gene name | parameters |	keep |
| --- | --- | --- | 
|gene1 |	[k_on,k_off,k_syn]_gene1	| True	|
|gene2 |	[k_on,k_off,k_syn]_gene2	| False	|

#### Example 
```
./txburstML.py UMI_counts.csv
```

The output of this command is a pickled pandas dataframe _UMI_counts_ML.pkl_ which can be loaded into python with the command
```
pd.read_pickle('UMI_counts_ML.pkl')
```

### txburstPL.py

_txburstPL.py_ calculates the confidence intervals parameters of bursting kinetics given transcript counts for an allele of a gene. Requires _txburstML.py_ to be run on the dataset to provide input for _txburstPL.py_.

#### Usage 

usage: txburstPL.py [-h] [--file file] [--njobs NJOBS] [--alpha ALPHA]
                    [--MLFile MLFILE]

Confidence intervals on parameters of bursting kinetics from scRNA-seq data

**optional arguments:**
  -h, --help       show this help message and exit
  --file file      .csv file file with allelic-resolution transcript counts
  --njobs NJOBS    Number of jobs for the parallelization, default 50
  --alpha ALPHA    Alpha significance level (default 0.05)
  --MLFile MLFILE  Maximum Likelihood file (required)
  
#### Output

_txburstML.py_ outputs a pickled pandas dataframe which for each gene has three arrays with the point estimate from _txburstML.py_ and the confidence intervals for burst frequency and size. For example:

|gene name | point estimates |	bf_CI | bs_CI |
| --- | --- | --- | --- |
|gene1 |	[k_on,k_off,k_syn]_gene1	| [bf_point,bf_lower,bf_upper]_gene1	| [bs_point,bs_lower,bs_upper]_gene1	 |
|gene2 |	[k_on,k_off,k_syn]_gene2	| [bf_point,bf_lower,bf_upper]_gene2	| [bs_point,bs_lower,bs_upper]_gene2	 |
  
#### Example
  ```
  ./txburstPL.py --file UMI_counts.csv --MLFile UMI_counts_ML.pkl
  ```
The output of this command is a pickled dataframe _UMI_counts_PL.pkl_ which can be loaded into python with the command
```
pd.read_pickle('UMI_counts_PL.pkl')
```

### txburstTEST.py

_txburstTEST.py_ performs hypothesis testing for differential burst frequency and size between two experiments. Requires _txburstML.py_ to be run on the dataset to provide input for _txburstTEST.py_.

#### Usage 

usage: txburstTEST.py [-h] [--file1 file1] [--file2 file2] [--njobs NJOBS]
                      [--ML1 ML1] [--ML2 ML2]

Hypothesis testing of differences in bursting kinetics from scRNA-seq data

**optional arguments:**
  -h, --help     show this help message and exit
  --file1 file1  .csv file 1 with allelic-resolution transcript counts
  --file2 file2  .csv file 2 with allelic-resolution transcript counts
  --njobs NJOBS  Number of jobs for the parallelization, default 50
  --ML1 ML1      Maximum Likelihood file 1 (required)
  --ML2 ML2      Maximum Likelihood file 2 (required)
#### Output

_txburstTEST.py_ outputs a pickled pandas dataframe with the two point estimates from _txburstML.py_ and p-values for the significance tests for differential bursting kinetics. For example:

|gene name | point estimates 1 | point estimates 2	|bf_pvalue | bs_pvalue |
| --- | --- | --- | --- | --- |
|gene1 |	[k_on,k_off,k_syn]_gene1_sample1	| [k_on,k_off,k_syn]_gene1_sample2	| burst frequency pvalue gene1	 | burst size pvalue gene1 |
|gene2 |	[k_on,k_off,k_syn]_gene2_sample1	| [k_on,k_off,k_syn]_gene2_sample2	| burst frequency pvalue gene2	 | burst size pvalue gene2 |
  
#### Example
  ```
  ./txburstTEST.py --file1 UMI_counts.csv --file2 UMI_counts_2.csv --ML1 UMI_counts_ML.pkl --ML2 UMI_counts_2_ML.pkl
  ```
The output of this command is a pickled dataframe _UMI_counts_PL.pkl_ which can be loaded into python with the command
```
pd.read_pickle('UMI_counts_TEST.pkl')
```
