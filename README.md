AuPairWise: A tool to assess biological replicability without using replicates
========

## Paper: 
Ballouz S, Gillis J (2016) AuPairWise: A Method to Estimate RNA-Seq Replicability through Co-expression. PLoS Comput Biol 12(4): e1004868. https://doi.org/10.1371/journal.pcbi.1004868

### Code 
- bin/run_on_sample_data.R : demonstration 
- bin/AuPairWise.R : a script to test the quality of your RNAseq experiment
- bin/helper.R : helper functions and libraries necessary for analyses

### Data and example run 
- sample/sample_brainspan.Rdata
- data/pairs.Rdata 

### Summary
 AuPairWise is a tool to estimate how replicable your expression experiment/study is without using explicit replicates. 
It uses a noise model to peturb a sample in your expression data, and then tests how this has affected the co-expression of gene pairs that are known to be coexpressed. The set of gene pairs are used as an indicator to determine the perturbed sample. The corresponding AUROC calculation is a measure of how "well" they do in this task. This is repeated over numerous runs and different levels of noise, and average AUROC is returned. The equivalent experiment is performed concurrently on a random set of pairs, in order to establish a null for your expression data. If the random pairs perform equivalently, there is more correlated noise in your experiment. This is reflected in the p-values that are calculated based on the distributions of the AUROCs of each run. 

![summary](https://github.com/sarbal/AuPairWise/blob/master/suppl/imgs/Fig9_new.png "Method summary")

###########################################################################################
### 1. Setting up the environment
###########################################################################################

You first need to label the directory you have saved AuPairWise in: eg. ``` masterdir="C:/AuPairWise/" ```
And label your experiment: eg. ``` name = "my_experiment" ```
All output will be saved in: ``` out = paste(masterdir,name,sep="/") ```

###########################################################################################
### 2. Expression data
###########################################################################################

#### a. To run AuPairWise, you need an expression dataset.
Currently, we have not implemented any pre-processing steps, so please make sure that the data is
set up as a matrix, with columns as your individual samples, and rows as genes.
The row names should be labelled by their gene entrez IDs.
The columns should be labelled by their sample IDs.
The expression dataset should be placed in the variable: ``` exprs ``` 

#### b. Gene pairs 
Another crucial element are the stoichiometric pairs, known to be coexpressed.
These pairs are located in the file ``` pairs.Rdata ```, labelled as ``` stoich.pairs ```.
They can be modified, but the variable must contain two columns, whereby each row
contains a pair of genes, labelled by their entrez gene IDs.

###########################################################################################
### 3. Running AuPairWise
###########################################################################################

Once the environment variables and the expression data is loaded, you can run the
script ``` run_AuPairWise_on_data.R ``` :
In an R session: ``` source("run_AuPairWise_on_data.R") ``` 
which checks for the necessary variables and runs the function: ``` summary = run_APW(exprs, out, stoich.pairs)``` 
The default repeats is set to 10, but we recommend 100 or 1000. However, this takes considerably 
more time. The process can be split up and parallelized but this is currently not implemented.


###########################################################################################
### 4. Output summary
###########################################################################################

#### a. Results:
The output of the ``` run_APW()``` function is a list that contains 6 elements. Here, we've named it summary.
The first element is a table with the average AUROCs for the noise factors that the model ran: ``` summary$aurocs ```
The first row of this table has the stoichiometric pair results, and the second row has the random pairs run.
The next two lists are the statistics for the above AUROCs (standard deviation and standard error): ``` summary$aurocs.sd ```  and ``` summary$aurocs.se```
The noise factors tested are listed in ```summary$n.factors```. These are the default values, and can 
be set in the ```run_APW()``` function, eg. ``` run_APW(exprs, out, stoich.pairs, n.factors=c(0,10,50) )```
The p-valaue of the Wilcoxon rank sum test between the stoichiometric pairs and the random pairs for each noise factor are in ```summary$pvals```
And finally, the estimated noise factors for any AUROC is given in: ``` summary$stats```

|AUROCs|0|1|2|5|10|15|20|25|50|100|
|---|---|---|---|---|---|---|---|---|---|---|
|Stoichiometric pairs|0.59|0.59|0.54|0.58|0.82|0.90|0.99|0.99|1.00|1.00|
|Random|0.59|0.49|0.55|0.59|0.74|0.79|0.88|0.95|1.00|1.00|


|P-values|0|1|2|5|10|15|20|25|50|100|
|---|---|---|---|---|---|---|---|---|---|---|
|Ranksum test|0.9096|0.2242|1.0000|0.9698|0.1381|0.0598|0.0025|0.0246|NaN|NaN|

#### b. Visualizing results:
You can view the summary results with the ``` plot_summary_results(summary)``` function.

![Sample output](https://github.com/sarbal/AuPairWise/blob/master/suppl/imgs/summary_encode.png "Sample output")
 
###########################################################################################
### 5. Intepreting results
###########################################################################################
A low noise factor implies that your experiment 'replicates' the coexpression patterns expected well
enough to be detected by small perturbations.
However, this should also be compared to the performances given by the random pairs.

###########################################################################################
### 6. Extras
###########################################################################################
To run on a different list of gene pairs, modify the stoich.pairs variable.
This a matrix of gene pairs, the first column has gene A and second column gene B.
These need to be entrez gene IDs, or must match the gene labels of your expresssion dataset.

You can also combine multiple runs, but this requires some careful labelling of your output files. 

