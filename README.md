AuPairWise: A tool to assess biological replicability without using replicates
========

# Code 
- bin/run_on_sample_data.R : demonstration 
- bin/AuPairWise.R : a script to test the quality of your RNAseq experiment
- bin/helper.R : helper functions and libraries necessary for analyses

# Data and example run 
- sample/sample_brainspan.Rdata
- data/pairs.Rdata 

# Summary
 AuPairWise is a tool to estimate how replicable your expression experiment/study is without using explicit replicates. 
It uses a noise model to peturb a sample in your expression data, and then tests how this has affected the co-expression of gene pairs that are known to be coexpressed. The set of gene pairs are used as an indicator to determine the peturbed sample. How "well" they do this corresponds to the AUROC calculated. This is repeated over numerous runs and different levels of noise. The average AUROC for each noise level is returned. The equivalent experiment is performed concurrently on a random set of pairs, in order to establish a null for your expression data. If the random pairs perform equivalently, there is more correlated noise in your experiment.

#################################################################################################
# 1. Setting up the environment
#################################################################################################

You first need to label the directory you have saved AuPairWise in: eg. masterdir="C:/AuPairWise/"
And label your experiment: eg. name = "my_experiment"
All output will be saved in: out = paste(masterdir,name,sep="/")

#################################################################################################
# 2. Expression data
#################################################################################################

## a. To run AuPairWise, you need an expression dataset.
Currently, we have not implemented any pre-processing steps, so please make sure that the data is
set up as a matrix, with columns as your individual samples, and rows as genes.
The row names should be labelled by their gene entrez IDs.
The columns should be labelled by their sample IDs.
The expression dataset should be placed in the variable: "exprs"

## b. Another crucial element are the stoichiometric pairs, known to be coexpressed.
These pairs are located in the file pairs.Rdata, labelled as stoich.pairs.
They can be modified, but the variable must contain two columns, whereby each row
contains a pair of genes, labelled by their entrez gene IDs.

#################################################################################################
# 3. Running AuPairWise
#################################################################################################

Once the environment variables and the expression data is loaded, you can run the
script run_AuPairWise_on_data.R:
In an R session: source("run_AuPairWise_on_data.R")
which checks for the necessary variables and runs the function: summary = run_APW(exprs, out, stoich.pairs)
The default repeats is set to 10, but we recommend 100 or 1000. However, this takes considerably 
more time. The process can be split up and parallelized but this is currently not implemented.


#################################################################################################
# 4. Output summary
#################################################################################################

## a. Results:
The output of the run_APW() function is a list that contains 6 elements. Here, we've named it summary.
The first element is a table with the average AUROCs for the noise factors that the model ran: summary$aurocs
The first row of this table has the stoichiometric pair results, and the second row has the random pairs run.
The next two lists are the statistics for the above AUROCs (standard deviation and standard error): summary$aurocs.sd and summary$aurocs.se
The noise factors tested are listed in summary$n.factors. These are the default values, and can 
be set in the run_APW() function, eg. run_APW(exprs, out, stoich.pairs, n.factors=c(0,10,50) )
The p-valaue of the Wilcoxon rank sum test between the stoichiometric pairs and the random pairs for each noise factor are in summary$pvals
And finally, the estimated noise factors for any AUROC is given in: summary$stats

## b. Visualizing results:
You can view the summary results with the plot_summary_results(summary) function.

#################################################################################################
# 5. Intepreting results
#################################################################################################
A low noise factor implies that your experiment 'replicates' the coexpression patterns expected well
enough to be detected by small perturbations.
However, this should also be compared to the performances given by the random pairs.

#################################################################################################
# 6. Extras
#################################################################################################
To run on a different list of gene pairs, modify the stoich.pairs variable.
This a 2D matrix, the first column has gene A and second column gene B.
These need to be entrez gene IDs, or must match the gene labels of your expresssion dataset.
eg. summary = run_APW(exprs, out, pairs )



