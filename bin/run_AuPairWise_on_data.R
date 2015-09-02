masterdir = ""

source(paste(masterdir,"/bin/helper.R",sep=""))
source(paste(masterdir,"/bin/AuPairWise.R",sep=""))
load(paste(masterdir,"/data/pairs.Rdata",sep=""))


# Load your RNA-seq experiment here
# Requires an expression table as a matrix, column names with sample IDs, and row names as gene entrez IDs.
# If your data is in expression format, please convert to a matrix.
# load(file)

name = "my_experiment"
out = paste(masterdir,name,sep="/")

summary = run_APW(exprs, out, stoich.pairs )

### Output summary:
#### The summary variable is a list that contains 6 elements
#### The first is a table with the average AUROCs for the noise factors that the model ran.
#### The first row has the stochiometric pair results, and the second row has the random pairs run.

# summary$aurocs

#### The next two lists are the statistics for the above AUROCs (standard deviation and standard error):

# summary$aurocs.sd
# summary$aurocs.se

#### The noise factors are listed in:

# summary$n.factors

#### The p-valaue of the Wilcoxon rank sum test between the stoichiometric pairs and the random pairs for each noise factor are in:

# summary$pvals

#### And finally, the estimated noise factors for any AUROC is given in:

# summary$stats


#### You can view the summary results with:
# plot_summary_results(summary)

####  Or output the plot to a .png file:
# plot_summary_results(summary, out)


#### To run on a different list of gene pairs, modify the stoich.pairs variable.
#### This a 2D matrix, the first column has gene A and second column gene B.
#### These need to be entrez gene IDs, or must match the gene labels of your expresssion dataset.
# summary = run_APW(exprs, out, pairs )



