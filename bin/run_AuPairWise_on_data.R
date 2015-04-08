masterdir = ""

source(paste(masterdir,"/bin/helper.R",sep=""))
source(paste(masterdir,"/bin/AuPairWise.R",sep=""))
load(paste(masterdir,"/data/pairs.Rdata",sep=""))


out = paste(masterdir,"",sep="")

# Load your RNA-seq experiment here 
# Requires an expression table as a matrix, column names with sample IDs, and row names as gene entrez IDs
# If your data is in expression format, please convert to a matrix.

run_APW(exprs, out, stoich.pairs )
