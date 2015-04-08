masterdir = ""
source(paste(masterdir,"/bin/helper.R",sep=""))
source(paste(masterdir,"/bin/AuPairWise.R",sep=""))

load(paste(masterdir,"/sample/sample_brainspan.Rdata",sep=""))
load(paste(masterdir,"/data/pairs.Rdata",sep=""))

out = paste(masterdir,"/sample/sample_brainspan",sep="")

run_APW(exprs, out, stoich.pairs )
