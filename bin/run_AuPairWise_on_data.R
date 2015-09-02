#################################################################################################
#	AuPairWise: A tool to assess biological replicability without using replicates		#
#################################################################################################

##  Written: Sara Ballouz
##  Date: September 29th, 2014

# Setup environment variables
if( missing(masterdir) ) {
   masterdir = "H:/AuPairWise/"
}
if( missing(name) ) {
   name = "my_experiment"
}
if( missing(out) ) {
   out = paste(masterdir,name,sep="/")
}


# Load helper files and stoichiometric pairs
source(paste(masterdir,"/bin/helper.R",sep=""))
source(paste(masterdir,"/bin/AuPairWise.R",sep=""))
load(paste(masterdir,"/data/pairs.Rdata",sep=""))

# Run AuPairWise
if( missing(exprs) ) {
   print ("Please load experiment!")
   exit;

} else {
	summary = run_APW(exprs, out, stoich.pairs )
}

