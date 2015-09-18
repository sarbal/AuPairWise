args <- commandArgs(trailingOnly = TRUE)
        dir  = args[1]
        indir   = args[2]
        bin  = args[3]
	file = args[4]
        j = as.numeric(args[5])
	flag    = args[6]
	ranked  = args[7]
	random  = args[8]
	n_samples = args[9]

        source(paste(bin,"helper_aupair_develp.r",sep=""))
	load(paste(dir,"pairs.Rdata", sep=""))

# Change noise factors here
	n.factors = c(0,1,2,5,7,10,15,20,25,50,100)
	n.repeats = 100

# Adjust list of pairs
	pairs = list()
	pairs$stoich = all_pairs(stoich.pairs)
	pairs$all = unique_all_pairs( pairs )
        pairs$labels  = c("Stoichiometric pairs")

	length = length(pairs$labels)

# Set flags and conditions
# ranked = FALSE
# random = FALSE
# n_samples = "multi"


# Select data
mat = exprs.rseq
A = (1:21)* 2 -1
B = 1:21 * 2


# Get data
if( flag == "micr"){
        load(paste(indir,"sample_brainspan.micr.Rdata", sep=""))

        dir =paste(indir, "/brainspan/", sep="")
	label="brainspan.micr"

	X = exprs
        filterout = colSums(X) != 0 

	X = X[,filterout]

        samples.list = colnames(X)
        genes.list = rownames(X)

      	out = paste(dir,label,j, sep="")


} else if( flag == "rna"){
        load(paste(indir,"sample_brainspan.Rdata", sep=""))

        dir =paste(indir, "/brainspan/", sep="")
	label="brainspan.rna"

	X = exprs
        filterout = colSums(X) != 0

	X = X[,filterout]

        samples.list = colnames(X)
        genes.list = rownames(X)

      	out = paste(dir,label,j, sep="")

}

# Store size of experiment (samples, genes, etc)

 	NN = dim(X)[2]    		# Number of samples
	N = dim(X)[1]     		# Number of genes/transcripts
	S = 1:NN          		# Indices for samples
	nS = NN                         # If subsampling, currently not implemented



# Adjust data (ranked, random or none)
if( ranked == T){
	X = apply(X, 2 ,rank, ties.method="average", na.last="keep")
	colnames(X) = samples.list
        rownames(X) = genes.list
        out=paste(out, "ranked", sep=".")

} else if( random == T){
	X.random = sapply( 1:N, function(i) X[i,sample(NN)], simplify=T)
        X.random = matrix(unlist(X.random) , ncol=NN, nrow=N, byrow=F)

	X = apply(X.random, 2 ,rank, ties.method="average", na.last="keep")
	colnames(X) = samples.list
        rownames(X) = genes.list
        out=paste(out, "randomized", sep=".")
}

# Get indices of pairs
        indices =  get_indices_stoich_pairs(pairs$all, genes.list)
	indices.stoich = get_indices_stoich_pairs(pairs$stoich, genes.list)
        nX = length(indices$x1)
        genes.stoich = sort(unique(c(indices.stoich$x1, indices.stoich$x2)))
        x = cbind( indices$x1, indices$x2)

	filter = filter_pairs(pairs, indices,length)


# 2fold validation, not implemented for replicate analysis
if( n_samples == "multi") {
      	out = paste(out, n_samples, sep=".")
        sS = sapply(1:n.repeats, function(i) sort(sample( NN, NN/2 )) )

	results.all = list()

	k = 1
	for (n.factor in n.factors) {
		tic()
		repeats = list()

		print(paste("Noise factor: ", n.factor))
                print("Samples: ")
		print(sS)

		repeats$noise = sapply((1:n.repeats), function(i) predict_sample_multi(X, sS[,i], n.factor, x , nS, nX, filter) , simplify=F)
		for ( j in 1:(length*2)){
			repeats$rocs[[j]]   = sapply( (1:n.repeats), function(i) get_roc_curve(repeats$noise[[i]][(1:nS)+(nS*j)] ,repeats$noise[[i]][(1:nS)+(nS*0)]), simplify=F)
	                repeats$aurocs[[j]] = sapply( 1:n.repeats, function(i) get_auc(repeats$rocs[[j]][[i]][,1], repeats$rocs[[j]][[i]][,2]))
	                repeats$avgroc[[j]] = get_avgroc_curve( repeats$rocs[[j]], n.repeats, nS+1)
	        }


	        temp = matrix( unlist(repeats$aurocs), nrow=(length*2), ncol=n.repeats, byrow=T)
	        rownames(temp) = array(rbind( pairs$labels, "Random") )

		repeats$stats = sapply( ((1:length)*2)-1, function(i) wilcox.test( temp[i,],temp[i+1,])$p.val )

		# Write out results
		#write.table( temp, file=paste(out, ".sample.", nS, ".noise.", n.factor,".avg.aurocs", sep=""), col.names=F)
	        #write.table( matrix( unlist(repeats$avgroc), nrow=nS+1, ncol=length*2, byrow=F), file=paste(out, ".sample.", nS, ".noise.", n.factor,".avg.aurocs.fpr.tpr", sep=""), col.names=F)
	        #write.table( matrix(unlist(repeats$noise), nrow=n.repeats, byrow=T), file=paste(out,  ".sample.", nS, ".noise.", n.factor,".labels.scores", sep=""), col.names=F)
                save(n.factors, n.repeats, repeats, file=paste(out,  ".sample.", nS, ".noise.", n.factor,".Rdata", sep=""))

		# Store results
		results.all[[k]]= repeats
		k = k + 1
		toc()
	}

	# Summary results
        write_out_summary(out, results.all, length, pairs, n.factors, n.repeats)


# n-fold validation
} else {
	results.all = list()

	k = 1
	for (n.factor in n.factors) {
		tic()
		repeats = list()

		print(paste("Noise factor: ", n.factor))
		shuff = sample(nS, n.repeats, replace=T)
	        subS =  t(sapply((1:n.repeats), function(i) sort(sample(NN, nS))))
	        
		repeats$noise = sapply((1:n.repeats), function(i) predict_sample(X[,subS[i,]], shuff[i], n.factor, x , nS, nX, filter) , simplify=F)
		for ( j in 1:(length*2)){
			repeats$rocs[[j]]   = sapply( (1:n.repeats), function(i) get_roc_curve(repeats$noise[[i]][(1:nS)+(nS*j)] ,repeats$noise[[i]][(1:nS)+(nS*0)]), simplify=F)
		        repeats$aurocs[[j]] = sapply( 1:n.repeats, function(i) get_auc(repeats$rocs[[j]][[i]][,1], repeats$rocs[[j]][[i]][,2]))
		        repeats$avgroc[[j]] = get_avgroc_curve( repeats$rocs[[j]], n.repeats, nS+1)
		}
		temp = matrix( unlist(repeats$aurocs), nrow=(length*2), ncol=n.repeats, byrow=T)
		rownames(temp) = array(rbind( pairs$labels, "Random") )

		repeats$stats = sapply( ((1:length)*2)-1, function(i) wilcox.test( temp[i,],temp[i+1,])$p.val )

		# Write out results
		#write.table( temp, file=paste(out, ".sample.", nS, ".noise.", n.factor,".avg.aurocs", sep=""), col.names=F)
	        #write.table( matrix( unlist(repeats$avgroc), nrow=nS+1, ncol=length*2, byrow=F), file=paste(out, ".sample.", nS, ".noise.", n.factor,".avg.aurocs.fpr.tpr", sep=""), col.names=F)
	        #write.table( matrix(unlist(repeats$noise), nrow=n.repeats, byrow=T), file=paste(out,  ".sample.", nS, ".noise.", n.factor,".labels.scores", sep=""), col.names=F)
                save(n.factors, n.repeats, repeats, file=paste(out,  ".sample.", nS, ".noise.", n.factor,".Rdata", sep=""))

		# Store results
		results.all[[k]]= repeats
		k = k + 1
		toc()
	}

	# Summary results
        write_out_summary(out, results.all, length, pairs, n.factors, n.repeats)

} 
