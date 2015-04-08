#########################################
#	AuPairWise for expression data  #
#########################################

##  Written: Sara Ballouz
##  Date: September 29th, 2014
##  Modified: April 8th 2015
## 	- Renamed

run_APW <- function(exprs, out, stoich.pairs,  n.factors=c(0,1,2,5,10,15,25,50,100), n.repeats=10, dist ="other", mode ="post" , ranked=FALSE ){

	# Filter data
	X = exprs

	# Remove samples with no expression data
        filterout = colSums(X) != 0
	X = X[,filterout]

	# Remove genes with too few counts
	X[which( log10(X) < -5 )] = 0

	# Set up variables
        NN = dim(X)[2]    		# Number of samples
	N = dim(X)[1]     		# Number of genes/transcripts
	S = 1:NN          		# Indices for samples
	nS = NN                         # If subsampling, currently not implemented

	# Visualize data so far
        plot_cummulative_counts(out, X)

        # Transform to log2
        # Med <- median(X, na.rm = T)
	# if (Med > 16) 
	X <- model.fx(X, log2)

	# Transform data
	if( ranked == T){
		X = apply(X, 2 ,rank, ties.method="average", na.last="keep")
		colnames(X) = samples.list
        	rownames(X) = genes.list
        	out=paste(out, "ranked", sep=".")
		plot_cummulative_counts(out, X)
	}

	# Properties of expression dataset
	m.X = rowMeans(X, na.rm=T) 	# Mean expression of genes/transcripts across samples
	sd.X = apply(X,1,sd, na.rm=T)	# SD of genes/transcripts expression across samples
	plot_expression_props(out, m.X, sd.X)

	# Update data
        genes.list = rownames(X)
	samples.list = colnames(X)

	# Adjust list of pairs
	pairs = list()
	pairs$stoich         = all_pairs(stoich.pairs)
	# pairs$ppin           = all_pairs(ppin.pairs)
	# pairs$yeastmouse.top = all_pairs(yeastmouse.top.pairs)
	# pairs$random 	     = all_pairs( cbind( sample( c(stoich.pairs[,1],ppin.pairs[,1])) , sample( c(stoich.pairs[,2],ppin.pairs[,2])))[1:dim(stoich.pairs)[1],] )
	pairs$random 	     = all_pairs( cbind( sample(stoich.pairs[,1]) , sample( stoich.pairs[,2]) ) )
	pairs$all            = unique_all_pairs( pairs )
	#pairs$labels         = c("Stoichiometric pairs", "Protein-protein interacting pairs, GO complexes","Yeast-mouse overlap, top", "Random pairs" )
        pairs$labels         = c("Stoichiometric pairs", "Random pairs" )

	length = length(pairs$labels)


	# Get indices of pairs
        indices =  get_indices_stoich_pairs(pairs$all, genes.list)
	indices.stoich = get_indices_stoich_pairs(pairs$stoich, genes.list)
        nX = length(indices$x1)
        genes.stoich = sort(unique(c(indices.stoich$x1, indices.stoich$x2)))
        x = cbind( indices$x1, indices$x2)

	filter = filter_pairs(pairs, indices,length)
	plot_expression_props(out, m.X, sd.X,genes.stoich)

        # Plot correlation distributions of pairs
	# plot_stoich_cors(out, length, filter, pairs, X)


	# Calculate AUROCs for each noise factor, using each pair set
	results.all = list()
	k = 1
	for (n.factor in n.factors) {
		tic()
		repeats = list()
		print("Noise factor: ", n.factor)
		shuff = sample(nS, n.repeats, replace=T)
	        subS =  t(sapply((1:n.repeats), function(i) sort(sample(NN, nS))))
	        repeats$noise = sapply((1:n.repeats), function(i) predict_sample(X[,subS[i,]], shuff[i], n.factor, x , nS, nX, filter) , simplify=F)
		for ( j in 1:length){
			repeats$rocs[[j]]   = sapply( (1:n.repeats), function(i) get_roc_curve(repeats$noise[[i]][(1:nS)+(nS*j)] ,repeats$noise[[i]][(1:nS)+(nS*0)]), simplify=F)
	                repeats$aurocs[[j]] = sapply( 1:n.repeats, function(i) get_auc(repeats$rocs[[j]][[i]][,1], repeats$rocs[[j]][[i]][,2]))
	                repeats$avgroc[[j]] = get_avgroc_curve( repeats$rocs[[j]], n.repeats, nS+1)
	        }

	        temp = matrix( unlist(repeats$aurocs), nrow=length, ncol=n.repeats, byrow=T)
	        rownames(temp) = pairs$labels

		# Write out results
		write.table( temp, file=paste(out, ".sample.", nS, ".noise.", n.factor,".avg.aurocs", sep=""), col.names=F)
	        write.table( matrix( unlist(repeats$avgroc), nrow=nS+1, ncol=length*2, byrow=F), file=paste(out, ".sample.", nS, ".noise.", n.factor,".avg.aurocs.fpr.tpr", sep=""), col.names=F)
	        write.table( matrix(unlist(repeats$noise), nrow=n.repeats, byrow=T), file=paste(out,  ".sample.", nS, ".noise.", n.factor,".labels.scores", sep=""), col.names=F)

		# Store results
		results.all[[k]]= repeats
		k = k + 1
		toc()
	}

	# Summary results
        write_out_summary(out, results.all, length, pairs, n.factors, n.repeats)

}


