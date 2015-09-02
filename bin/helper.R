############################################
#	Helper functions for AuPairWise    #
############################################


##  Written: Sara Ballouz
##  Date: September 29th, 2014

##  Updates:
##  Date: September 2nd, 2015


## Necessary libraries
library(MASS)
library(Matrix)
library(scales)
library(gplots)
library(zoo)



lm.studentized  <- function(x,y)
{
	z = lm(y ~ x )
        z = rstudent(z)
	return( rank(abs(z)) )
}

lm.function  <- function(x,y)
{
	z = lm(y ~ x )
	return( rank(abs(z$residuals)) )
}

get_value <- function( x1, x2, y1,y2, x) 
{
	m = (y2 - y1) / (x2 - x1 )
	y = y1 + m *( x - x1)
	return(y)
}

model.fx <- function(x, fx){
	x = fx(x)
	replace = ( x == -Inf)
	x[replace] = 0
	return (x)
}


get_approx_expression <- function( temp.x, temp.x.r, temp.x.s, max.r, min.r)
{
	if( temp.x.s < min.r ){
		x0 = 1
        	x1 = 2
        	xn = min.r
       	} else if( temp.x.s > max.r ){
               	x0 = max.r-1
        	x1 = max.r
        	xn = max.r

	}  else {
		x0 = floor(temp.x.s)
        	x1 = ceiling(temp.x.s)
        	xn = temp.x.s
	}
        new.x = get_value(  x0, x1, temp.x[which( x0 ==  temp.x.r )] ,temp.x[which( x1 ==  temp.x.r )], xn)
	return (new.x)
}


shuffle <- function( X, s, n.factor, dist="other", mode="post"){

        nX1 = dim(X)[1]
	nS1 = dim(X)[2]
	X.r = t(apply ( X,1,rank, ties.method="random", na.last="keep") )   # method breaks down with "average" tied ranks

	if( n.factor != 0){
		if( dist == "gauss" ) {
                	noise = rnorm(nX1)
                	X.s = (X.r[,s]*((100-n.factor)/100) + nS1*noise*(n.factor/100))
                	# X.s = (X.r[,s] + nS1*noise*(n.factor/100))
                	# X.s = X.r[,s] + nS1*noise*(n.factor/100)


		} else {
			noise = runif(nX1, min=-n.factor/100, max=n.factor/100)
			X.s = X.r[,s] + nS1*noise
		}


        	X.e = (sapply( 1:nX1, function(i) get_approx_expression(X[i,], X.r[i,], X.s[i], nS1, 1), simplify=T))  # this does not work with tied ranks
                if( mode == "post") {
                        o1 = order(X.e)
                        o2 = order(X[,s])
                        X[o2,s] = X[o1,s]
		} else {
                        X[,s] = X.e
                }
 	 }
	 return(X)

}


predict_sample <- function( X, s,n.factor, k, nS, nK, filter, dist="other", mode="post" ){

	nX = dim(X)[1]
	X.new  = shuffle(X,s, n.factor, dist,mode)
	X.r    = t(apply (X.new,1,rank, ties.method="average", na.last="keep") )
	Z.ranks = t(sapply( 1:nK, function(i) lm.studentized( X.r[k[i,1],],X.r[k[i,2],]), simplify=T ))

	k2 = t(apply(all_pairs ( cbind(sample(nX,nK*2), sample(nX,nK*2))), 1, as.numeric))
     	Z.ranks.r = t(sapply( 1:nK, function(i) lm.studentized( X.r[k2[i,1],],X.r[k2[i,2],]), simplify=T ))

	ranks = matrix(unlist(Z.ranks) , ncol=nS, nrow=nK, byrow=F)
	all = ranks == nS

	scores    = rank(colSums(all), ties.method="random")
	labels    = scores*0
	labels[s] = 1

        ranks.r = matrix(unlist(Z.ranks.r) , ncol=nS, nrow=nK, byrow=F)
	all.r = ranks.r == nS

	scores.r    = rank(colSums(all.r), ties.method="random")

        scores = c(scores,scores.r)

	if( !missing(filter)){
		scores.filters =  list()

		for (i in 1:(length(filter))){
                	temp = rank(colSums(all[filter[[i]],]))
                        temp.r = rank(colSums(all.r[filter[[i]],]))

                        scores.filters = c(scores.filters, temp, temp.r)
                }
                scores = unlist(scores.filters)
	}

	return( c(labels,scores) )

}

predict_sample_with_replicates <- function( X,Y, s, n.factor, nX, nS, dist="other", mode="post" ){

	X.new  = shuffle(X,s, n.factor, dist,mode)
	X.r    = t(apply (X.new,1,rank, ties.method="average", na.last="keep") )
        Y.r    = t(apply (Y,1,rank, ties.method="average", na.last="keep") )

	Z.ranks = t(sapply( 1:nX, function(i) lm.studentized( X.r[i,],Y.r[i,]), simplify=T ))

	ranks = matrix(unlist(Z.ranks) , ncol=nS, nrow=nX, byrow=F)
	all = ranks == nS

	scores    = rank(colSums(all), ties.method="random")
	labels    = scores*0
	labels[s] = 1

	return( c(labels,scores) )

}

plot_samples_noise_aurocs <- function( data, S )
{
	range = data$n.factors[-1]     # remove 0 factor
	cols = colorpanel(nS, "red", "purple", "blue")

	plot(log10(range), data$aurocs[-1] ,ylim=c(0.4,1), type="l", lwd=3, col=0, xlab="Noise factor", ylab="AUROC", axes=F)
        axis(2)
        axis(1, lab=range, at=log10(range) )

        for (i.s in 1:length(S)) {
        	fit = glm( data$aurocs.s.comb[-1,i.s] ~ range)
                points( log10(range), data$aurocs.s.comb[-1,i.s], pch=19, col=cols[i.s])
                segments( log10(range), data$aurocs.s.comb[-1,i.s]-data$aurocs.s.comb.se[-1,i.s], log10(range), data$aurocs.s.comb[-1,i.s]+data$aurocs.s.comb.se[-1,i.s], col=makeTransparent(cols[i.s]))
		lines( log10(range), fit$fitted, lwd=2, col=makeTransparent(cols[i.s]))
	}
}





get_replicability <- function( data, AUROC ){
	range = data$n.factors[-1]   # remove 0 factor
	fit = glm( data$aurocs[-1] ~ range)
        predictions = predict(fit, data.frame(range=1:100),type = "response")
        max.pred = max(which(predictions <= AUROC))
        min.pred = min(which(predictions >= AUROC))
        n.pred = get_value_x( max.pred, min.pred, predictions[max.pred], predictions[min.pred], AUROC)

        plot(log10(1:100), predictions, ylim=c(0.4,1), type="l", lwd=3, col="lightgrey", xlab="Noise factor", ylab="AUROC", axes=F)
        axis(2)
        axis(1, lab=data$n.factors, at=log10(data$n.factors) )
        abline(h= data$aurocs[1], lty=2, col="lightgrey", lwd=2)
        points( log10(data$n.factors), data$aurocs, pch=19)
        abline( h = AUROC, lwd=3, lty=3, col=2)
        abline( v = log10(n.pred), lwd=3, lty=3, col=2)
	segments( log10(data$n.factors), data$aurocs-data$aurocs.se, log10(data$n.factors), data$aurocs+data$aurocs.se)

        return(n.pred)

}



get_indices_stoich_pairs <- function( stoich.pairs, genes.list )
{
	if( missing(stoich.pairs) ){ }
	if( missing(genes.list)  ) { }

        m = match( stoich.pairs[,1], genes.list )
	p1 = !is.na(m)
	x1 = m[p1]

	m = match( stoich.pairs[p1,2], genes.list )
	p2 = !is.na(m)
	x2 = m[p2]

	m = match( stoich.pairs[p1,1][p2], genes.list )
	p3 = !is.na(m)
	x1 = m[p3]

	indices = list()
	indices$p1 = p1
	indices$p2 = p2
	indices$x1 = x1
	indices$x2 = x2

	return (indices)
}



plot_cummulative_counts <- function(out, X)
{

	xN    = dim(X)[2]
	Xtemp = apply(X, 2, sort, decreasing=T )
	XSum  = apply(Xtemp, 2, cumsum)
	col=colorpanel(xN, "red","purple","blue")
	maxx = max(XSum)

        png( paste( out, ".all.csum", ".png", sep=""))
	plot(XSum[,1], col=0, ylim =c(0, maxx), xlab="Cummulative gene count", ylab="Expression levels covered", bty="n")
	for (i in 1:xN){
		lines(XSum[,i], col=col[i])
	}
	dev.off()

        png( paste( out, ".csum", ".png", sep=""))
	plot(XSum[,1]/max(XSum[,1]), col=0, ylim =c(0, 1), xlab="Cummulative gene count", ylab="Proportion of expression levels covered", bty="n")
	for (i in 1:xN){
		lines(XSum[,i]/max(XSum[,i]), col=col[i])
	}
        abline(h=0.5,lty=2, col="lightgrey")
	abline(h=0.75,lty=3, col="lightgrey")
	abline(h=0.95,lty=2, col="lightgrey")
	dev.off()


}

plot_expression_props <- function( out, m.X, sd.X, genes.stoich)
{
	if( missing(genes.stoich)){
		png( paste( out, ".exprs.png", sep=""))
	} else {
        	png( paste( out, ".exprs.pairs.png", sep=""))
	}

	plot( (m.X), (sd.X), pch=19, cex=0.1, col=makeTransparent(1), xlab="Average expression", ylab="SD", bty="n")

	if( !missing(genes.stoich)){
		points( m.X[genes.stoich], sd.X[genes.stoich], pch=19, col=makeTransparent("red"))
	}
	dev.off()

}

plot_expression_density <- function( out, X, xlab)
{

	png( paste( out, ".density.png", sep="") )
	plot_density("", xlab, 0.1, X)
	dev.off()
}


plot_density <- function(labels, xlab, b, a)
{
        n = dim(a)[2]

	l = density( (unlist(a)), bw=b, na.rm=T)
        all = cbind(l$x, l$y)

        for (i in 1:n){
                l = density((a[,i]), bw=b)
                all = cbind(all, l$x, l$y)
        }

        # Setup plot variables
        odd = (1:(dim(all)[2]/2))*2 - 1
        even = (1:(dim(all)[2]/2))*2
        ymax=max( all[,even])
        ymin=min( all[,even])
        xmin=min( all[,odd])
        xmax=max( all[,odd])

    	cols = colorpanel(n, "red", "blue", "green")

        legcols = list()

        plot(-xmax,-ymax, type="l", lwd=2, col="grey", ylim=c(ymin,ymax),xlim=c(xmin,xmax), ylab="Density", xlab=xlab)
        j = 1
	for (i in odd){
                lines( all[,c(i,i+1)], lwd=0.1, col=cols[j])
		legcols = append(legcols, cols[j])
                j = j + 1
        }
}

plot_stoich_cors <- function(out, length, filter, pairs, X){
	cols = colorpanel(length, "red", "purple", "blue")
	stats = list()
	h.all = list()
	m = "spearman"
        cors = diag( cor( t(X[indices$x1,]), t(X[indices$x2,]), method=m))
        r1 = sample(indices$x1)
        r2 = sample(indices$x2)

	cors.r = diag(cor( t(X[r1,]), t(X[r2,]), method=m))

	# Total
       	h.total = get_density(hist(cors.r, plot=F))
       	m.total = mean(cors.r, na.rm=T)
        y.max = max(h.total[,2])

	for( i in 1:length ){
		h = get_density( hist(cors[filter[[i]]], breaks=50, plot=F ))
		m = mean(cors[filter[[i]]], na.rm=T)
		stats[[i]] = m
		h.all[[i]] = h
		if( y.max < max(h[,2]) ){ y.max=max(h[,2]) }
	}

	png( paste( out, ".hist.cors.pairs.png", sep=""))
	plot(h.total, type="l", xlim=c(-1,1), ylim=c(0,y.max), col=makeTransparent(1,5), lwd=4, xlab="Correlation coefficients", ylab="Density" )
	abline(v = m.total, lty=2, col=makeTransparent(1,5))
	for( i in 1:length ){
		lines( h.all[[i]], col=cols[i], lwd=2)
		abline(v = stats[[i]], lty=2, col=makeTransparent(cols[i]) )
	}
	legend( "topleft", legend=pairs$labels, col=cols, lwd=2, bty="n")

	dev.off()
}


write_out_summary <- function(out, results.all, length, pairs, n.factors, n.repeats, AUROC.default = 0.8, nn=100){

	# Summary results
	aurocs = matrix(0, nrow=length*2, ncol=length(n.factors), dimnames=list( array(rbind( pairs$labels, "Random") ) ,n.factors))
	aurocs.sd = aurocs
	aurocs.se = aurocs
        pvals = matrix( 0, nrow=length, ncol=length(n.factors), dimnames=list( pairs$labels,n.factors) )

	i = 1
	for( n.factor in n.factors){
        	data = matrix(unlist(results.all[[i]]$aurocs), ncol=n.repeats , byrow=T )
		aurocs[,i] = rowMeans( data )
                aurocs.sd[,i] = apply(data, 1,sd)
                aurocs.se[,i] = aurocs.sd[,i]/sqrt(dim(data)[2])
                pvals[,i] = sapply( (1:length)*2 -1, function(j) wilcox.test(data[j,], data[j+1,])$p.val )
                i = i + 1

        }

        data = list()
        data$aurocs = aurocs
        data$aurocs.se = aurocs.se
        data$aurocs.sd = aurocs.sd
        data$n.factors = n.factors
        data$pvals = pvals

	write.table (data$aurocs, file = paste(out,".avg.aurocs.summary", sep=""))
	write.table (data$aurocs.sd, file = paste(out,".avg.aurocs.summary", sep=""), append=T, col.names=F, row.names=paste(array(rbind( pairs$labels, "Random") ), "- SD"))
	write.table (data$aurocs.se, file = paste(out,".avg.aurocs.summary", sep=""), append=T, col.names=F, row.names=paste(array(rbind( pairs$labels, "Random") ), "- SE"))


	# Predictions
	stats = matrix(NA, ncol=nn, nrow=length*2, dimnames=list( array(rbind( pairs$labels, "Random"), 1:nn /nn ) )
	range = n.factors[-1]


	for ( i in 1:(length*2) ){
		temp = data$aurocs[i,][-1]
                temp.se = data$aurocs.se[i,][-1]
		temp.sd = data$aurocs.sd[i,][-1]

		predictions = approx( n.factors, data$aurocs[i,], n=nn)

	        for (j in (nn/2):nn){
	        	AUROC = j/nn
			max.pred = max(which(predictions$y <= AUROC))
	        	min.pred = min(which(predictions$y > AUROC))
	        	n.pred = get_value_x( predictions$x[max.pred], predictions$x[min.pred], predictions$y[max.pred], predictions$y[min.pred], AUROC)
   	        	stats[i,j] = n.pred
	        }

        }

	data$stats = stats
	# write.table (data$stats, file = paste(out,".avg.aurocs.predictions", sep="")  )

        save(data, results.all, file=paste(out,".avg.aurocs.Rdata", sep="") )
        return(data)
}


plot_summary_results <- function(data, out, AUROC.default=0.8){
	if( !missing(out) ){
	        png( paste(out, ".predictions.png", sep="")  )
	}

        col.def = which(colnames(data$stats)==AUROC.default)

	n = dim(data$aurocs)[1]

	plot( log10(data$n.factors), data$aurocs[1,], ylim=c(0.4,1), type="l", lwd=3, col=0, xlab="Noise factor", ylab="AUROC", axes=F)
        axis(2)
        axis(1, lab=data$n.factors, at=log10(data$n.factors) )
        cols = makeTransparent(colorpanel(n, "black", "lightgrey"),150)

	for( i in 1:n){
        	lines( log10(data$n.factors), data$aurocs[i,], lwd=3, col=cols[i])
		segments(log10(data$n.factors), data$aurocs[i,]-data$aurocs.se[i,], log10(data$n.factors),data$aurocs[i,]+data$aurocs.se[i,],col=cols[i])
        	points( log10(data$n.factors), data$aurocs[i,], pch=19,col=cols[i])
        }

	abline( h = AUROC.default, lwd=3, lty=3, col=2)
        abline( v = log10(data$stats[1,col.def]), lwd=3, lty=3, col=2)
	text(log10(max(data$n.factors)/2), 0.5, paste("For an AUROC of:", AUROC.default, "\nthe estimated noise factor is:\n", round(data$stats[1,col.def],1), "%") )
        legend("topleft", legend=rownames(data$aurocs), col=cols, lwd=3, pch=19)

	if( !missing(out) ){
	       	dev.off()
	}
}


get_gene_corrs <- function(X, m )
{
	gene.corr = cor(t(X), method=m)
	#na = is.na(gene.corr)
	#gene.corr[na] = 0
	#diag(gene.corr) = 1
	return(gene.corr)
}


# Transparent colors

makeTransparent <- function(someColor, alpha=100)
{
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
	blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Basic linear algebra

# Given x and two points
get_value <- function( x1, x2, y1,y2, x) {
	m = (y2 - y1) / (x2 - x1 )
	y = y1 + m *( x - x1)
	return(y)
}

# Given y and two points
get_value_x <- function( x1, x2, y1,y2, y) {
	m = (y2 - y1) / (x2 - x1 )
	x = x1 + (y - y1)/m 
        return(x)
}

## Formats the density distribution from the histogram function
get_density <- function(hist)
{
        x = sort(rep(hist$breaks,2))
        y = matrix(rbind(hist$density, hist$density))
        y = c(0,y,0)

        return(cbind(x,y))
}


## Formats the counts distribution from the histogram function
get_counts <- function(hist)
{
        x = sort(rep(hist$breaks,2))
        y = matrix(rbind(hist$counts, hist$counts))
        y = c(0,y,0)

        return(cbind(x,y))
}

# List of pairs to a matrix
make_network_from_data<- function(pairs, genes){
	n = length(genes)
	net = matrix(0, ncol=n, nrow=n)

	m = match( as.numeric(pairs[,1]), genes  )
	p1 = !is.na(m)
	m1 = m[p1]

	m = match( as.numeric(pairs[p1,2]), genes )
	p2 = !is.na(m)
	m2 = m[p2]

	net[cbind(m1,m2)] = 1
	net[cbind(m2,m1)] = 1
	diag(net) = 0

	colnames(net) = genes
        rownames(net) = genes
	return(net)
}

# List of pairs to a matrix
make_network_from_data<- function(data, listA, listB){
	nr = length(listA)
	nc = length(listB)

	net = matrix(0, ncol=nc, nrow=nr)

	m = match( (data[,1]), listA  )
	p1 = !is.na(m)
	m1 = m[p1]

	m = match( (data[p1,2]), listB )
	p2 = !is.na(m)
	m2 = m[p2]

	net[cbind(m1,m2)] = 1

	colnames(net) = listB
        rownames(net) = listA
	return(net)
}

# Tic toc functions
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
	type <- match.arg(type)
	assign(".type", type, envir=baseenv())
	if(gcFirst) gc(FALSE)
	tic <- proc.time()[type]
	assign(".tic", tic, envir=baseenv())
	invisible(tic)
}

toc <- function()
{
	type <- get(".type", envir=baseenv())
	toc <- proc.time()[type]
	tic <- get(".tic", envir=baseenv())
	print(toc - tic)
	invisible(toc)
}


## Calculates the ROC for a given list of scores(or ranks) and the known labels
roc_score <- function(scores,labels)
{
	negatives = which(labels == 0, arr.ind=T)
	scores[negatives] <- 0

	p  = sum(scores)
	nL = length(labels)
	np = sum(labels)
        nn = nL - np

	roc  <-  (p/np - (np+1)/2)/nn

	return(roc)
}

# Calculates and plots the ROC for a given list of scores(or ranks) and the known labels
plot_roc2 <- function(scores,labels) {
	o = order(scores, decreasing=T)

	h1 = labels[o]
        h2 = !labels[o]

	auroc  <-  roc_score(scores, labels)
	tpr = c(0,cumsum(h1)/sum(h1))
	fpr = c(0,cumsum(h2)/sum(h2))
	plot(fpr,tpr, xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), type="l" )
	lines( c(0,1), c(0,1), col="grey")
	text(labels=auroc, 0.5, 1 )
        pval = wilcox.test(scores[labels], scores[!labels])
	return(cbind(fpr,tpr))
}

# Calculates and plots the ROC for a given list of scores(or ranks) and the known labels
plot_roc <- function(scores,labels, file)
{
        rocs = get_roc_curve(scores, labels)
	auroc  <-  roc_score(scores, labels)
	tpr = rocs$tpr
        fpr = rocs$fpr
	png(file)
	plot(fpr,tpr,  type="l", lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="FPR", ylab="TPR")
	lines( c(0,1), c(0,1), col="grey")
	dev.off()
        pval = wilcox.test(scores[labels], scores[!labels])
	return(cbind(fpr,tpr))
}


get_roc_curve <- function(scores,labels){
        #dups = rev(duplicated(rev(scores)))
        #cutoffs = c(Inf, scores[!dups])
	o = order(scores, decreasing=T)
        cutoffs = c(Inf, length(scores):1)
        roc = sapply( (1:length(cutoffs)), function(i) calc_rates(scores[o], labels[o],cutoffs[i]) , simplify=F)
        rocs = matrix(unlist(roc), byrow=T, ncol=2)
        colnames(rocs) = c("fpr", "tpr")
        return(rocs)
}

get_avgroc_curve <- function( rocs,n.repeats, n){

	sum = matrix(0, ncol=2, nrow = n)
	colnames(sum) = c("tpr", "fpr")
	for( i in 1:n.repeats) {
		sum = rocs[[i]] + sum
	}
	sum = sum/n.repeats
        return(sum)
}



get_auc <- function(x,y){
        o = order(x)
        auc <- sum(diff(x[o])*rollmean(y[o],2))
        return (auc)
}

make_hist_matrix <- function(x1,y1){
	x = sort(rep(x1,2))
        y = matrix(rbind(y1, y1))
        y = c(0,y,0)

        return(cbind(x,y))
}

calc_rates <- function(scores, labels, thresh){
        tp = sum((scores >= thresh) * labels)
        fp = sum ((scores >= thresh) * !labels)
        p  = sum(labels)
        n = sum(!labels)
        #fn = p - tp
        #tn = n - fp
        fpr = (fp/n)
        tpr = (tp/p)

        #return(cbind(tp,fp, tn, fn, p, n, fpr,tpr))
        return(cbind(fpr,tpr))
}


all_pairs <- function(list){
	list = list[(list[,1] != list[,2]),]
	list = sort(c( paste( list[,1],list[,2]), paste( list[,2], list[,1])))
        list = matrix ( unlist(strsplit(as.character(list), " ") ), ncol=2, nrow= length(list), byrow=T )
	list = unique( t(apply( list, 1, sort) ))
	return(list)
}


unique_all_pairs <- function(pairs){
	list = pairs[[1]]
	if( length(pairs) < 2 ) { return(list) }

	for( i in 2:length(pairs)){
        	list = rbind(list,pairs[[i]])
	}
	list = unique(list)
	return(list)
}

filter_pairs <- function(pairs, indices, length ){

	temp = pairs$all[indices$p1,][indices$p2,]
	indexed = paste(temp[,1], temp[,2])

	filter = list()
	for( i in 1:length){
        	m = match(indexed, paste(pairs[[i]][,1], pairs[[i]][,2]))
		filt = !is.na(m)
		filter[[i]] = filt
	}
	return(filter)
}



get_expression_levels <- function( indices, X, fx){
	
	if( missing(fx) ){
        	xp1 = X[indices$x1,]
        	xp2 = X[indices$x2,]
	} else {
		xp1 = model.fx(X[indices$x1,], fx)     # protein 1, in rep X
		xp2 = model.fx(X[indices$x2,], fx)     # protein 2, in rep X
	}
	results = list()
	results$xp1 = xp1
	results$xp2 = xp2
	return(results)
}

