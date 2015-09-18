#################################################################################################
#	AuPairWise: A tool to assess biological replicability without using replicates		#
#################################################################################################

# Code for figures and analysis

##############################################################################################
# Figure 1 : Tissue specific replication
##############################################################################################

load(file="H:/AuPairWise/suppl/corrs-examples.Rdata")
i = 1
j = 4

# To plot comparison between liver expression within same experiment
plot( log2(exprs[,i]), log2(exprs[,i+1]),cex=0.5, main=round(cor(exprs[,i], exprs[,i+1], use="pair", method="p"),3),sub=round(cor(exprs[,i], exprs[,i+1], use="pair", method="s"),3),
xlab=paste(samples[i], "Liver mRNA - log2 FPKM") , ylab=paste(samples[i+1], "Liver mRNA - log2 FPKM"), bty="n", pch=19, col=makeTransparent(1))


# To plot comparison between liver expression and kidney expression within same experiment
plot( log2(exprs[,i]), log2(exprs[,j]),cex=0.5, main=round(cor(exprs[,i], exprs[,j], use="pair", method="p"),3),sub=round(cor(exprs[,i], exprs[,j], use="pair", method="s"),3),
xlab=paste(samples[i], "Liver mRNA - log2 FPKM") , ylab=paste(samples[j], "Kidney mRNA - log2 FPKM"), bty="n", pch=19, col=makeTransparent(1))

# To plot comparison between liver expression and kidney expression of different experiments
plot( log2(exprs[f3,i]), log2(rsem[f4,i]), cex=0.5,main=round(cor(exprs[f3,i], rsem[f4,i], use="pair", method="p"),3),sub=round(cor(exprs[f3,i], rsem[f4,i], use="pair", method="s"),3),
xlab=paste(samples[i],"Liver mRNA - log2 FPKM"), ylab=paste(colnames(rsem)[i], "Liver - log2 RPKM"), bty="n", pch=19, col=makeTransparent(1))

# To plot comparison  between liver expression and liver expression from the RNA-seq Atlas
plot( log2(exprs[f1,i]), log2(atlas$liver[f2]), cex=0.5,main=round(cor(exprs[f1,i], atlas$liver[f2], use="pair", method="p"),3),sub=round(cor(exprs[f1,i], atlas$liver[f2], use="pair", method="s"),3) ,xlab=paste(samples[i], "Liver mRNA - log2 FPKM"), ylab="RNA-seq atlas, liver - log2 RPKM", bty="n", pch=19, col=makeTransparent(1))


##############################################################################################
# Loading ENCODE data for figures 2-5
##############################################################################################

load(file="H:/AuPairWise/sample/sample_ENCODE.Rdata")

# The dataset is split into a sample and it's techincal replicate.
# We select for them here. There is no reason for one replicate to belong to either group A or B, so they are interchangeable.
A = (1:21)* 2 -1
B = 1:21 * 2

# We select the data from the RSEQtools estimation, which has gene's RPKM values.
# exprs.cuff: cufflinks (FPKM)
# exprs.htseq: HTseq (raw counts)

X = exprs.rseq[,A]
Y = exprs.rseq[,B]

# Filter the samples with missing data and its replicate
filterout = colSums(X) != 0 & colSums(Y) != 0
X = X[,filterout]
Y = Y[,filterout]

# Properties of dataset
NN = dim(X)[2]    # number of samples
N = dim(X)[1]     # number of genes
S = 1:NN          # indices for samples

# Gene-level replicability calculation ie. correlation between X and Y
rhoXY.p   = sapply( 1:N, function(i) cor(X[i,], Y[i,], method="p"), simplify=T)
rhoXY.s = sapply( 1:N, function(i) cor(X[i,], Y[i,], method="s"), simplify=T)

##############################################################################################
# Figure 2 : Simpson's paradox
##############################################################################################



##############################################################################################
# Figure 3:  genomic coverage and sample size
##############################################################################################

# Pearson version
m = "p"
rhoXY = rhoXY.p

# Spearman version
m = "s"
rhoXY = rhoXY.s

# Sample-sample density
rho = cor(X,Y, method=m)
rho.h = get_density( hist(rho, plot=F) )

# Gene-gene density
rhoXY.h =  get_density( hist(rhoXY, breaks=50, plot=F) )

# Plot correlation densities of sample-sample and gene-gene correlations
xmin =  min(rho.h[,1], rhoXY.h[,1], na.rm=T)
xmax = max(rho.h[,1], rhoXY.h[,1], na.rm=T)
plot (rho.h, xlab="Correlations coefficients",lwd=2, xlim=c(xmin,xmax),type="l", ylab="Density of gene/transcript correlations", bty="n", main =round(mean(rhoXY, na.rm=T),3), sub = round(mean(rho, na.rm=T),3))
polygon(rho.h, col=makeTransparent( "purple", 50 ))
par(new=T)
plot(rhoXY.h,axes=F, lwd=2,type="l", xlim=c(xmin,xmax), xlab="", ylab="")
polygon(rhoXY.h, col=makeTransparent("blue", 50 ))
axis(4 )
mtext("Density of sample correlations", side=4, padj=3)
abline( v=mean(rhoXY, na.rm=T), col=makeTransparent("purple", 50 ))
abline( v=mean(rho, na.rm=T),col=makeTransparent("blue", 50 ))

# Calculate expected p-values/rhos for given N
MAX=1000
pvals_adj  = matrix(0, ncol=MAX, nrow=length(rhoXY))
test = rhoXY
for( i in 3:MAX){
        pvals = calc_pvals( rhoXY, i)
        pvals.adj = p.adjust(pvals, method ="holm")
        pvals_adj[,i] =  pvals.adj
}
        pvals_test = (pvals_adj <= 0.05)
        scores_test = colSums(pvals_test, na.rm=T)
        scores_test[1] = 0
        scores_test[2] = 0
        pvals.adj = pvals_adj[,NN]
        P = length(pvals.adj[!is.na(pvals.adj)])
        filt.p = log10((sort(pvals.adj[!is.na(pvals.adj)]) )) < log10(0.05)
        sum(filt.p)/P


# Plot fraction of replicated transcripts
plot (log10(1:MAX), scores_test/P, xlab="Number of samples required (log10)", main="", lwd=2, type="l", ylab="Fraction of transcripts significantly replicated", bty="n")
abline(v = log10(NN), lwd=3,col=colors[round(NN/3)], lty=3)
abline(h = sum(filt.p)/P, lwd=3,col=colors[round(NN/3)], lty=3)
#abline(h = 0.5, lwd=3,col="lightgrey")
#abline(h = 0.9, lwd=3, col="lightgrey")
#abline(v = log10(min(which(scores_test/P > 0.5))), lwd=3,col="lightgrey")
#abline(v = log10(min(which(scores_test/P > 0.9))), lwd=3,col="lightgrey")


load(file="H:/AuPairWise/suppl/GEO_series_2014.Rdata")

# Calculate the fraction of replicable genes per experiment
h = get_counts(hist(scores_test[all_exp_parsed[,2]]/P, plot=F, breaks=100 ))
h[,2] = h[,2]/length(all_exp_parsed[,2])

# Plot coverage
plot(h, main="", lwd=2, type="l", xlab="Fraction of genes significantly replicated", bty="n", ylab="Proportion of experiments" )
polygon(h, col=makeTransparent("blue"))



##############################################################################################
# Figure 4:  self versus coexp
##############################################################################################

# WARNING! These produce large matrices - and are a little slow

m = "p"
all.corrs =  cor( t(X),t(Y), method=m)
all.corrs.rank = apply(all.corrs, 1, rank, na.last="keep", ties.method="average")
rhoXY = diag(all.corrs)
rhoXY.ranks = diag(all.corrs.rank) /  max(all.corrs.rank, na.rm=T)
M =  max(all.corrs.rank, na.rm=T)
#save(rhoXY , rhoXY.ranks, all.corrs, all.corrs.rank, M, file="rhoXY.pearson.Rdata")

m = "s"
all.corrs =  cor( t(X),t(Y), method=m)
all.corrs.rank = apply(all.corrs, 1, rank, na.last="keep", ties.method="average")
rhoXY = diag(all.corrs)
rhoXY.ranks = diag(all.corrs.rank) /  max(all.corrs.rank, na.rm=T)
M =  max(all.corrs.rank, na.rm=T)
#save(rhoXY , rhoXY.ranks, all.corrs, all.corrs.rank, M, file="rhoXY.spearman.Rdata")


hxY = hist(as.numeric(rhoXY.ranks[,1]), plot=F, breaks=100)
frac = ( max(cumsum(hxY$counts)) - cumsum(hxY$counts)) / max(cumsum(hxY$counts))

# Plot fractions
plot(rhoXY.ranks[,1],  (rhoXY[,1]), bty="n", pch=19, col=makeTransparent("blue",5), xlab="Normalized rank of correlation coefficients between replicates", ylab="Normalized rank of correlation coefficients")
par(new=T)
plot(  hxY$breaks[-1],frac, xlab="", ylab="", lwd=2, type="l", bty="n", axes=F)
axis(4 )
mtext("Fraction of correlations tied or greater to self-replicate", side=4, padj=3)
abline(v=1, lty=3, lwd=3, col=makeTransparent(1))
abline(h=sum(rhoXY.ranks == 1, na.rm=T)/ max(cumsum(hxY$counts)), lwd=3, lty=3, col=makeTransparent(1))

##############################################################################################
# Figure 5 : MAQC mapping
##############################################################################################

# Choose "replicability" threshold
p = (rhoXY < 0.9)

# Calculate means, SDs, etc.
m.X = rowMeans(X)
m.Y = rowMeans(Y)

sd.X = apply(X,1,sd)
sd.Y = apply(Y,1,sd)

# Very coarsely split conditions by treatment versus not treated
fc1 = cols[,7] == "None"
fc2 = cols[,7] != "None"

# Calculate fold change for the X replicate
X.fc1 = X[,fc1[A][filterout]]
X.fc2 = X[,fc2[A][filterout]]
m.X.fc1 = rowMeans(X.fc1)
m.X.fc2 = rowMeans(X.fc2)

fc = abs( log2( m.X.fc2/m.X.fc1) )

# MAQC filters
minX = sort(m.X)[dim(X)[1] / 3 ]
fc[fc==-Inf] = NA
fc[fc==Inf] = NA
filt.fc = fc > 2
filt.exp = m.X > minX


# Plot MAQC comparison
plot(  log2(m.X), fc, pch=19, cex=0.1, xlab="Average expression - log2", ylab="Fold change - log2", col=makeTransparent(1))
abline(h = 1, col=4, lty=2, lwd=3)
abline(h = 1, col=4, lty=2, lwd=3)
abline(v = log2(minX), col=4, lty=3, lwd=3)

# Plot with replicability
points( log2(m.X)[ p], (fc[p]), col=makeTransparent(2), cex=0.1, pch=19)

# Plot fractions
b1 = c(-15:15)
b2 = c(0:20)

h1 = get_counts(hist(log2(m.X)[p], plot=F, breaks=b1 ))
h2 = get_counts(hist(log2(m.X)[!p] , plot=F, breaks=b1 ))
h3 = get_counts(hist(fc[p], plot=F, breaks=b2))
h4 = get_counts(hist(fc[!p] , plot=F, breaks=b2 ))

h5 = get_counts(hist(log2(m.X)[filt.fc&p], plot=F, breaks=b1 ))
h6 = get_counts(hist(log2(m.X)[filt.fc&!p] , plot=F, breaks=b1 ))
h7 = get_counts(hist(fc[filt.exp&p], plot=F, breaks=b2 ))
h8 = get_counts(hist(fc[filt.exp&!p] , plot=F, breaks=b2 ))

x  = (h1/(h1 + h2))[,2]
y  = (h3/(h3 + h4))[,2]
x2 = (h5/(h5 + h6))[,2]
y2 = (h7/(h7 + h8))[,2]

x[is.na(x)] = 0
y[is.na(y)] = 0
x2[is.na(x2)] = 0
y2[is.na(y2)] = 0

bx = h1[,1]
by = h3[,1]

# Plot mean expression density
plot(bx, x,lwd=3,type="l", , col=1,  xlab="Average expression - log2", ylab="Fraction poorly replicable",  ylim=range(x,x2), xlim=range(bx), bty="n")
lines(bx, x2, lwd=3, col="lightgrey")
abline(v = log2(minX), col=4, lty=3, lwd=3)

# Plot FC density
plot( y, by, lwd=3,type="l", ylab="Fold change - log2", xlab="Fraction poorly replicable", col=1, xlim=range(y,y2), ylim=range(by) , bty="n")
lines(y2, by, lwd=3,col="lightgrey")
abline(h = 2, col=4, lty=2, lwd=3)
abline(h = 1, col=4, lty=2, lwd=3)


# Calculate mean correlations
hm = tapply( rhoXY,  round(log2(m.X)), mean, na.rm=T)
hfc = tapply( rhoXY,  round(fc), mean, na.rm=T)
hm.fc = tapply( rhoXY[filt.fc], round(log2(m.X)[filt.fc]),  mean, na.rm=T)
hfc.m = tapply( rhoXY[filt.exp], round(fc[filt.exp]), mean, na.rm=T)

# Remove NaNs
hfc = hfc[-length(hfc)]
hm = hm[-1]

x = sort(rep( as.numeric(rownames(hfc)) ,2))
y = matrix(rbind(hfc , hfc))
x2 = sort(rep( as.numeric(rownames(hfc.m)) ,2))
y2 = matrix(rbind(hfc.m , hfc.m))
ymin=min(y, y2, na.rm=T)

# Plot mean correlation density, across FC
plot( c(x, x[length(x)]) , c(ymin,y[-length(y)],ymin), lwd=3, type="l", ylab="Mean correlations", xlab="FC - log2", pch=19, ylim=c(ymin,1) )
lines( c(x2, x2[length(x2)]) , c(ymin,y2[-length(y2)],ymin), lwd=3, col="lightgrey")
abline(v=2, col=4, lty=3, lwd=3)
abline(v=1, col=4, lty=3, lwd=3)


x = sort(rep( as.numeric(rownames(hm)) ,2))
y = matrix(rbind(hm , hm))
x2 = sort(rep( as.numeric(rownames(hm.fc)) ,2))
y2 = matrix(rbind(hm.fc , hm.fc))
ymin=min(y, y2, na.rm=T)

# Plot mean correlation density, across mean expression
plot( c(x, x[length(x)]) , c(ymin,y[-length(y)],ymin), lwd=3, type="l", ylab="Mean correlations", xlab="Mean expression - log2", pch=19, ylim=c(ymin,1) )
lines( c(x2, x2[length(x2)]) , c(ymin,y2[-length(y2)],ymin), lwd=3, col="lightgrey")
abline(v=log2(minX), col=4, lty=3, lwd=3)



##############################################################################################
#  Figure: AuPairWise on sample of RNA-seq experiments
##############################################################################################

load(file="H:/AuPairWise/suppl/gemma_rnaseq_aurocs.Rdata")

# Calculate predicted noise factor for AUROC shifts
nn = 100
i = 1
k = 1
stats.list = list()
for ( exp in unique(finalQC[,1])){
        stats = matrix(NA, ncol=nn, nrow=1 )
	filtE = finalQC[,1] == exp
	aurocs =  as.numeric(finalQC[filtE,4])
        noises =  as.numeric(finalQC[filtE,3])

	predictions = approx( noises, aurocs, n=nn)

	for (j in (nn/2):nn){
	        AUROC = j/nn
		max.pred = max(which(predictions$y <= AUROC))
	        min.pred = min(which(predictions$y > AUROC))
	        n.pred = get_value_x( predictions$x[max.pred], predictions$x[min.pred], predictions$y[max.pred], predictions$y[min.pred], AUROC)
   	        stats[i,j] = n.pred
	}
	stats.list[[k]] = stats
	k = k + 1
}

preds.all = matrix(unlist(stats.list), ncol=100, byrow=T)
colnames(preds.all) =  1:nn /nn
sizes = unique(finalQC[,1:2])

# For a given AUROC, plot the density of noise factors
AUROC.default = 0.8
col.def = which(colnames(preds.all)==AUROC.default)
h.p = get_density( hist(preds.all[,col.def], plot=F) )
cc = cor.test( preds.all[,col.def],  log10(as.numeric(sizes[,2])) , use="pair", method="s")
plot( h.p, type="l", lwd=3, xlab ="Noise factor", ylab="Freq", bty="n", main=AUROC.default )
legend("topright", legend=c( round(cc$estimate,2), format(round(cc$p.val,6), scientific=T)) ,pch=19, col=0)
par(new=T)
plot(preds.all[,col.def],  log10(as.numeric(sizes[,2])), axes=F, pch= 19, col=makeTransparent(1), xlab="", ylab="")
axis(4, at=(0:3), lab = 10^(0:3))


# Plot the pvalues of a particular noise factor
# eg.
filtN = finalQC[,3] == "25"
plot( log10(as.numeric(finalQC[filtN,2])),  -log10( as.numeric(finalQC[filtN,6])), pch=19, xlab ="Sample size", ylab="-log 10 pval", axes=F)
axis(1, at=(0:3), lab = 10^(0:3))
axis(2)




##############################################################################################
#  Figure: Brainspan RNA-seq and microarray results
##############################################################################################
colors = makeTransparent(colorpanel(15, "red", "blue"), 150))


file="Y:/human/RNAseq/brainspan/qc_results_noise_update/brainspan.rna.sub.summary.Rdata")
load(file)

# ROC examples
# eg.
ns = 5
i = 1
j = 1
k = 3
        plot(c(0,1), c(0,1),  col="grey", type="l", lwd=2, xlim=c(0,1), ylim=c(0,1),xlab="FPR", ylab="TPR", bty="n")
	for( n in 1:length(n.factors))  {
                lines(results.all[[ns]][[i]]$rocs[[j]][[n]][[k]] ,lwd=3,col=(colors[n]))
        }
        legend( "bottomright", legend= n.factors, col=colors, lwd=3, bty="n")

# Average AUROCS across sample sizes
# eg.
k = 3
aurocs.ns = t(sapply( 1:length(nSs), function(ns) aurocs.all[[k]][ns,]))
aurocs.ns.se = t(sapply( 1:length(nSs), function(ns) aurocs.se[[k]][ns,]))
aurocs.ns.sd = t(sapply( 1:length(nSs), function(ns) aurocs.sd[[k]][ns,]))

        plot( 0,0 , col=0, pch=19, bty="n", xlim=range(log10(nSs)), ylim=c(0.3,1),axes=F,xlab="Sample size", ylab="AUROC")
	axis(2)
        axis(1, lab=nSs, at=log10(nSs) )
        for ( i in 1:length(n.factors) ){
               	points( log10(nSs), aurocs.ns[,i] , col=colors[i], pch=19)
		segments(log10(nSs), aurocs.ns[,i] - aurocs.ns.se[,i], log10(nSs), aurocs.ns[,i]+aurocs.ns.se[,i],col=colors[i])
		fit = glm( aurocs.ns[,i] ~ log10(nSs), family=binomial(logit) )
                lines(  log10(nSs), fit$fitted, lwd=2, col=makeTransparent(colors[i]))
        }

	legend( "topright", legend= n.factors, col=colors, lwd=2, bty="n")

# Estimated noise for each AUROC shift
# eg.
k = 3
aurocs.n = t(sapply( 1:length(n.factors), function(n) aurocs.all[[k]][,n]))
aurocs.n.se = t(sapply( 1:length(n.factors), function(n) aurocs.se[[k]][,n]))
aurocs.n.sd = t(sapply( 1:length(n.factors), function(n) aurocs.sd[[k]][,n]))

j=1
AUROCs= c(0.55,0.6,0.7,0.8,0.9 )
n.preds = matrix(0, ncol = length(AUROCs), nrow = length(nSs) )
for (AUROC in AUROCs){
	for( i in 1:length(nSs)) {
	        data = data.frame( n.factors = n.factors, aurocs = aurocs.n[,i], aurocs.se=aurocs.n.se[,i], aurocs.sd=aurocs.n.sd[,i] )
	        fit = glm( aurocs ~ n.factors, data=data, family=binomial(logit))
	        predictions = predict(fit, data.frame(n.factors = 1:100),type = "response")
	        max.pred = max(which(predictions < AUROC))
	        min.pred = min(which(predictions > AUROC))
	        n.pred = get_value_x( max.pred, min.pred, predictions[max.pred], predictions[min.pred], AUROC)
	        n.preds[i,j] = c( n.pred)
	}
	j = j+1
}

plot( 0, 0, ylim=range(log10(n.factors[-1])),xlim=range(log10(nSs)),type="l", lwd=3, col="lightgrey", xlab="Sample size", ylab="Noise", axes=F)
axis(2,  lab=n.factors, at=log10(n.factors) )
axis(1, lab=nSs, at=log10(nSs) )
for(j in 1:5 ){
	lines(log10(nSs), log10(n.preds[,j]), lwd=3, col=colors[j] )
}
legend( "topright", legend=AUROCs, col=colors, lwd=2, bty="n")



# Plot RNA-seq brainspan results
plot_results("H:/AuPairWise/suppl/brainspan/brainspan.rna.summary.Rdata")
plot_distributions("H:/AuPairWise/suppl/brainspan/brainspan.rna.summary.Rdata")

# Get table of pvalues
load("H:/AuPairWise/suppl/brainspan/brainspan.rna.summary.Rdata")
pvals.rnseq = t(sapply(1:length(n.factors),  function(i) c(  wilcox.test( results.all[[i]]$aurocs[[1]], results.all[[i]]$aurocs[[2]])$p.val, mean(results.all[[i]]$aurocs[[1]]) , mean(results.all[[i]]$aurocs[[2]]) )))


# Plot microarray brainspan results
plot_results("H:/AuPairWise/suppl/brainspan/brainspan.micr.summary.Rdata")
plot_distributions("H:/AuPairWise/suppl/brainspan/brainspan.micr.summary.Rdata")

# Get table of pvalues
load("H:/AuPairWise/suppl/brainspan/brainspan.micr.summary.Rdata")
pvals.micr = t(sapply(1:length(n.factors),  function(i) c(  wilcox.test( results.all[[i]]$aurocs[[1]], results.all[[i]]$aurocs[[2]])$p.val, mean(results.all[[i]]$aurocs[[1]]) , mean(results.all[[i]]$aurocs[[2]]) )))



##############################################################################################
#  Figure: ENCODE results
##############################################################################################

# RPKM
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.rep.summary.Rdata")
# TPM
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.normalized.summary.Rdata")
# TMM
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.TMM.summary.Rdata")
# VST
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.VST.summary.Rdata")
# Ranked
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.ranked.summary.Rdata")
# Randomized
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.randomized.summary.Rdata")
# Replicates
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.XY.rep.summary.Rdata")
# Replicates randomized
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.between_reps_rseq.XY.rep.randomized.summary.Rdata")

# eg.
plot_results(file)
plot_distributions(file)


# Centiles results
file=("H:/AuPairWise/suppl/ENCODE/ENCODE.rseq.post.sample.18.all.aurocs.percentiles.spearman.abs.setsize.2665.Rdata")
load(file)
indices.nfactors  = 1:length(n.factors.sub)
indices.nfactors  = c(1,6)

plot( R,percentiles[R,1], pch=19, ylim=c(0.4,1), col=0, xlab="Centile", ylab="Average AUROCs" ) #main=paste(methods[im], "samples", nS[is])
for( j in indices.nfactors ){
    points( R, percentiles[,j], pch=19, col=makeTransparent(colors[j]))
    #smooth = convolve_nets(R, percentiles[,j],20)
    #lines(smooth,col=makeTransparent(colors[j]))
    segments( R, percentiles[,j]-percentiles.se[,j], R, percentiles[,j]+percentiles.se[,j],col=makeTransparent(colors[j]))
    #abline( h = mean( percentiles[,j]), col=makeTransparent(colors[j]))

}


##############################################################################################
#  Functions: plotting results of AuPairWise
##############################################################################################

plot_results <- function(file,J=2){
        load(file)
        range = n.factors[-1]
        plot( log10(n.factors), aurocs.all[1,], pch= 19, col=0, ylim=c(0.4,1 ), xlab="Noise factor", ylab="AUROC", axes=F)
        axis(2)
        axis(1, lab=range, at=log10(range) )
        for ( j in 1:J){
            points(log10(n.factors), aurocs.all[j,], pch=19, col=colors[j])
            lines(log10(n.factors), aurocs.all[j,], lwd=2, col=colors[j])
            segments( log10(n.factors), aurocs.all[j,]-aurocs.se[j,], log10(n.factors),aurocs.all[j,]+aurocs.se[j,],col=colors[j])
        }
}

plot_result <- function(file, j=1){
        load(file)
        range = n.factors[-1]
        plot( log10(n.factors), aurocs.all, pch= 19, col=0, ylim=c(0.4,1 ), xlab="Noise factor", ylab="AUROC", axes=F)
        axis(2)
        axis(1, lab=range, at=log10(range) )
        points(log10(n.factors), aurocs.all, pch=19, col=colors[j])
        lines(log10(n.factors), aurocs.all, lwd=2, col=colors[j])
        segments( log10(n.factors), aurocs.all-aurocs.se, log10(n.factors),aurocs.all+aurocs.se,col=colors[j])
}


plot_distributions <- function(file){
        load(file)
	for( j in 1:length(n.factors) ) {
		aurocs1 = results.all[[j]]$aurocs[[1]]
    		aurocs2 = results.all[[j]]$aurocs[[2]]
    		h2 = get_density(hist(aurocs2, plot=F))
		h1 = get_density(hist(aurocs1, plot=F))
		m2 = round(mean(aurocs2),2)
		m1 = round(mean(aurocs1),2)
		pval = wilcox.test( aurocs2,aurocs1 )$p.val
		    
		plot(h1, type="l", lwd=3, xlim=c(0,1), sub = pval, main= n.factors[j], bty="n", xlab="Average AUROCs", ylab="Freq/density")
		polygon(h1, col=makeTransparent(1))
		polygon(h2, col=makeTransparent("grey"))
		lines(h2, lwd=3, col="lightgrey")
		abline(v=m1, lwd=3, lty=2)
		abline(v=m2, lwd=3,lty=2, col="lightgrey")
		legend("topleft", legend = c(m1,m2), col=c(1,"lightgrey"), lwd=3, lty=2)
	}
}


