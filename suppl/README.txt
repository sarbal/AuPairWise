
###################################################################
###################################################################
###################################################################

corrs-examples.Rdata:  Data for correlation figures
Variables in file:
- exprs: SRP000225 - Cufflinks exprsession estimation
- samples: sample names for SRP000225 experiment
- rsem: GSE26109 - RSEM exprsession estimation
- atlas: RNAseq Atlas
- all.corrs: pearson correlation coefficients between the samples
- f1, f2: match between exprs and atlas ( exprs[f1,] vs atlas[f2,] )
- f3, f4: match between exprs and rsem  ( exprs[f3,] vs rsem[f4,] )
- f5, f6: match between rsem and atlas  ( rsem[f5,] vs atlas[f6,] )

# To plot panels of figure in paper:
i = 1
j = 4
pdf("corrs-examples.pdf")
plot( log2(exprs[,i]), log2(exprs[,i+1]),cex=0.5, main=round(cor(exprs[,i], exprs[,i+1], use="pair", method="p"),3),sub=round(cor(exprs[,i], exprs[,i+1], use="pair", method="s"),3),
xlab=paste(samples[i], "Liver mRNA - log2 FPKM") , ylab=paste(samples[i+1], "Liver mRNA - log2 FPKM"), bty="n", pch=19, col=makeTransparent(1))
plot( log2(exprs[,i]), log2(exprs[,j]),cex=0.5, main=round(cor(exprs[,i], exprs[,j], use="pair", method="p"),3),sub=round(cor(exprs[,i], exprs[,j], use="pair", method="s"),3),
xlab=paste(samples[i], "Liver mRNA - log2 FPKM") , ylab=paste(samples[j], "Kidney mRNA - log2 FPKM"), bty="n", pch=19, col=makeTransparent(1))
plot( log2(exprs[f3,i]), log2(rsem[f4,i]), cex=0.5,main=round(cor(exprs[f3,i], rsem[f4,i], use="pair", method="p"),3),sub=round(cor(exprs[f3,i], rsem[f4,i], use="pair", method="s"),3),
xlab=paste(samples[i],"Liver mRNA - log2 FPKM"), ylab=paste(colnames(rsem)[i], "Liver - log2 RPKM"), bty="n", pch=19, col=makeTransparent(1))
plot( log2(exprs[f1,i]), log2(atlas$liver[f2]), cex=0.5,main=round(cor(exprs[f1,i], atlas$liver[f2], use="pair", method="p"),3),sub=round(cor(exprs[f1,i], atlas$liver[f2], use="pair", method="s"),3) ,xlab=paste(samples[i], "Liver mRNA - log2 FPKM"), ylab="RNA-seq atlas, liver - log2 RPKM", bty="n", pch=19, col=makeTransparent(1))
dev.off()

