##############################################################################################
# Supplementary data and code for AuPairWise paper



##############################################################################################
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


##############################################################################################
gemma_rnaseq_aurocs.Rdata:  Results of RNA-seq experiments from GEMMA
Variables in file:
- resultsQC.list
- finalQC


##############################################################################################
GEO_series_2014.Rdata: GSE ids and sample size of SRA experiments in GEO dated until 2014
Variables in file:
- all_exp_parsed


##############################################################################################
AuPairWise_analysis.r
Code to draw figures in the paper



##############################################################################################
bin/
..
example_qsub.txt
helper_aupair_develp.r
run_aupairwise_brainspan_qsub.r
run_aupairwise_brainspan_qsub.sh
run_aupairwise_encode_qsub.r
run_aupairwise_encode_qsub.sh


Scripts to run ENCODE and brainspan analyses of the paper

