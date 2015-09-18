#!/bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
#  specify that the job requires 20GB of memory
#$ -l m_mem_free=20G
#$ -v OMP_NUM_THREADS

export OMP_NUM_THREADS=16

R='/opt/hpc/pkg/R-2.15/bin/R'
dir='/sonas-hs/gillis/hpc/data/sballouz/'
bin='/sonas-hs/gillis/hpc/home/sballouz/'

in='/human/RNAseq/ENCODE/qc'
file='ENCODE.between_reps_rseq'
j=$1
flag=$2
rank=$3
rand=$4
nfold=$5

$R --no-save --args $dir $in $bin $file $j $flag $rank $rand $nfold < $bin/run_aupairwise_encode_qsub.r