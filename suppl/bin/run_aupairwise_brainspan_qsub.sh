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

in='/human/brainspan/qc_results_extra'
file='brainspan'
j=$1
flag=$2

$R --no-save --args $dir $in $bin $file $j $flag < $bin/run_aupairwise_brainspan.r