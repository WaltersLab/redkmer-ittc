#!/bin/bash
#PBS -N redk_plot5
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=12:mem=250gb:tmpspace=500gb
#PBS -e /work/nikiwind/
#PBS -o /work/nikiwind/

module purge
module load R

source $PBS_O_WORKDIR/redkmer.cfg
export R_LIBS="/home/nikiwind/localRlibs"

cd $PBS_O_WORKDIR/Rscripts/
R CMD BATCH --no-save --no-restore R_kmers_plot5.R > ${CWD}/reports/R_kmers_plot5.Rout

cp $PBS_O_WORKDIR/Rscripts/R_kmers_plot5.Rout ${CWD}/reports/
