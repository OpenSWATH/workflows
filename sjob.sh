#!/bin/bash
#$ -N SLOG
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::

snakemake --snakefile Snakefile.library -j 100 --cluster-config envs/lsf.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0"
#snakemake --snakefile Snakefile.library -j 24 --cluster-config envs/lsf.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem=32G,time=8:0:0"

