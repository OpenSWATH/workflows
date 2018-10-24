#!/bin/bash
#$ -N SLOG
#$ -cwd
#$ -pe smp 1
#$ -l mem=8G,time=8::

snakemake --snakefile Snakefile.openswath -j 100 --latency-wait=30 --cluster-config envs/lsf.json --cluster "qsub -cwd -pe smp {cluster.nCPUs} -l mem={cluster.memory}M,time=8:0:0"

