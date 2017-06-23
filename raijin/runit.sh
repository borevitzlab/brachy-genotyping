#!/bin/bash
#PBS -P xe2
#PBS -q express
#PBS -l ncpus=16
#PBS -l walltime=24:00:00
#PBS -l other=gdata1
#PBS -l mem=63G
#PBS -l jobfs=100G
#PBS -l wd

. <(grep module raijin/jobscript.sh)

snakemake --unlock

snakemake                        \
    -j 16                        \
    --rerun-incomplete           \
    --keep-going                 \
    varcall                      \
    >data/log/snakemake.log 2>&1 \

