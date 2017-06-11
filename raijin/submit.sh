#!/bin/bash
#PBS -P xe2
#PBS -q express
#PBS -l ncpus=16
#PBS -l walltime=24:00:00
#PBS -l other=gdata1
#PBS -l mem=31G
#PBS -l jobfs=100G
#PBS -l wd

logdir=raijin/log
if [ -d $logdir ]
then
  pushd $logdir >/dev/null
  if [ -n "$(ls *.ER *.OU 2>/dev/null)" ]
  then
    tar cf `date +%y%m%d_%H%M%S`.tar *.OU *.ER 2>/dev/null
    rm *.OU *.ER
  fi
  popd >/dev/null
else
  mkdir -p $logdir
fi

. <(grep module raijin/jobscript.sh)

QSUB="qsub -q {cluster.queue} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem}"
QSUB="$QSUB -l wd -o $logdir -e $logdir -P xe2"

snakemake                                \
    -j 3000                              \
    --cluster-config raijin/cluster.yaml \
    --js raijin/jobscript.sh             \
    --local-cores 16                     \
    --rerun-incomplete                   \
    --keep-going                         \
    --cluster "$QSUB" $@ >raijin/log/snakemake.log
