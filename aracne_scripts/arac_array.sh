#!/usr/bin/env bash

#SBATCH --time=08:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --job-name=aracne-array
#SBATCH --output=/home/exacloud/lustre1/CompBio/manningh/slurm/log-files/aracne-array_%j.out
#SBATCH --error=/home/exacloud/lustre1/CompBio/manningh/slurm/log-files/aracne-array_%j.err


cd /home/exacloud/lustre1/CompBio/manningh/PyCharmProjects/smmartVIPER/viper_scripts/
source activate py351Crayon

# ALLOW COMMAND LINE ARGS
OPTIND=1

while getopts "a:e:o:t:" opt; do
    case "$opt" in
    a)  ARACNE=$OPTARG;;
    e)  EFILE=$OPTARG;;
    o)  OUTDIR=$OPTARG;;
    t)  TFFILE=$OPTARG;;
    esac
done

shift $((OPTIND-1))

# print slurm array task ID
#echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

## run X bootstraps
# Note that with the Tatlow TCGA BRCA, each bootstrap takes ~1 hour
java -Xmx5G -jar ${ARACNE} -e ${EFILE} -o ${OUTDIR} --tfs ${TFFILE} --pvalue 1E-8