#!/usr/bin/env bash

# IMPORTANT: this script should be run as:
#   nohup ./run_aracne.sh -allofyouroptions > ~/lustrehome/slurm/log-files/mylog.log 2>&1 &
#   echo $! > save_pid.tx
# this will run it in the background and prevent the hangup signal if you close your terminal or travel

# to kill a job, options are:
#   1.
#       kill -9 'cat save_pid.tx'
#
#   2.
#       ps aux | grep run_aracne.sh
#       kill -9 thejobIDincol2
#
# or, to kill all slurm jobs: scancel -u manningh


# set default number of reps and study, modifiable if provided to getopts by user.
REPS=20
STUDY=TCGA

while getopts "c:s:r:" opt; do
    case "$opt" in
    c)  COHORT=$OPTARG;;    # (i.e. BRCA)
    s)  STUDY=$OPTARG;;     # (i.e. TCGA)
    r)  REPS=$OPTARG;;      # number of aracne runs to be initiated and consolidated (i.e. 20)
    esac
done

shift $((OPTIND-1))

REPS=$(seq ${REPS})
ARACNE=~/lustrehome/Aracne.jar
OUTDIR=/home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/ARACNE_nets/${STUDY}_${COHORT}/
EFILE=${OUTDIR}${STUDY}_${COHORT}_bygene_prepped_ensembl.txt
TFFILE=${OUTDIR}pancancer-tfgenes_2014_04_17_filtered_${COHORT}_ensembl.txt

# tatlow expression to be formatted into ${EFILE}
TATEX=/home/exacloud/lustre1/CompBio/tatlow_matrices/${STUDY}_${COHORT}_bygene.txt.gz

cd ~/lustrehome/PyCharmProjects/smmartVIPER/aracne_scripts/
source activate py351Crayon

mkdir ${OUTDIR}

python preprocess.py -ei ${TATEX} -eo ${EFILE} -ro ${TFFILE}

# RUN ARACNE
# 1. calculate threshold with a fixed seed
# todo: INCLUDE SRUN ARGS???
# todo: why not something more like this:
# https://stackoverflow.com/questions/41900600/slurm-sbatch-job-array-for-the-same-script-but-with-different-input-arguments-ru?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
srun java -Xmx5G -jar ${ARACNE} -e ${EFILE} -o ${OUTDIR} --tfs ${TFFILE} --pvalue 1E-8 --seed 2027 --calculateThreshold

# 2. run X reproducible bootstraps using an "array" which is really a string...
# Note that with the Tatlow TCGA BRCA, each bootstrap takes ~1 hour
ARR=""
for I in ${REPS}
do
    SUBM_STATEMENT=$(sbatch arac_array.sh -a ${ARACNE} -e ${EFILE} -o ${OUTDIR} --t ${TFFILE})
    SUBM_ID=$(echo ${SUBM_STATEMENT} | awk '{print $4}')
    ARR="${ARR} ${SUBM_ID}"
done

# periodically check whether all STEP2 jobs have completed
STEP2DONE="FALSE"
while [ "${STEP2DONE}" != "TRUE" ];
do
    for J in ${ARR};
    do
        if [ "$(sacct -j ${J} -n --format=State -X -P)" = "COMPLETED" ]; then
        ARR=${ARR//${J}/}
        fi
    done
    sleep 30
    # if there's nothing left in the array besides the spacing, all jobs have completed
    if [ "${ARR//[[:blank:]]/}" = "" ]; then
        STEP2DONE="TRUE"
    fi
done

# 3. consolidate bootstraps in the output folder
# note that no lines are added to network.txt if --nobonferroni tag isn't present and the network isn't large enough (or not enough reps?)
srun java -Xmx5G -jar ${ARACNE} -o ${OUTDIR} --consolidate --nobonferroni

# zip the expression file back up
if [ ${EFILE} ]; then
    gzip ${EFILE}
fi
