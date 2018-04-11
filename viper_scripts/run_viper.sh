#!/usr/bin/env bash

#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=20000
#SBATCH --job-name=smmart_viper
#SBATCH --output=/home/exacloud/lustre1/CompBio/manningh/slurm/log-files/smmart_viper_%j.out
#SBATCH --error=/home/exacloud/lustre1/CompBio/manningh/slurm/log-files/smmart_viper_%j.err
#SBATCH --verbose


# STEP1: GET ORIENTED
cd /home/exacloud/lustre1/CompBio/manningh/PyCharmProjects/smmartVIPER/viper_scripts/
source activate py351Crayon

# STEP2: ALLOW COMMAND LINE ARGS
OPTIND=1

# note that COHORT is the TCGA cohort background, regardless of study
while getopts "c:s:u:i:" opt; do
    case "$opt" in
    c)  COHORT=$OPTARG;;    # (i.e. BRCA)
    s)  STUDY=$OPTARG;;     # (i.e. TCGA)
    u)  UHR=$OPTARG;;       # (i.e. 111111_NN555555_0000_AAAA1AAAA1_RNA111111AA_UHR)
    i)  SOI=$OPTARG;;       # sample of interest (i.e. 111111_NN55555_0000_AAAA1AAAA1_RNA111111AA_DNA.11.00000)
    esac
done

shift $((OPTIND-1))

# STEP3: SPECIFY VARIABLES
DATE=`date +%Y%m%d`

# specify dirs
BASE_DIR=/home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/
OUTDIR=${BASE_DIR}VIPER_results/${STUDY}_${COHORT}_${SOI}/

# get context specific ARACNE network and the expression used to generate it
NETWORK=${BASE_DIR}ARACNE_nets/${STUDY}_${COHORT}/network.txt
ORIG_EXPR=/home/exacloud/lustre1/CompBio/tatlow_matrices/${STUDY}_${COHORT}_bygene.txt.gz
EXPR_ARAC=${BASE_DIR}ARACNE_nets/${STUDY}_${COHORT}/${STUDY}_${COHORT}_bygene_prepped_ensembl.txt

# if that expression file is gzipped, make note of that
if [[ ! -f ${EXPR_ARAC} ]]; then
    EXPR_ARAC=${BASE_DIR}ARACNE_nets/${STUDY}_${COHORT}/${STUDY}_${COHORT}_bygene_prepped_ensembl.txt.gz
fi

# specify location of new expression matrix (background expr + smmart sample + UHR, normalized)
NEW_EXPR=${OUTDIR}pre_viper_expr.tsv

# if outdir doesn't exist, create it
if [[ ! -e ${OUTDIR} ]]; then
    echo "Making outdir"
    mkdir ${OUTDIR}
fi

# write README.txt
READFL=${OUTDIR}README.txt

RDMSGS=("Date of run: ${DATE}"
"R dependencies may be viewed in RsessionInfo.txt in this dir."
"COHORT (for background expr) is ${COHORT}."
"STUDY of sample of interest is ${STUDY}."
"ARACNE NETWORK is ${NETWORK}; regulon is derived from this file."
"ARACNE NETWORK was derived from ${EXPR_ARAC}. A list of the regulators used is in that dir."
"Expression network given to VIPER is pre_viper.tsv in this dir (includes SOI and UHR if present)")

touch ${READFL}
for ((i = 0; i < ${#RDMSGS[@]}; i++))
do
    echo "${RDMSGS[$i]}" >> ${READFL}
done

# STEP4: PREPROCESS DATA (Will add specified samples to tatlow matrix, normalize, and print expr into outdir)
# if NEITHER sample of interest nor uhr have been set...
echo "Preparing to run preprocess.py"
if [[ -z ${SOI} ]] && [[ -z ${UHR} ]]; then
    echo "Neither Sample of Interest nor UHR provided" >> ${READFL}
    python preprocess.py -b ${ORIG_EXPR}
# if sample of interest has been set but no UHR is provided...
elif [[ -z ${UHR} ]]; then
    echo "Sample of Interest is ${SOI}. UHR not provided." >> ${READFL}
    python preprocess.py -b ${ORIG_EXPR} -n ${NEW_EXPR} -s ${SOI}
# if BOTH sample of interest and UHR have been set...
else
    echo "Sample of Interest is ${SOI}" >> ${READFL}
    echo "UHR is ${UHR}" >> ${READFL}
    python preprocess.py -b ${ORIG_EXPR} -n ${NEW_EXPR} -s ${SOI} -u ${UHR}
fi

# STEP5: RUN VIPER (convert network to regulon, use expr without pdata instead of xset)
echo "Preparing to run_viper.R"
Rscript run_viper.R -d ${OUTDIR} -e ${EXPR_ARAC} -n ${NETWORK} -s ${NEW_EXPR}
# Rscript run_viper.R -d /home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/VIPER_results/TCGA_BRCA_170413_NS500681_0121_AHJC5WBGX2_RNA170210CC_DNA.17.00536/ -e /home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/ARACNE_nets/TCGA_BRCA/TCGA_BRCA_bygene_prepped_ensembl.txt.gz -n /home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/ARACNE_nets/TCGA_BRCA/network.txt -s /home/exacloud/lustre1/share_your_data_here/precepts/aracne_viper_analyses/kallisto/VIPER_results/TCGA_BRCA_170413_NS500681_0121_AHJC5WBGX2_RNA170210CC_DNA.17.00536/pre_viper_expr.tsv.gz


# STEP6: POSTPROCESS DATA (MAP IDs TO HUGO; PROVIDE THAT MAP)
python postprocess.py -d ${OUTDIR}


# STEP7: CLEAN UP
if [ ${NEW_EXPR} ]; then
    gzip ${NEW_EXPR}
fi
