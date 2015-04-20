RUN_WALKER

##########################################################
##########################################################

## Names/Paths for Output
DATE=$1 # Today's date plus an option identifier
HOME_DIR=$2 # Directory all sub-folders will be created in

## Parameters
SNP_OR_BASE=$3 # "SNP" or "BASE": Bin regions by number of SNPs or number of Bases?
NUM_UNITS=$4 # Number of Units (defined above) to be pooled into region
WHICH_TESTS=$5 # List of Test ID's that will trigger running a given test (e.g., "BURD","SKAT","LASSO","ELNET",etc...)

## Files
VAR_FILE=$6
PHENO_FILE=$7 # Which Phenotype File are you using?
PHENO_TYPE=$8 # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=$9 # Path to Covariate File or "F"
COVS=${10} # Which Covariates to Include?
EIG_VEC=${11} # Output from Plink's --pca command (MAF>1%) or "F"
PC_COUNT=${12} # How many PCs to Include as Covariates?
START_STEP=${13} # Which Step do you want to start on?

##########################################################
## SPECIFY INPUT PARAMETERS ##############################
#### CHR22, bin 1000 SNPs ################################

## Names/Paths for Output
DATE=20150408
HOME_DIR=/projects/janssen/Walker

## Parameters
SNP_OR_BASE=SNP
NUM_UNITS=1000
WHICH_TESTS=`echo BURD,ELNET`

## Files
VAR_FILE=/projects/janssen/Walker/TEMP_CHR22.bed

PHENO_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt
PHENO_TYPE=C
# COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20141229_Full_Table.txt
COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo DAS_BL_MN SEX`
EIG_VEC=/projects/janssen/ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=2
START_STEP=6
CLEAN_UP=F

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

##########################################################
## SPECIFY INPUT PARAMETERS ##############################
#### CHR22, bin 1 Megabase ###############################

## Names/Paths for Output
DATE=20150409
HOME_DIR=/projects/janssen/Walker

## Parameters
SNP_OR_BASE=BASE
NUM_UNITS=1000000
WHICH_TESTS=`echo BURD,ELNET`

## Files
VAR_FILE=/projects/janssen/Walker/TEMP_CHR22.bed

PHENO_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt
PHENO_TYPE=C
# COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20141229_Full_Table.txt
COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo DAS_BL_MN SEX`
EIG_VEC=/projects/janssen/ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=2
START_STEP=1
CLEAN_UP=F

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`
##########################################################
## RUN COMMAND ###########################################

/projects/janssen/Walker/SCRIPTS/Walker.sh \
${DATE} \
${HOME_DIR} \
${SNP_OR_BASE} \
${NUM_UNITS} \
${WHICH_TESTS} \
${VAR_FILE} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT} \
${START_STEP} \
${CLEAN_UP}









