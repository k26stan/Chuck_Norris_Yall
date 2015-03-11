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

# ANNOTS=$4 # Path to Annotation File
PHENO_FILE=$5 # Which Phenotype File are you using?
PHENO_TYPE=$6 # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=$7 # Path to Covariate File or "F"
COVS=$8 # Which Covariates to Include?
EIG_VEC=$9 # Output from Plink's --pca command (MAF>1%) or "F"
PC_COUNT=${10} # How many PCs to Include as Covariates?
START_STEP=${11} # Which Step do you want to start on?

##########################################################
##########################################################

DATE=20150310a
HOME_DIR=/projects/janssen/Walker

## Parameters
SNP_OR_BASE=SNP
NUM_UNITS=1000
WHICH_TEST=`echo BURD, ELNET`

## Files
VAR_FILE=/projects/janssen/Walker/TEMP_BED.bed

PHENO_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt
PHENO_TYPE=C
# COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20141229_Full_Table.txt
COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo DAS_BL_MN AGE SEX`
EIG_VEC=F
PC_COUNT=2
START_STEP=1

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`











