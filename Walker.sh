## Write Script to Walk Through Genome & Test Regions ##
## February 6, 2015 ##
## Kristopher Standish ##

######################################################
## OVERVIEW/GOALS ####################################
######################################################

## Steps ##

# Specify whether grouping is by number of Bases or number of SNPs
# Specify number of Bases/SNPs (X)
# Loop through Variant File, pulling out X SNPs at a time
  # Save to variant ID file
# Select Tests to be run on each Region
  # Run gamut of tests
  # Collect test results (p-value, AUC, R2...some measure of quality)
# Compile results for each region
# Plot this Shiz

##########################################################################
## 1 ## Set up Paths #####################################################
##########################################################################
 # Use Bash
 # Take in arguments, set up directories/paths for files/tools
echo \### 1 - `date` \###
echo \### Define Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

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

###########################################################
## Constant Paths ##

## Public Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt

## Custom Scripts

###########################################################
## Pull some Info out of Parameters ##

## Get Names of Specific Files
TEMP=(${VAR_FILE//\// })
VAR_FILE_NAME=${TEMP[${#TEMP[@]} - 1]} # Get Name of Variant File
TEMP=(${PHENO_FILE//\// })
PHENO=${TEMP[${#TEMP[@]} - 1]} # Get Name of Phenotype File

## Specify list of Covariates to include (for command and for filename)
if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

## Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi

## Specify extensions for Continuous vs Binary Phenotype
if [ $PHENO_TYPE = "C" ]
then
SUFFIX=linear
else
SUFFIX=logistic
fi

## Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
echo $OUT_DIR
OUT_DIR=${HOME_DIR}/${DATE}_${PHENO%%.txt}_${COVS_FILENAME}
mkdir ${OUT_DIR}
cd ${OUT_DIR}

## Specify Output File Paths & Names
 # Maybe

## Specify a File to which to Write Updates
UPDATE_FILE=${OUT_DIR}/Update.txt

## Done
if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Define Set Variables and Paths - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 2 ## Determine Variant Format and Adjust ##############################
##########################################################################
 # Use Bash
 # Convert Variant File into BED format
   # If File is VCF, make a PED/MAP file, then BED file
   # If File is PED, make BED file
   # If File is BED, skip ahead
if [ "$START_STEP" -le 2 ]; then
echo \### 2 - `date` \###
echo \### Determine/Adjust Variant File Formats \###
echo `date` "2 - Determine/Adjust Variant File Formats" >> ${UPDATE_FILE}

## Determing File Type and Convert to .bed File
if [ ${VAR_FILE: -4} == ".vcf" ] ; then
        ${VCF_TOOLS} --plink --vcf ${VAR_FILE} --out ${OUT_DIR}_${VAR_FILE_NAME%%.vcf}
        VAR_FILE=${OUT_DIR}_${VAR_FILE_NAME%%.vcf}.ped
fi
if [ ${VAR_FILE: -4} == ".ped" ] ; then
        ${PLINK} --make-bed --file ${VAR_FILE%%.ped} --out ${VAR_FILE%%.ped}
        VAR_FILE=${VAR_FILE%%.ped}.bed
fi

## Done
echo `date` "2 - Determine/Adjust Variant File Formats - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 ## Create Variant Lists for Specified Regions #######################
##########################################################################
 # Use Bash
 # Create Directory for Coordinate Files and Write Coord Files for each Region
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Create Variant Lists \###
echo `date` "3 - Create Variant Lists" >> ${UPDATE_FILE}

## Make Directory to write New SNP lists to
OUT_DIR_3=${OUT_DIR}/SNP_Lists
mkdir OUT_DIR_3

## Specify Variant List
VAR_LIST=${VAR_FILE%%.bed}.bim

##########################################################
## If Regions are determined by number of SNPs
if [ $SNP_OR_BASE == "SNP" ]
then

## Run Python Script to Write New Variant Lists (output to $OUT_DIR_3)
python 3-Split_Var_List_SNP.py ${VAR_LIST} ${OUT_DIR_3} ${NUM_UNITS}

else
##########################################################
## If Regions are determined by number of BASES

fi










































## Done
echo `date` "12 - Filter Candidates and Build Model - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
#####################################################################
## END OF DOC #######################################################
#####################################################################