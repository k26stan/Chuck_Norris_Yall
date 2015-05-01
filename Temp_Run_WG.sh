#!/bin/bash
#PBS -N 20150430_WG_1000SNP
#PBS -q hotel
#PBS -v QOS=10
#PBS -l walltime=10:00:00
#PBS -m abe
#PBS -l nodes=1:ppn=1
#PBS -o 20150430_WG_1000SNP.run.oe
#PBS -j oe
#PBS -M k26stan@yahoo.com
#PBS -V
#PBS -A schork-group

echo "<startTime>"`date`"</startTime>"
echo "<output>"

## Overall Names/Paths for Output
DATE_1=20150430
TOP_DIR=/projects/janssen/Walker/${DATE_1}_WG_1000SNP
FULL_VAR_FILE=/projects/janssen/VCFs/PLINK/BED_FULL.ALL
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink


## Make New Directory for this little Experiment
mkdir ${TOP_DIR}
cd ${TOP_DIR}

##########################################################
## Loop through Chromosomes
for chr in `seq 1`; do
echo \############## RUN CHR ${chr} \################

## Pull out Chromosome
mkdir ${TOP_DIR}/Chroms/
${PLINK} --make-bed --chr ${chr} --bfile ${FULL_VAR_FILE} --out ${TOP_DIR}/Chroms/BED_${chr}

DATE=${DATE_1}_${chr}
HOME_DIR=${TOP_DIR}

## Parameters
SNP_OR_BASE=SNP
NUM_UNITS=1000
WHICH_TESTS=`echo BURD,ELNET`

## Files
VAR_FILE=${TOP_DIR}/Chroms/BED_${chr}.bed

PHENO_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt
PHENO_TYPE=C
COV_FILE=/projects/janssen/ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo DAS_BL_MN`
EIG_VEC=/projects/janssen/ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=2
START_STEP=1
CLEAN_UP=T

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

done # Close "chr" Loop

echo "</output>"
echo "<exitStatus>"$?"</exitStatus>"
echo "<stopTime>"`date`"</stopTime>"
qstat -f $PBS_JOBID | grep Job
qstat -f $PBS_JOBID | grep resources