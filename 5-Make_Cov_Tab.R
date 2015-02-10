# Make Covariate File w/ EigenVectors #

# Make_Cov_Tab.R <Path/To/EigenVectors> <Path/To/Covariate_Table> <Path/To/Phenotype> <Path/To/Output>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/Walker/20150209a_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_PC1_PC2/PCs.eigenvec","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/COV.txt","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt", "/projects/janssen/Walker/20150209a_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_PC1_PC2/Cov_w_PCs.txt")
PathToVec <- LINE[1]
PathToCov <- LINE[2]
PathToPheno <- LINE[3]
PathToNewCov <- LINE[4]

print(paste( "Vec:", PathToVec ))
print(paste( "Cov:", PathToCov ))
print(paste( "Pheno:", PathToPheno ))
print(paste( "NewCov:", PathToNewCov ))

## Load Tables
VEC <- read.table(PathToVec,header=T)
COV <- read.table(PathToCov,sep="\t",header=T)
PHE <- read.table(PathToPheno, sep="\t",header=T)
print("Done Loading Files")

## Merge Tables
VEC.2 <- VEC[,c(1,3:ncol(VEC))] ; colnames(VEC.2)[1] <- "FID"
PHE.2 <- PHE[,c(1,3:ncol(PHE))] ; colnames(PHE.2)[1] <- "FID"
MRG.1 <- merge(x=COV,y=VEC.2,by="FID")
MRG.2 <- merge(x=MRG.1,y=PHE.2,by="FID")
print("Done Merging Files")

## Write Table
write.table(MRG.2,PathToNewCov,sep="\t",row.names=F,col.names=T,quote=F)

######### END OF DOC ###########
