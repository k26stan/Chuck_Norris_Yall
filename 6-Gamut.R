## R Script to Load Genotype Data & Run Tests ##
## Will load GT data for each bin and run whichever tests are specified ##
## To be used in Walker.sh Script ##
## Feburary 9, 2015 ##
## Kristopher Standish ##

## Usage ##
# Rscript 6-Gamut.py <Path/To/SNP_Lists.txt> <Path/To/Cov_Pheno.txt> <Which_Tests>

library( glmnet )

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/Walker/20150209a_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150209a_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,AGE,SEX,PC1,PC2")
# LINE <- c("/projects/janssen/Walker/20150310a_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150310a_LT8_DEL_MNe_MN_DAS_BL_MN_AGE_SEX_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,AGE,SEX,PC1,PC2")
PathToList <- LINE[1]
PathToPheno <- LINE[2]
Which_Tests <- LINE[3]
Which_Covs <- LINE[4]

# Specify other Paths
PathToOut <- gsub( "Cov_w_PCs.txt","", PathToPheno )
PathToGT <- gsub( "SNP_Lists.txt","", PathToList )

# Check for proper parsing
print(paste( "List:", PathToList ))
print(paste( "Pheno:", PathToPheno ))
print(paste( "Tests:", Which_Tests ))
print(paste( "Covariates:", Which_Covs ))
print(paste( "Output:", PathToOut ))
print(paste( "Genotype Files:", PathToGT ))

# Split Which_Covs and Which_Tests
Which_Tests.spl <- strsplit( Which_Tests, "," )[[1]]
Which_Covs.spl <- strsplit( Which_Covs, "," )[[1]]

###############################################################
## LOAD DATA ##################################################
###############################################################

## Open List of Genotype Files
FILE_LIST <- as.character( read.table( PathToList, header=F )[,1] )

## Open Phenotype/Covariate File
PHENO_FILE <- read.table( PathToPheno, sep="\t", header=T )

###############################################################
## GET ORGANIZED ##############################################
###############################################################

## Get Sample List
SAMPS <- PHENO_FILE$FID
N_SAMPS <- length(SAMPS)

## Get Number of Files in File List
N_FILES <- length(FILE_LIST)

## Pull out Phenotype
PHE <- PHENO_FILE$Pheno ; names(PHE) <- PHENO_FILE$FID
 
## Pull out Covariates
COV <- PHENO_FILE[,Which_Covs.spl] ; rownames(COV) <- PHENO_FILE$FID

## Throw these into single data.frame
CP <- data.frame( Pheno=PHE, COV )

## Fit Phenotype vs Covariates & Calculate Residuals
RESIDS <- residuals( lm( Pheno ~ ., data=CP ) )

###############################################################
## CREATE OBJECTS TO COMPILE RESULTS ##########################
###############################################################

# ## Compile Stats on Cohort for each Region
# COH <- array( ,c(N_SAMPS,5) )
# colnames(COH) <- c("FID","Pheno","Resid","")
# COH <- data.frame( rownames(CP), )

## Single-Locus Results/Stats
 # SNP name, MAF, P-Val (linear model)
GWAS <- array( ,c(0,5) )
colnames(GWAS) <- c("SNP","MAF","WT","P_Assoc","P_HW")

## Regional Test Results/Stats
REG <- array( ,c(N_FILES,8) )
colnames(REG) <- c("First","Last","SNPs","Rare_SNPs","BURD","BURD.wt","BURD.rare","BURD.wt.rare")
rownames(REG) <- gsub( ".txt","", FILE_LIST, fixed=T )

## Burden Test Stats
BURD <- BURD.wt <- BURD.rare <- BURD.wt.rare <- array( ,c(N_SAMPS, N_FILES) )
colnames(BURD) <- colnames(BURD.wt) <- colnames(BURD.rare) <- colnames(BURD.wt.rare) <- gsub( ".txt","", FILE_LIST, fixed=T )
rownames(BURD) <- rownames(BURD.wt) <- rownames(BURD.rare) <- rownames(BURD.wt.rare) <- SAMPS

###############################################################
## LOOP THROUGH BINS ##########################################
###############################################################

## Loop through FILE_LIST and load GT data
start_time <- proc.time()
for ( f in 1:N_FILES ) {
	# Specify File Name & Index Name
	file <- FILE_LIST[f]
	file_name <- paste( PathToGT, gsub( "txt", "raw", file ), sep="")
	name <- gsub( ".txt", "", file, fixed=T )
	print(paste( "Running:",f,"of",N_FILES,"-",name,"-",round(proc.time()-start_time,3)[3] ))
	# Load Genotype File
	GT <- read.table( file_name, header=T )
	GT.samps <- as.character( GT[,"FID"] )
	GT <- GT[ , 7:ncol(GT) ]
	rownames(GT) <- GT.samps
	GT <- GT[ SAMPS, ]
	# Load Hardy-Weinberg File
	HW <- read.table( gsub("raw","hwe",file_name), header=T )
	N_SNPS <- nrow(HW)
	## Compile Stats for this set
	GWAS.temp <- array( ,c( N_SNPS, ncol(GWAS) ) )
	colnames(GWAS.temp) <- colnames(GWAS)
	 # SNP Data
	GWAS.temp[,"SNP"] <- as.character( HW[,"SNP"] )
	GWAS.temp[,"P_HW"] <- as.numeric( HW[,"P"] )
	 # MAF
	MAF <- 2 * colMeans(GT) # apply( GT, 2, mean )
	MAF[ which(MAF>.5) ] <- 1 - MAF[ which(MAF>.5) ]
	GWAS.temp[,"MAF"] <- MAF
	 # Calculate Weight based on MAF
	b1 <- 1 ; b2 <- 10
	WTS <- dbeta( MAF, b1, b2 )
	GWAS[,"WT"] <- round( WTS, 5 )
	## Regional Stats
	REG[f,"First"] <- as.character( HW$SNP[1] )
	REG[f,"Last"] <- as.character( HW$SNP[N_SNPS] )
	HW.threshold <- 1e-30
	RM.HWE <- which( HW$P < HW.threshold )
	RM.MAF <- which( MAF < .01 )
	REG[f,"SNPs"] <- N_SNPS - length(RM.HWE)
	REG[f,"Rare_SNPs"] <- length(which( MAF<.01 )) - length(which( MAF[RM.HWE]<.01 ))

	#####################################################
	## Run Single-Locus Test on Relevant SNPs ###########
	 # Get rid of some SNPs
	print("Compiling Single-Locus Stats")
	RM <- union( RM.MAF, RM.HWE )
	KP <- setdiff( 1:N_SNPS, RM )
	 # Loop through Remaining SNPs
	for ( k in KP ) {
		# CP.temp <- data.frame( )
		MOD <- lm( RESIDS[SAMPS] ~ GT[ SAMPS, k ] )
		GWAS.temp[k,"P_Assoc"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	}
	GWAS <- rbind( GWAS, GWAS.temp )

	#####################################################
	## Run Burden Test on Region ###########
	## (Removing HWE Violations) ###########
	print("Compiling Burden Stats")
	 # Normal Burden Test
	if ( length(RM.HWE)>0 ) {
		BURD.samp <- rowSums(GT[ , -RM.HWE ])
	}else{ BURD.samp <- rowSums(GT) }
	BURD[,f] <- BURD.samp
	MOD <- lm( RESIDS[SAMPS] ~ BURD.samp[SAMPS] )
	REG[f,"BURD"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	 # Weighted Burden Test
	if ( length(RM.HWE)>0 ) {
		BURD.wt.samp <- rowSums( WTS[-RM.HWE] * GT[ , -RM.HWE ] )
		}else{ BURD.wt.samp <- rowSums( WTS * GT ) }
	BURD.wt[,f] <- BURD.wt.samp
	MOD <- lm( RESIDS[SAMPS] ~ BURD.wt.samp[SAMPS] )
	REG[f,"BURD.wt"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	# Rare Variants Burden Test
	RARE <- which( MAF < .01 )
	WHICH <- setdiff( RARE, RM.HWE )
	BURD.rare.samp <- rowSums(GT[ , WHICH ])
	BURD.rare[,f] <- BURD.rare.samp
	MOD <- lm( RESIDS[SAMPS] ~ BURD.rare.samp[SAMPS] )
	REG[f,"BURD.rare"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	 # Weighted Rare Variants Burden Test
	BURD.wt.rare.samp <- rowSums( WTS[WHICH] * GT[ , WHICH ] )
	BURD.wt.rare[,f] <- BURD.wt.rare.samp
	MOD <- lm( RESIDS[SAMPS] ~ BURD.wt.rare.samp[SAMPS] )
	REG[f,"BURD.wt.rare"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]

}

###############################################################
## RE-FORMAT OUTPUTS ##########################################
###############################################################

## GWAS Table
# GWAS <- data.frame( )

###############################################################
## PLOT RESULTS ###############################################
###############################################################

## Plot GWAS Results
XLIM <- c( 1, nrow(GWAS) )
YLIM <- c( 0, -log10( min(as.numeric(GWAS[,"P_Assoc"]),na.rm=T) ) )
jpeg( paste(PathToOut,"Pl_Gamut_GWAS.jpeg",sep=""), width=2000,height=1000,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Single-Locus Association",xlab="Position",ylab="-log10(P)" )
abline( h=seq(0,YLIM[2],1), lty=2,col="grey50",lwd=1 )
abline( h=-log10(.05/1e6), lty=2,col="firebrick2",lwd=1 )
points( 1:nrow(GWAS),-log10(as.numeric(GWAS[,"P_Assoc"])), col="slateblue3", pch="+" )
dev.off()

## Plot Burden Test Results
XLIM <- c( 1, nrow(REG) )
YLIM <- c( 0, -log10( min(as.numeric(REG[,5:8]),na.rm=T) ) )
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2")
jpeg( paste(PathToOut,"Pl_Gamut_REG_BURD.jpeg",sep=""), width=2000,height=1000,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Burden Test",xlab="Window",ylab="-log10(P)" )
abline( h=seq(0,YLIM[2],1), lty=2,col="grey50",lwd=1 )
points( 1:nrow(REG),-log10(as.numeric(REG[,"BURD"])), col=COLS[1], pch="+" )
points( 1:nrow(REG),-log10(as.numeric(REG[,"BURD.wt"])), col=COLS[2], pch="+" )
points( 1:nrow(REG),-log10(as.numeric(REG[,"BURD.rare"])), col=COLS[3], pch="+" )
points( 1:nrow(REG),-log10(as.numeric(REG[,"BURD.wt.rare"])), col=COLS[4], pch="+" )
legend( "topleft", pch="+",col=COLS, legend=colnames(REG)[5:8] )
dev.off()

## Plot Variant Counts
XLIM <- c( 1, nrow(REG) )
YLIM <- c( 0, max(as.numeric(REG[,3:4]),na.rm=T) )
COLS <- c("tomato2","cadetblue2","darkorchid2")
jpeg( paste(PathToOut,"Pl_Gamut_REG_VARS.jpeg",sep=""), width=2000,height=1000,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Variant Counts",xlab="Window",ylab="# Variants" )
abline( h=seq(0,YLIM[2],100), lty=2,col="grey50",lwd=1 )
points( 1:nrow(REG),as.numeric(REG[,"SNPs"]), col=COLS[1], pch="+" )
points( 1:nrow(REG),as.numeric(REG[,"Rare_SNPs"]), col=COLS[2], pch="+" )
points( 1:nrow(REG),as.numeric(REG[,"SNPs"])-as.numeric(REG[,"Rare_SNPs"]), col=COLS[3], pch="+" )
legend( "bottomleft", pch="+",col=COLS, legend=c(colnames(REG)[3:4],"Common_SNPs") )
dev.off()

###############################################################
## WRITE TABLES ###############################################
###############################################################

write.table( GWAS, paste(PathToOut,"Gamut_GWAS.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

write.table( REG, paste(PathToOut,"Gamut_REG.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

write.table( BURD, paste(PathToOut,"Gamut_BURD.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
write.table( BURD.wt, paste(PathToOut,"Gamut_BURD.wt.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
write.table( BURD.rare, paste(PathToOut,"Gamut_BURD.rare.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
write.table( BURD.wt.rare, paste(PathToOut,"Gamut_BURD.wt.rare.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )










###############################################################
## END OF DOC #################################################
###############################################################














