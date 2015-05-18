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
# LINE <- c("/projects/janssen/Walker/20150408_LT8_DEL_MNe_MN_DAS_BL_MN_SEX_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150408_LT8_DEL_MNe_MN_DAS_BL_MN_SEX_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,SEX,PC1,PC2")
# LINE <- c("/projects/janssen/Walker/20150409_LT8_DEL_MNe_MN_DAS_BL_MN_SEX_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150409_LT8_DEL_MNe_MN_DAS_BL_MN_SEX_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,SEX,PC1,PC2")
# LINE <- c("/projects/janssen/Walker/20150424_LT8_DEL_MNe_MN_DAS_BL_MN_SEX_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150424_LT8_DEL_MNe_MN_DAS_BL_MN_SEX_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,SEX,PC1,PC2")
# LINE <- c("/projects/janssen/Walker/20150430_WG_1000SNP/20150430_3_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150430_WG_1000SNP/20150430_3_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,PC1,PC2")
# LINE <- c("/projects/janssen/Walker/20150508_WG_100KB/20150508_21_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/Genotype_Files/SNP_Lists.txt","/projects/janssen/Walker/20150508_WG_100KB/20150508_21_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/Cov_w_PCs.txt","BURD,SKAT,ELNET","DAS_BL_MN,PC1,PC2")
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
GWAS_COLNAMES <- c("SNP","MAF","WT","P_Assoc","P_HW")
GWAS <- array( ,c(0,length(GWAS_COLNAMES)) )
colnames(GWAS) <- GWAS_COLNAMES

## Regional Test Results/Stats
REG_COLNAMES <- c("First","Last","Mid","SNPs","Rare_SNPs","BURD","BURD.wt","BURD.rare","BURD.wt.rare","F.temp","MULT.p","MULT.rare.f","MULT.rare.p")
REG <- array( ,c(N_FILES,length(REG_COLNAMES)) )
colnames(REG) <- REG_COLNAMES
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
# for ( f in 1:500 ) {
# for ( f in 1330:N_FILES ) {
	# Specify File Name & Index Name
	file <- FILE_LIST[f]
	file_name <- paste( PathToGT, gsub( "txt", "raw", file ), sep="")
	name <- gsub( ".txt", "", file, fixed=T )
	position <- strsplit( name, "_" )[[1]][2]
	print(paste( "Running:",f,"of",N_FILES,"-",name,"-",round(proc.time()-start_time,3)[3] ))
	# Load Genotype File
	GT <- read.table( file_name, header=T )
	GT.samps <- as.character( GT[,"FID"] )
	if ( ncol(GT)==7 ) { 
		REG[f,"First"] <- colnames(GT)[7]
		REG[f,"Last"] <- colnames(GT)[7]
		REG[f,"Mid"] <- position
		REG[f,"SNPs"] <- 1
		next
	}
	GT <- GT[ , 7:ncol(GT) ]
	rownames(GT) <- GT.samps
	GT <- GT[ SAMPS, ]
	# GT.maf.0 <- which( colSums(GT[SAMPS,])==0 )
	# GT <- which
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
	MAF <- colMeans(GT, na.rm=T) / 2 # apply( GT, 2, mean )
	SD <- apply( GT,2,sd, na.rm=T)
	MAF[ which(MAF>.5) ] <- 1 - MAF[ which(MAF>.5) ]
	GWAS.temp[,"MAF"] <- MAF
	 # Calculate Weight based on MAF
	b1 <- 1 ; b2 <- 10
	WTS <- dbeta( MAF, b1, b2 )
	GWAS.temp[,"WT"] <- round( WTS, 5 )
	## Regional Stats
	REG[f,"First"] <- as.character( HW$SNP[1] )
	REG[f,"Last"] <- as.character( HW$SNP[N_SNPS] )
	REG[f,"Mid"] <- position
	HW.threshold <- 1e-20
	RM.HWE <- which( HW$P < HW.threshold | is.na(MAF) | SD==0 | is.na(SD) )
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
		BURD.samp <- rowSums(GT[ , -RM.HWE ], na.rm=T)
	}else{ BURD.samp <- rowSums(GT, na.rm=T) }
	BURD[,f] <- BURD.samp
	MOD <- lm( RESIDS[SAMPS] ~ BURD.samp[SAMPS] )
	REG[f,"BURD"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	 # Weighted Burden Test
	if ( length(RM.HWE)>0 ) {
		BURD.wt.samp <- rowSums( WTS[-RM.HWE] * GT[ , -RM.HWE ], na.rm=T )
		}else{ BURD.wt.samp <- rowSums( WTS * GT, na.rm=T ) }
	BURD.wt[,f] <- BURD.wt.samp
	MOD <- lm( RESIDS[SAMPS] ~ BURD.wt.samp[SAMPS] )
	REG[f,"BURD.wt"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	## Rare Variants Burden Test
	RARE <- which( MAF < .01 )
	WHICH <- setdiff( RARE, RM.HWE )
	if ( length(WHICH)>0 ) {
		BURD.rare.samp <- rowSums( data.frame(GT[ , WHICH ]) )
		BURD.rare[,f] <- BURD.rare.samp
		MOD <- lm( RESIDS[SAMPS] ~ BURD.rare.samp[SAMPS] )
		REG[f,"BURD.rare"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
		 # Weighted Rare Variants Burden Test
		BURD.wt.rare.samp <- rowSums( WTS[WHICH] * data.frame(GT[ , WHICH ]) )
		BURD.wt.rare[,f] <- BURD.wt.rare.samp
		MOD <- lm( RESIDS[SAMPS] ~ BURD.wt.rare.samp[SAMPS] )
		REG[f,"BURD.wt.rare"] <- summary(MOD)$coefficients[ 2,"Pr(>|t|)" ]
	}else{
		print( "IN HERE!!!" )
		BURD.rare.samp <- rowSums(GT[ , WHICH ])
		BURD.rare[,f] <- BURD.rare.samp
		REG[f,"BURD.rare"] <- NA
		 # Weighted Rare Variants Burden Test
		BURD.wt.rare.samp <- rowSums(GT[ , WHICH ])
		BURD.wt.rare[,f] <- BURD.wt.rare.samp
		REG[f,"BURD.wt.rare"] <- NA
	}

	#####################################################
	## Run Mult-Regression Test on Region ###########
	print("Running Multiple Regression")
	if ( length(RM.HWE)>0 ) {
		MULT <- merge( data.frame(RESIDS),GT[ , -RM.HWE ], by="row.names" )
	}else{ MULT <- merge( data.frame(RESIDS),GT, by="row.names" ) }
	MULT.2 <- lm( RESIDS ~ ., data=MULT[,2:ncol(MULT)] )
	MULT.3 <- summary(MULT.2)
	F.temp <- MULT.3$fstatistic
	MULT.p <- pf(F.temp[1],F.temp[2],F.temp[3],lower.tail=F)
	REG[f,"MULT.p"] <- MULT.p
	REG[f,"F.temp"] <- F.temp[1]

	if ( length(WHICH)>0 ) {
		MULT <- merge( data.frame(RESIDS),data.frame(GT[ , WHICH ],row.names=rownames(GT)), by="row.names" )
		MULT.2 <- lm( RESIDS ~ ., data=MULT[,2:ncol(MULT)] )
		MULT.3 <- summary(MULT.2)
		F.temp <- MULT.3$fstatistic
		MULT.p <- pf(F.temp[1],F.temp[2],F.temp[3],lower.tail=F)
		REG[f,"MULT.rare.p"] <- MULT.p
		REG[f,"MULT.rare.f"] <- F.temp[1]
	}else{
		REG[f,"MULT.rare.p"] <- NA
		REG[f,"MULT.rare.f"] <- NA
	}
	# MULT <- merge( data.frame(RESIDS),GT[,WHICH], by="row.names" )


	#####################################################
	## Save Tables every few Loops ###########	
	if ( f %% round(nrow(REG)/10,0) == 0 ) {
		print( "WRITING TABLES!!" )
		## Write GWAS Results Table
		write.table( GWAS, paste(PathToOut,"Gamut_GWAS.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

		## Write Regression Results Table
		write.table( REG, paste(PathToOut,"Gamut_REG.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

		## Write Burden Stats Table
		write.table( BURD, paste(PathToOut,"Gamut_BURD.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
		write.table( BURD.wt, paste(PathToOut,"Gamut_BURD.wt.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
		write.table( BURD.rare, paste(PathToOut,"Gamut_BURD.rare.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
		write.table( BURD.wt.rare, paste(PathToOut,"Gamut_BURD.wt.rare.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
	}

} # Close File Loop

###############################################################
## WRITE TABLES ###############################################
###############################################################

## Write GWAS Results Table
write.table( GWAS, paste(PathToOut,"Gamut_GWAS.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

## Write Regression Results Table
write.table( REG, paste(PathToOut,"Gamut_REG.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

## Write Burden Stats Table
write.table( BURD, paste(PathToOut,"Gamut_BURD.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( BURD.wt, paste(PathToOut,"Gamut_BURD.wt.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( BURD.rare, paste(PathToOut,"Gamut_BURD.rare.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( BURD.wt.rare, paste(PathToOut,"Gamut_BURD.wt.rare.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )

###############################################################
## RE-LOAD FILES FOR MANUAL PLOTTING ##########################
###############################################################

# for file in `ls Gamut_*`; do mv ${file} ${file%%.txt}.1-1340.txt; done

# for file in `ls Gamut_*`; do mv ${file} ${file%%.txt}.651-1000.txt; done
# for file in `ls Gamut_*`; do mv ${file} ${file%%.txt}.1001-End.txt; done
# for file in `ls Gamut_*1-650.651-1000.1001-End.txt`; do mv ${file} ${file%%1-650.651-1000.1001-End.txt}.txt; done
# for file in `ls Gamut_*..txt`; do mv ${file} ${file%%..txt}.1-650.txt; done
# for file in `ls Gamut_*651-1000.1001-End.txt`; do mv ${file} ${file%%651-1000.1001-End.txt}651-1000.txt; done

# ######### 1-650 ###############
# ## GWAS Results Tables
# GWAS.1 <- read.table( paste(PathToOut,"Gamut_GWAS.1-650.txt",sep=""), sep="\t",header=T )

# ## Regression Results Table
# REG.1 <- read.table( paste(PathToOut,"Gamut_REG.1-650.txt",sep=""), sep="\t",header=T )

# ## Burden Stats Table
# BURD.1 <- read.table( paste(PathToOut,"Gamut_BURD.1-650.txt",sep=""), sep="\t",header=T )
# BURD.wt.1 <- read.table( paste(PathToOut,"Gamut_BURD.wt.1-650.txt",sep=""), sep="\t",header=T )
# BURD.rare.1 <- read.table( paste(PathToOut,"Gamut_BURD.rare.1-650.txt",sep=""), sep="\t",header=T )
# BURD.wt.rare.1 <- read.table( paste(PathToOut,"Gamut_BURD.wt.rare.1-650.txt",sep=""), sep="\t",header=T )

# ######### 651-1000 ###############
# ## GWAS Results Tables
# GWAS.2 <- read.table( paste(PathToOut,"Gamut_GWAS.651-1000.txt",sep=""), sep="\t",header=T )

# ## Regression Results Table
# REG.2 <- read.table( paste(PathToOut,"Gamut_REG.651-1000.txt",sep=""), sep="\t",header=T )

# ## Burden Stats Table
# BURD.2 <- read.table( paste(PathToOut,"Gamut_BURD.651-1000.txt",sep=""), sep="\t",header=T )
# BURD.wt.2 <- read.table( paste(PathToOut,"Gamut_BURD.wt.651-1000.txt",sep=""), sep="\t",header=T )
# BURD.rare.2 <- read.table( paste(PathToOut,"Gamut_BURD.rare.651-1000.txt",sep=""), sep="\t",header=T )
# BURD.wt.rare.2 <- read.table( paste(PathToOut,"Gamut_BURD.wt.rare.651-1000.txt",sep=""), sep="\t",header=T )

# ######### 1001-End ###############
# ## GWAS Results Tables
# GWAS.3 <- read.table( paste(PathToOut,"Gamut_GWAS.1001-End.txt",sep=""), sep="\t",header=T )

# ## Regression Results Table
# REG.3 <- read.table( paste(PathToOut,"Gamut_REG.1001-End.txt",sep=""), sep="\t",header=T )

# ## Burden Stats Table
# BURD.3 <- read.table( paste(PathToOut,"Gamut_BURD.1001-End.txt",sep=""), sep="\t",header=T )
# BURD.wt.3 <- read.table( paste(PathToOut,"Gamut_BURD.wt.1001-End.txt",sep=""), sep="\t",header=T )
# BURD.rare.3 <- read.table( paste(PathToOut,"Gamut_BURD.rare.1001-End.txt",sep=""), sep="\t",header=T )
# BURD.wt.rare.3 <- read.table( paste(PathToOut,"Gamut_BURD.wt.rare.1001-End.txt",sep=""), sep="\t",header=T )

# ######### 1001-End ###############
# ## GWAS Results Tables
# GWAS <- rbind( GWAS.1, GWAS.2, GWAS.3 )

# ## Regression Results Table
# REG <- rbind( REG.1[1:650,], REG.2[651:1000,], REG.3[1001:nrow(REG.3),] )

#####################################################

# ######### 1001-End ###############
# ## GWAS Results Tables
# # GWAS.1 <- read.table( paste(PathToOut,"Gamut_GWAS.1-500.txt",sep=""), sep="\t",header=T )
# GWAS.2 <- read.table( paste(PathToOut,"Gamut_GWAS.txt",sep=""), sep="\t",header=T )
# GWAS.3 <- read.table( paste(PathToOut,"Gamut_GWAS.500-End.txt",sep=""), sep="\t",header=T )

# # ## Regression Results Table
# # REG.1 <- read.table( paste(PathToOut,"Gamut_REG.1-500.txt",sep=""), sep="\t",header=T )
# REG.2 <- read.table( paste(PathToOut,"Gamut_REG.txt",sep=""), sep="\t",header=T )
# REG.3 <- read.table( paste(PathToOut,"Gamut_REG.500-End.txt",sep=""), sep="\t",header=T )

# ######### 1001-End ###############
# ## GWAS Results Tables
# WHICH <- grep( GWAS.2$SNP[nrow(GWAS.2)], GWAS.3$SNP )
# GWAS <- rbind( GWAS.2, GWAS.3[WHICH:nrow(GWAS.3),] )

# ## Regression Results Table
# WHICH <- tail(which(!is.na(REG.2[,1])),1)
# # tail(which(is.na(REG.3[,1])))
# REG <- rbind( REG.2[1:WHICH,], REG.3[(WHICH+1):nrow(REG.3),] )

#####################################################

######### 1001-End ###############
## GWAS Results Tables
# GWAS.1 <- read.table( paste(PathToOut,"Gamut_GWAS.1-500.txt",sep=""), sep="\t",header=T )
# GWAS.2 <- read.table( paste(PathToOut,"Gamut_GWAS.1-1340.txt",sep=""), sep="\t",header=T )
# GWAS.3 <- read.table( paste(PathToOut,"Gamut_GWAS.txt",sep=""), sep="\t",header=T )

# # ## Regression Results Table
# # REG.1 <- read.table( paste(PathToOut,"Gamut_REG.1-500.txt",sep=""), sep="\t",header=T )
# REG.2 <- read.table( paste(PathToOut,"Gamut_REG.1-1340.txt",sep=""), sep="\t",header=T )
# REG.3 <- read.table( paste(PathToOut,"Gamut_REG.txt",sep=""), sep="\t",header=T )

# ######### 1001-End ###############
# ## GWAS Results Tables
# WHICH <- grep( GWAS.2$SNP[nrow(GWAS.2)], GWAS.3$SNP )
# GWAS <- rbind( GWAS.2, GWAS.3[WHICH:nrow(GWAS.3),] )

# ## Regression Results Table
# WHICH <- tail(which(!is.na(REG.2[,1])),1)
# # tail(which(is.na(REG.3[,1])))
# REG <- rbind( REG.2[1:WHICH,], REG.3[(WHICH+1):nrow(REG.3),] )


PathToOut <- "20150508_WG_100KB/20150508_22_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/"
GWAS <- read.table( "20150508_WG_100KB/20150508_22_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/Gamut_GWAS.txt",sep="\t",header=T)
REG <- read.table( "20150508_WG_100KB/20150508_22_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/Gamut_REG.txt",sep="\t",header=T)

###############################################################
## PLOT RESULTS ###############################################
###############################################################

## Plot GWAS Results
GWAS.plot <- GWAS[ which(!is.na(GWAS$P_Assoc)), ]
XLIM <- c( 1, nrow(GWAS.plot) )
YLIM <- c( 0, -log10( min(as.numeric(GWAS.plot[,"P_Assoc"]),na.rm=T) ) )
jpeg( paste(PathToOut,"Pl_Gamut_GWAS.jpeg",sep=""), width=2000,height=1000,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Single-Locus Association",xlab="Position",ylab="-log10(P)" )
abline( h=seq(0,YLIM[2]+3,1), lty=2,col="grey50",lwd=1 )
abline( h=-log10(.05/1e6), lty=2,col="firebrick2",lwd=1 )
points( 1:nrow(GWAS.plot),-log10(as.numeric(GWAS.plot[,"P_Assoc"])), col="slateblue3", pch="+" )
dev.off()

## Plot Burden Test Results
XLIM <- range( as.numeric(REG[,"Mid"]), na.rm=T ) # c( 1, nrow(REG) )
YLIM <- c( 0,7 )
YLIM <- c( 0, -log10( min(as.numeric( REG[,grep("BURD",colnames(REG))] ),na.rm=T) ) )
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2","chocolate2","purple2")
jpeg( paste(PathToOut,"Pl_Gamut_REG_BURD.jpeg",sep=""), width=2000,height=1000,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Burden Test",xlab="Window",ylab="-log10(P)" )
abline( h=seq(0,YLIM[2],1), lty=2,col="grey50",lwd=1 )
points( REG[,"Mid"],-log10(as.numeric(REG[,"BURD"])), col=COLS[1], pch="+" )
points( REG[,"Mid"],-log10(as.numeric(REG[,"BURD.wt"])), col=COLS[2], pch="+" )
points( REG[,"Mid"],-log10(as.numeric(REG[,"BURD.rare"])), col=COLS[3], pch="+" )
points( REG[,"Mid"],-log10(as.numeric(REG[,"BURD.wt.rare"])), col=COLS[4], pch="+" )
points( REG[,"Mid"],-log10(as.numeric(REG[,"MULT.p"])), col=COLS[5], pch="+" )
points( REG[,"Mid"],-log10(as.numeric(REG[,"MULT.rare.p"])), col=COLS[6], pch="+" )
legend( "bottomleft", pch="+",col=COLS, legend=c(colnames(REG)[grep("BURD",colnames(REG))],"MULT.p","MULT.rare.p"),ncol=3 )
dev.off()

## Plot Variant Counts
XLIM <- range( as.numeric(REG[,"Mid"]), na.rm=T ) # c( 1, nrow(REG) )
# YLIM <- c( 0,1.2*REG[1,"SNPs"] )
# YLIM <- c( 0,max(REG[,"SNPs"],na.rm=T) )
YLIM <- c( 0, max(as.numeric(REG[,grep("SNPs",colnames(REG))]),na.rm=T) )
COLS <- c("tomato2","cadetblue2","darkorchid2")
jpeg( paste(PathToOut,"Pl_Gamut_REG_VARS.jpeg",sep=""), width=2000,height=1000,pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM, main="Variant Counts",xlab="Window",ylab="# Variants" )
abline( h=seq(0,YLIM[2],100), lty=2,col="grey50",lwd=1 )
points( REG[,"Mid"],as.numeric(REG[,"SNPs"]), col=COLS[1], pch="+" )
points( REG[,"Mid"],as.numeric(REG[,"Rare_SNPs"]), col=COLS[2], pch="+" )
points( REG[,"Mid"],as.numeric(REG[,"SNPs"])-as.numeric(REG[,"Rare_SNPs"]), col=COLS[3], pch="+" )
legend( "bottomleft", pch="+",col=COLS, legend=c(colnames(REG)[grep("SNPs",colnames(REG))],"Common_SNPs") )
dev.off()

## QQ Plot for Burden Tests
LIM <- c( 0,7 )
LIM <- c( 0, -log10( min(as.numeric( REG[,grep("BURD",colnames(REG))] ),na.rm=T) ) )
EXP <- (1:nrow(REG)+1) / (nrow(REG)+1)
N.REG.rare <- length(which(!is.na(REG[,"BURD.rare"])))
EXP.rare <- (1:N.REG.rare+1) / (N.REG.rare+1)
COLS <- c("firebrick2","gold2","chartreuse2","deepskyblue2","chocolate2","purple2")
jpeg( paste(PathToOut,"Pl_Gamut_REG_QQ.jpeg",sep=""), width=1400,height=1400,pointsize=30 )
plot( 0,0,type="n", xlim=LIM,ylim=LIM, xlab="-log10(EXP)",ylab="-log10(OBS)",main="QQ Plot for Bins" )
abline( h=seq(0,10,1),lty=2,col="grey50" ) ; abline( v=seq(0,10,1),lty=2,col="grey50" )
abline( 0,1, lty=1,col="black" )
abline( h=-log10(.05/nrow(REG)), lty=2,col="firebrick2",lwd=2)
COL_NAMES <- c("BURD","BURD.wt","BURD.rare","BURD.wt.rare","MULT.p","MULT.rare.p")
for ( i in 1:6 ) {
	WHICH_COL <- COL_NAMES[i]
	EXP <- (1:length(sort(REG[,WHICH_COL]))+1) / (length(sort(REG[,WHICH_COL]))+1)
	points( -log10(EXP), -log10(sort(as.numeric(REG[,WHICH_COL]))), col=COLS[i],pch="+",type="o" )
}
# EXP.MULT.p <- (1:length(sort(REG[,"BURD"]))+1) / (length(sort(REG[,"BURD"]))+1)
# points( -log10(EXP), -log10(sort(as.numeric(REG[,"BURD"]))), col=COLS[1],pch="+",type="o" )
# points( -log10(EXP), -log10(sort(as.numeric(REG[,"BURD.wt"]))), col=COLS[2],pch="+",type="o" )
# points( -log10(EXP.rare), -log10(sort(as.numeric(REG[,"BURD.rare"]))), col=COLS[3],pch="+",type="o" )
# points( -log10(EXP.rare), -log10(sort(as.numeric(REG[,"BURD.wt.rare"]))), col=COLS[4],pch="+",type="o" )
# EXP.MULT.p <- (1:length(sort(REG[,"MULT.p"]))+1) / (length(sort(REG[,"MULT.p"]))+1)
# points( -log10(EXP.MULT.p), -log10(sort(as.numeric(REG[,"MULT.p"]))), col=COLS[5],pch="+",type="o" )
# EXP.MULT.rare.p <- (1:length(sort(REG[,"MULT.rare.p"]))+1) / (length(sort(REG[,"MULT.rare.p"]))+1)
# points( -log10(EXP.MULT.rare.p), -log10(sort(as.numeric(REG[,"MULT.rare.p"]))), col=COLS[6],pch="+",type="o" )
legend( "topleft", pch="+",col=COLS, legend=c(colnames(REG)[grep("BURD",colnames(REG))],"MULT.p","MULT.rare.p") )
dev.off()

###############################################################
## MULTIPLE REGRESSION w/ BURDEN STATS ########################
###############################################################

# MG.BURD.wt.rare <- merge( data.frame(RESIDS), BURD.wt.rare, by="row.names" )
# SORT.BURD.wt.rare <- order( as.numeric(REG[,"BURD.wt.rare"]), decreasing=F )
# How_Many <- 10
# CAND.BURD.wt.rare <- SORT.BURD.wt.rare[1:How_Many] # SORT.BURD.wt.rare[2] # 
# MOD.BURD.wt.rare <- lm( RESIDS ~ ., data=MG.BURD.wt.rare[,c(2,CAND.BURD.wt.rare+2)] )
# summary(MOD.BURD.wt.rare)

# Ps.wt.rare <- R2.adj.wt.rare <- numeric(How_Many)
# for ( i in 1:How_Many ) {
# 	CAND.BURD.wt.rare <- SORT.BURD.wt.rare[1:i] # SORT.BURD.wt.rare[2] # 
# 	MOD.BURD.wt.rare <- lm( RESIDS ~ ., data=MG.BURD.wt.rare[,c(2,CAND.BURD.wt.rare+2)] )
# 	F.temp <- summary(MOD.BURD.wt.rare)$fstatistic
# 	Ps.wt.rare[i] <- pf(F.temp[1],F.temp[2],F.temp[3],lower.tail=F)
# 	R2.adj.wt.rare[i] <- summary(MOD.BURD.wt.rare)$adj.r.squared
# }



# MG.BURD.wt.rare <- merge( data.frame(RESIDS), BURD.wt.rare, by="row.names" )
# SORT.BURD.wt.rare <- order( as.numeric(REG[,"BURD.wt.rare"]), decreasing=F )
# How_Many <- 10
# CAND.BURD.wt.rare <- SORT.BURD.wt.rare[1:How_Many] # SORT.BURD.wt.rare[2] # 
# MOD.BURD.wt.rare <- lm( RESIDS ~ ., data=MG.BURD.wt.rare[,c(2,CAND.BURD.wt.rare+2)] )
# summary(MOD.BURD.wt.rare)
 

###############################################################
## END OF DOC #################################################
###############################################################














