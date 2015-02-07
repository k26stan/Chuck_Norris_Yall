## Python Script to Split Variant IDs into Separate Lists by Region ##
## Split by Number of BASEs ##
## To be used in Walker.sh Script ##
## Feburary 6, 2015 ##
## Kristopher Standish ##

## Usage ##
# Rscript 3-Split_Var_List_SNP.py <Path/To/Variant_List.bim> <Path/To/Out/Dir> <Number_of_SNPs>
# Rscript 3-Split_Var_List_SNP.py /projects/janssen/VCFs/PLINK/BED_FULL.ALL.bim /projects/janssen/Walker/20150206_NAME/ 1000
# Rscript 3-Split_Var_List_SNP.py ${VAR_FILE%%.bed}.bim ${OUT_DIR} ${NUM_UNITS}

from sys import argv # standard packages to allow imput params
from sys import exit
import gzip
import datetime

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

# Checks for the correct number of input params, displays usage info if needed
try:
    script, Var_List, PathToOut, Num_SNPs_Str = argv
except ValueError:
    print "\nScript used to walk through variants and create lists for each region."
    exit("cmd_usage$> script, Var_List, PathToOut, Num_SNPs_Str \n")

print "Variant List: %s" % Var_List
print "Write Path: %s" % PathToOut
print "Number of SNPs: %s" % Num_SNPs_Str
Num_SNPs = int( Num_SNPs_Str )

###############################################################
## LOAD DATA ##################################################
###############################################################

## Open Position file
var = open(Var_List)
print "Var_List open - %s" % datetime.datetime.now().time().isoformat()

###############################################################
## SPLIT VARIANT LIST #########################################
###############################################################

## Specify Naming Parameters for first File List
which_chrom = 1

## Specify/Split First Line in Variant List File
var_line = var.next()
split_line = var_line.strip().split()
line_chr = split_line[0]

## Skip all variants listed prior to Chromsome 1 (e.g., chr=0)
while line_chr == 0:
	var_line = var.next()
	split_line = var_line.strip().split()
	line_chr = split_line[0]

## Specify First File Name
line_coord = int( split_line[3] )
file_name = ( "%s/Chr%s_%s.txt" ) % ( PathToOut, line_chr, line_coord )

## Open File to Write to
w = open( file_name, 'w' )

## Write first SNP to file
w.write( split_line[1] + "\n" )

## Start running through interesting part of genome
 # Write SNP ID's under certain circumstances:
   # 1) Number of SNPs written is still less than Num_SNPs
   # 2) SNP is on same chromosome as last SNP
   # 3) Next SNP is within 1 Mb of last SNP

# Specify Write Count
count = 1
for var_line in var:
	# Specify Previous Chromosome and Coordinates
	line_chr_prev = line_chr
	line_coord_prev = line_coord
	# Split Current Line
	split_line = var_line.strip().split()
	# Get Chromosome and Coordinate
	line_chr = split_line[0]
	line_coord = int( split_line[3] )
	# Check Write Criteria
	if ( line_chr == line_chr_prev ) and ( line_coord-line_coord_prev < 1000000 ) and ( count < Num_SNPs ):
		# Write to File
		w.write( split_line[1] + "\n" )
		count += 1
	else:
		# Close old file, open new
		w.close()
		file_name = ( "%s/Chr%s_%s.txt" ) % ( PathToOut, line_chr, line_coord )
		w = open( file_name, 'w' )
		# Write to File
		w.write( split_line[1] + "\n" )
		count = 1

## Close All Open File Threads
var.close()
w.close()



###############################################################
## END OF DOC #################################################
###############################################################



















