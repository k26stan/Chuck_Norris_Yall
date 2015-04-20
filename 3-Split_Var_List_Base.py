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
import math

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

# Checks for the correct number of input params, displays usage info if needed
try:
    script, Var_List, PathToOut, Num_Bases_Str = argv
except ValueError:
    print "\nScript used to walk through variants and create lists for each region."
    exit("cmd_usage$> script, Var_List, PathToOut, Num_Bases_Str \n")

print "Variant List: %s" % Var_List
print "Write Path: %s" % PathToOut
print "Number of Bases: %s" % Num_Bases_Str
Num_Bases = float( Num_Bases_Str )

###############################################################
## LOAD DATA ##################################################
###############################################################

## Open Position file
var = open(Var_List)
print "Var_List open - %s" % datetime.datetime.now().time().isoformat()

## Create Rounding Function for Cordinates
def rounddown(x):
	return int(math.floor(x / Num_Bases)) * Num_Bases

def roundup(x):
	return int(math.ceil(x / Num_Bases)) * Num_Bases

###############################################################
## SPLIT VARIANT LIST #########################################
###############################################################

## Specify/Split First Line in Variant List File
var_line = var.next()
split_line = var_line.strip().split()
line_chr = split_line[0]

## Skip all variants listed prior to Chromsome 1 (e.g., chr=0)
while line_chr == 0:
	var_line = var.next()
	split_line = var_line.strip().split()
	line_chr = split_line[0]

## Specify Coordinate of First Variant
line_coord = float( split_line[3] )

## Specify Bin Start/End Points
bin_start = int( rounddown( line_coord ) )
bin_end = int( roundup( line_coord ) - 1 )

## Specify File Name for Output
file_name = ( "%s/Chr%s_%s-%s.txt" ) % ( PathToOut, line_chr, bin_start, bin_end )

## Open File to Write to
w = open( file_name, 'w' )

## Write first SNP to file
w.write( split_line[1] + "\n" )

## Start running through interesting part of genome
 # Write SNP ID's under certain circumstances:
   # SNPs are on Same Chromosome
   # SNPs are in same bin of size Num_Bases

# Specify Write Count
for var_line in var:
	# Specify Previous Chromosome and Coordinates
	line_chr_prev = line_chr
	line_coord_prev = line_coord
	# Split Current Line
	split_line = var_line.strip().split()
	# Get Chromosome and Coordinate
	line_chr = split_line[0]
	line_coord = float( split_line[3] )
	## Check Write Criteria
	 # Same Chromosome?
	if line_chr != line_chr_prev:
		# Close old file, open new
		w.close()
		# Specify Bin and File Name for New Output
		bin_start = int( rounddown( line_coord ) )
		bin_end = int( roundup( line_coord ) - 1 )
		file_name = ( "%s/Chr%s_%s-%s.txt" ) % ( PathToOut, line_chr, bin_start, bin_end )
		# Open File to Write to
		w = open( file_name, 'w' )
	if line_coord < bin_end:
		# Write to File
		w.write( split_line[1] + "\n" )
	else:
		# Close old file, open new
		w.close()
		# Specify Bin and File Name for New Output
		bin_start = int( rounddown( line_coord ) )
		bin_end = int( roundup( line_coord ) - 1 )
		file_name = ( "%s/Chr%s_%s-%s.txt" ) % ( PathToOut, line_chr, bin_start, bin_end )
		# Open File to Write to
		w = open( file_name, 'w' )
		# Write to File
		w.write( split_line[1] + "\n" )


## Close All Open File Threads
var.close()
w.close()



###############################################################
## END OF DOC #################################################
###############################################################



















