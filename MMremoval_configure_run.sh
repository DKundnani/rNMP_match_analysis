#! usr/bin/bash
#Author: Deepali L. Kundnani
#Institute: Georgia institute of Technology(GT)
#Description: Configure run file for batch of RiboseMap generated bed files to filter rNMPs that match at the rNMP location while removing tailing nucleotide mismatches upstrema of rNMP 


################################## Define Variables ---------TO BE CONFIGURED BY THE USER----------------
echo ".....ASSIGNING VARIABLES...."

pymatch='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/rNMP_match_analysis/match_analysis.py' #match analysis script
reference='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa' #Fasta file for reference genome
genome=$(echo $reference.fai) #This assumes you have respective genome.sizes in the same folder with etx as .fai to the fasta file. You can edit if needed
polyN='5' # number of polyN's to check including the ribo position Default=5: Will cleanup upto AAAArA mismatches in dA tialing library
bedloc='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/Yeast_all/bed/Shortnuclrepeat' #location with bed files - Make sure you have the bed files only required for analysis, or it will analyze all the files
fastqloc='/storage/coda1/p-fstorici3/0/shared/raw_reads' #Make sure the basename of files to be analysed is same as the basename of '.fq' files. Sample.bed used Sample.fq for analysis
resultsloc='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/Yeast_all/bed/Shortnuclrepeat/all_polyNfiltered/' #location you want the results in 
RrmpolyN='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/rNMP_match_analysis/removepolyN.R' #removepolyN.R script 
Ntail='A' #Nucleotide used for dN tailing (Default-'A')


################################## Setup environment ---------TO BE CONFIGURED BY THE USER----------------
echo ".....ACTIVATING ENVIROMENT...."
source activate r_env
source /storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/rNMP_match_analysis/mismatch_cleanup.sh #run to activate following functions

################################### DO NOT CHANGE/CONFIGURE BELOW THIS LINE ###################################
################################### Running the functions ---------TO BE RUN AS IT IS & STEP WISE----------------
echo ".....RUNNING FUNCTIONS...."
#Run the funcitons one by one and wait for them to finish OR just run this script

mismatch_analysis 		# Will run mismatch_analysis and generate files with prefix 'MManalysis' giving the reference and fastq sequences for 'polyN' bases upstream of rNMP 
filter_matches			# Filter matches from the previous results
remove_polyN 			# Will filter polyN
cal_mm_percent chrM			# Will calculate mismatch percentages in the input bed files and filtered files and separate chrM from other chromosomes

##################################
