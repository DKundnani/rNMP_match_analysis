#! usr/bin/bash
#Author: Deepali L. Kundnani
#Institute: Georgia institute of Technology(GT)
#Description: Filters Ribosemap generated bedfile from polyN found in ribose-seq sequencing data(upstream / 5' end from ribo) that doesn't match to the corresponding ribo position in reference genome (N=nucletoide used for dN tailing)

#pymatch='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/rNMP_match_analysis/match_analysis.py' #match analysis script
#reference='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/reference/sacCer2/sacCer2.fa' #Fasta file for reference genome
#genome=$(echo $reference.fai) #This assumes you have respective genome.sizes in the same folder with etx as .fai to the fasta file
#polyN='5' # # number of polyN's to check including hte ribo position
#bedloc='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/results/coordinate/test' #location with bed files - Make sure you have the bed files only required for analysis, or it will analyze all the files
#fastqloc='/storage/coda1/p-fstorici3/0/shared/raw_reads/' #Make sure the basename of files to be analysed is same as the basename of '.fq' files. Sample.bed used Sample.fq for analysis
#resultsloc='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/AGS/results/coordinate/testresults' #location you want the results in 
#RrmpolyN='/storage/home/hcoda1/5/dkundnani3/p-fstorici3-0/rich_project_bio-storici/bin/rNMP_match_analysis/removepolyN.R' #removepolyN.R script 
#Ntail='A' #Nucleotide used for dN tailing (Default-'A')

mismatch_analysis() {
	mkdir $resultsloc
    for bed in $(ls $bedloc/*.bed)
    do
    fastq=$(echo $fastqloc/$(basename $bed .bed).fq)
    file=MManalysis_$(basename $bed)
    python3 $pymatch $bed $reference $fastq -n $polyN -o $resultsloc/MManalysis
    done
    
}

filter_matches() {
	
    for file in $(ls $resultsloc/MManalysis*)
    do
    cat $file | awk BEGIN'{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,substr($4,length($4),1),substr($5,length($5),1)}' > $resultsloc/poly_$(basename $file)
    cat $resultsloc/poly_$(basename $file) | awk -v OFS="\t" -F"\t" '{if($7==$8) {print $0} }' > $resultsloc/matches_$(basename $file)
    done

}

remove_polyN() {
    for file in $(ls $resultsloc/matches_*)
    do
	Rscript $RrmpolyN -f $file -c 4 -t $Ntail -n $polyN -o $resultsloc
    done
}

cal_mm_percent() {
    echo -e 'Library\tA\tC\tG\tT\tTotal' > $resultsloc/chrM_mismatch_percent.txt
	echo -e 'Library\tA\tC\tG\tT\tTotal' > $resultsloc/nucl_mismatch_percent.txt
	
	for bed in $(ls $bedloc/*.bed)
    do
		bedtools getfasta -tab -s -fi $reference -bed $bed | grep chrM > $resultsloc/chrM_intemp
		bedtools getfasta -tab -s -fi $reference -bed $resultsloc/*$(basename $bed .bed)*final.bed | grep chrM > $resultsloc/chrM_outtemp
		bedtools getfasta -tab -s -fi $reference -bed $bed | grep -v chrM > $resultsloc/nucl_intemp
		bedtools getfasta -tab -s -fi $reference -bed $resultsloc/*$(basename $bed .bed)*final.bed | grep -v chrM > $resultsloc/nucl_outtemp

		for N in A C G T; do
			before=$(cat $resultsloc/chrM_intemp | grep -ie $N$ | wc -l)
			after=$(cat $resultsloc/chrM_outtemp | grep -ie $N$ | wc -l)
			eval $(echo mm_$N)=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100| bc -l)
		done

		before=$(cat $resultsloc/chrM_intemp | wc -l)
		after=$(cat $resultsloc/chrM_outtemp| wc -l)
		mm_total=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100 | bc -l)
		echo -e $(basename $bed .bed)'\t'$mm_A'\t'$mm_C'\t'$mm_G'\t'$mm_T'\t'$mm_total >> $resultsloc/chrM_mismatch_percent.txt
	  
	
		for N in A C G T; do
			before=$(cat $resultsloc/nucl_intemp | grep -ie $N$ | wc -l)
			after=$(cat $resultsloc/nucl_outtemp | grep -ie $N$ | wc -l)
			eval $(echo mm_$N)=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100| bc -l)
		done
		before=$(cat $resultsloc/nucl_intemp | wc -l)
		after=$(cat $resultsloc/nucl_outtemp| wc -l)
		mm_total=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100 | bc -l)
		echo -e $(basename $bed .bed)'\t'$mm_A'\t'$mm_C'\t'$mm_G'\t'$mm_T'\t'$mm_total >> $resultsloc/nucl_mismatch_percent.txt
	
    done
	rm chrM_outtemp
	rm nucl_outtemp
}


