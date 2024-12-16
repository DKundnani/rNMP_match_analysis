#! usr/bin/bash
#Author: Deepali L. Kundnani
#Institute: Georgia institute of Technology(GT)
#Description: Filters Ribosemap generated bedfile from polyN found in ribose-seq sequencing data(upstream / 5' end from ribo) that doesn't match to the corresponding ribo position in reference genome (N=nucletoide used for dN tailing)

#pymatch='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/bin/rNMP_match_analysis/match_analysis.py' #match analysis script
#reference='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/reference/puc18/puc18.fa' #Fasta file for reference genome
#genome=$(echo $reference.fai) #This assumes you have respective genome.sizes in the same folder with etx as .fai to the fasta file
#polyN='5' # # number of polyN's to check including hte ribo position
#bedloc='/storage/coda1/p-fstorici3/0/shared/bed/temp' #location with bed files - Make sure you have the bed files only required for analysis, or it will analyze all the files
#fastqloc='/storage/coda1/p-fstorici3/0/shared/raw_reads/' #Make sure the basename of files to be analysed is same as the basename of '.fq' files. Sample.bed used Sample.fq for analysis
#resultsloc='/storage/coda1/p-fstorici3/0/shared/mm_filtered_bed/puc18' #location you want the results in 
#RrmpolyN='/storage/coda1/p-fstorici3/0/dkundnani3/rich_project_bio-storici/bin/rNMP_match_analysis/removepolyN.R' #removepolyN.R script 
#Ntail='A' #Nucleotide used for dN tailing (Default-'A')

mismatch_analysis() {
	mkdir $resultsloc
    for bed in $(ls $bedloc/*.bed)
    do
    fastq=$(echo $fastqloc/$(basename $bed .bed).fq)
    #python3 $pymatch $bed $reference $fastq -d 0 -n $polyN -o $resultsloc/MM &
	python3 $pymatch $bed $reference $fastq -d 0 -n $polyN -o $resultsloc/MM & #element bio
    done
    wait
}

single_filter_matches() {
    cat $1 | awk BEGIN'{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,substr($4,length($4),1),substr($5,length($5),1)}' > $resultsloc/poly_$(basename $1)
	cat $resultsloc/poly_$(basename $1) | awk -v OFS="\t" -F"\t" '{if($7!=$8) {print $0} }' > $resultsloc/mismatches_$(basename $1)
    cat $resultsloc/poly_$(basename $1) | awk -v OFS="\t" -F"\t" '{if($7==$8) {print $0} }' > $resultsloc/matches_$(basename $1)
	for N in A C G T; do
	awk -v n="$N" -F"\t" '{if($7==n) {print $0} }' $resultsloc/matches_$(basename $file) > $resultsloc/matches/${N}_matches_$(basename $file)
	awk -v n="$N" -F"\t" '{if($7==n) {print $0} }' $resultsloc/mismatches_$(basename $file) > $resultsloc/mismatches/${N}_mismatches_$(basename $file)
	done
}

filter_matches() {
    for file in $(ls $resultsloc/MM_*)
    do
		mkdir $resultsloc/matches $resultsloc/mismatches
		single_filter_matches $file &
    done
	wait
}

remove_polyN() {
	conda activate r_env
	for file in $(ls $resultsloc/matches_*)
    do
		Rscript $RrmpolyN -f $file -c 4 -t $Ntail -n $polyN -o $resultsloc &
    done
	wait
   	mkdir $resultsloc/final
   	cp $resultsloc/*final* $resultsloc/final/
	rename 'matches_MM_' '' $resultsloc/final/*
	rename '_final.bed' '' $resultsloc/final/*
}

single_rNMP_mmanalysis() {
	fastq=$(echo $fastqloc/$(basename $1 .bed).fq)
    python3 $pymatch $1 $reference $fastq -n 1 -o $resultsloc/single
    python3 $(dirname $pymatch)/count_matched_bed.py $1 -o $resultsloc/$(basename $1 .bed).counts
	python3 $(dirname $pymatch)/draw_match_heatmap.py $resultsloc/$(basename $1 .bed).counts --hide_cbar -o $resultsloc/

}

rNMP_mmvis() {
    for bed in $(ls $bedloc/*.bed)
    do
		single_rNMP_mmanalysis $bed &
    done
	wait

}

mkdir temp
cal_mm_percent(){
    #echo -e 'Library\tA\tC\tG\tT\tTotal' > temp/chrM_mismatch_percent.txt
	#echo -e 'Library\tA\tC\tG\tT\tTotal' > temp/nucl_mismatch_percent.txt

	for bed in $(ls $bedloc/*.bed)
    do
    echo $bed
    inbed=$bed
    #outbed=$resultsloc/*$(basename $bed .bed)*final.bed
	outbed=$resultsloc/matches_MM_*$(basename $bed) #element bio
		cal_single_mm_percent ${1}
    rm temp/*temp
   done
   wait

}

cal_single_mm_percent() {
	bedtools getfasta -tab -s -fi $reference -bed $inbed | grep ${1} > temp/${1}_intemp
	bedtools getfasta -tab -s -fi $reference -bed $outbed | grep ${1} > temp/${1}_outtemp
	bedtools getfasta -tab -s -fi $reference -bed $inbed | grep -v ${1} > temp/nucl_intemp
	bedtools getfasta -tab -s -fi $reference -bed $outbed | grep -v ${1} > temp/nucl_outtemp
	for N in A C G T; do
		before=$(cat temp/${1}_intemp | grep -ie $N$ | wc -l)
		after=$(cat temp/${1}_outtemp | grep -ie $N$ | wc -l)
		eval $(echo mm_$N)=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100| bc -l)
	done
	before=$(cat temp/${1}_intemp | wc -l)
	after=$(cat temp/${1}_outtemp| wc -l)
	mm_total=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100 | bc -l)
	echo -e $(basename $inbed .bed)'\t'$mm_A'\t'$mm_C'\t'$mm_G'\t'$mm_T'\t'$mm_total >> temp/${1}_mismatch_percent.txt
	for N in A C G T; do
		before=$(cat temp/nucl_intemp | grep -ie $N$ | wc -l)
		after=$(cat temp/nucl_outtemp | grep -ie $N$ | wc -l)
		eval $(echo mm_$N)=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100| bc -l)
	done
	before=$(cat temp/nucl_intemp | wc -l)
	after=$(cat temp/nucl_outtemp| wc -l)
	mm_total=$(echo $(echo $(echo $before-$after|bc -l)/$before | bc -l)*100 | bc -l)
	echo -e $(basename $inbed .bed)'\t'$mm_A'\t'$mm_C'\t'$mm_G'\t'$mm_T'\t'$mm_total >> temp/nucl_mismatch_percent.txt
}
