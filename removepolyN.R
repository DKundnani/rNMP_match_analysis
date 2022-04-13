#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="file with matched rNMP with longer sequence output ",metavar="character"),
  make_option(c("-c", "--colnum"), type="integer", default=4, 
              help="column which has sequence more than 1 ribo",metavar="integer"),
  make_option(c("-t", "--tail"), type="character", default=NULL, 
              help="nucleotide used for dN tailing (N)",metavar="character"),
  make_option(c("-n", "--num"), type="integer", default=4, 
              help="number of poly N's to be filtered upstream of ribos",metavar="integer"),
  make_option(c("-o", "--output_folder"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Please specify annotated file.n", call.=FALSE)
}
if (is.null(opt$colnum)){
  print_help(opt_parser)
  stop("Please specify the column to be used to draw concensus.n", call.=FALSE)
}


#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")
suppressMessages(library(ggseqlogo, quietly = T))
suppressMessages(library(magrittr, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(stringr, quietly = T))



#2.Defining functions if any
writeLines("\n...Defining required functions...\n")
plot_meme <- function(df,i,s) {
    df=df[df$V1 != 'chrM',]  #Can be changeg!!!!!!!!!!!!!!!!!!!!!!!!!
    if (nrow(df)>5*10^5) {
        set.seed(12345)  
        df=df[sample(1:nrow(df), 5*10^5), ]
    }
    jpeg(paste(paste(opt$output_folder,"/",sep=""),paste(str_split(file, pattern="[.]")[[1]][1],"kmer",i,s,"fastq.jpeg", sep="_"), sep=""), width =5, height = 3, unit ='in', res=1000)
    g<-ggseqlogo(df[c], method = 'bits')
    print(g)
    dev.off()

    jpeg(paste(paste(opt$output_folder,"/",sep=""),paste(str_split(file, pattern="[.]")[[1]][1],"kmer",i,s,"ref.jpeg", sep="_"), sep=""), width =5, height = 3, unit ='in', res=1000)
    g<-ggseqlogo(df[c+1], method = 'bits')
    print(g)
    dev.off()
    
}

#3. Preprocessing input files
writeLines("\n...Processing Input files...\n")
file=opt$file
bed<-read.table(file,sep="\t",header=F)
bed$V5=toupper(bed$V5)
c=opt$colnum
tail=opt$tail
x=opt$num

#4.Main code
 #Limit of filtration
filter=bed
kmer1=paste(replicate(1,tail), collapse = "")
filter=mutate(filter,V7 = substr(V4, nchar(V4)-1+1, nchar(V4)))
filter=mutate(filter,V8 = substr(V5, nchar(V5)-1+1, nchar(V5)))

for (i in seq(1,x)) {
    kmer=paste(replicate(i,tail), collapse = "")
    filter=mutate(filter,V9 = substr(V4, nchar(V4)-i+1, nchar(V4)))
    filter=mutate(filter,V10 = substr(V5, nchar(V5)-i+1, nchar(V5)))
    kmer.df=filter[filter$V9==kmer,]
    if (i>1) {
        m_kmer=kmer.df[kmer.df$V9==kmer.df$V10,]
        status='match'
        #plot_meme(m_kmer,i-1,status)
        print(paste(status, i,nrow(m_kmer),sep=" "))

        mm_kmer=kmer.df[kmer.df$V9!=kmer.df$V10,]
        status='mismatch'
        #plot_meme(mm_kmer,i-1,status)
        print(paste(status, i,nrow(mm_kmer),sep=" "))
        filter=filter[!(rownames(filter) %in% rownames(mm_kmer)),]
    } 
    status='filtered'
    #plot_meme(filter,i-1,status)
    print(paste(status,i,nrow(filter),sep=" "))
    status='filtered_rA'
    kmer1.df=filter[filter$V7==kmer1,]
    #plot_meme(kmer1.df,i-1,status)
    print(paste(status, i,nrow(kmer1.df),sep=" "))
}

write.table(filter[1:6], paste(str_split(file, pattern="[.]")[[1]][1],kmer,i-1,"final.bed", sep="_"), sep='\t',col.names=F, row.names=F, quote=F)

#5.Plotting

writeLines("\nJob Finished.\n")
writeLines(paste("Output Folder:", paste(opt$output_folder,"/", sep="")))
