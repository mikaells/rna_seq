#clean out existing stuff to start fresh
rm(list = ls())

#0. SETUP THE LIBRARIES

#Uncomment if you dont have the packages

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Rsubread")

#read in the libraries we need
library("Rsubread")


#####
#1. SETUP THE SEQUENCING FILES 
#####

#get the directory of the raw files
fastq_directory <- "fastqs_subsets//"
#get the adresses of the files we need
#remember that all samples have two files
reads1 <- list.files(path = file.path(fastq_directory), pattern = "*1.fastq$", full.names = TRUE)
reads2 <- list.files(path = file.path(fastq_directory), pattern = "*2.fastq$", full.names = TRUE)

#do we have all 12 files? Should return 12 for both
cat(paste("reads1 =",length(reads1),"\nreads2 =", length(reads2)))


#####
#2. MAKE THE INDEX TO ALIGN AGAINST
#####

#first, we make a folder/directory for the index
dir.create("indx")

#Then we build the actual index from the whole genome file
#and we will put it in the indx/-folder we just made
buildindex(basename="indx/NC_000913_WG", #the name of our index
           reference="NC_000913_WG.fasta", #what file to make it from?
           memory=2000, #max memory is surely enough for a smallish bacterial genome
           gappedIndex = F) #?

#####
#3. MAKE THE MAP BETWEEN POSITION AND GENE
#####

# now we map the reference
#The .gff3 file has all that info, so we will turn it into a SAF-file and use that
saf=flattenGTF(
  GTFfile =  "NC_000913.gff3", 
  GTF.featureType = "gene",
  GTF.attrType = "gene",
  # the option specifying the merging algorithm
  method = "merge")


#####
#4. ALIGN THE SEQUENCES TO INDEX
#####

#we now align our sequencing files to the reference
#e.g. each read pair will be aligned to the genome
#note that the alignments are saved as external files rather than as R-objects
#****
#LOOK AT THE nthreads option
#if your computer only has 4 threads, its gonna die with 12. Set it at your threads minus 1.
#****
#will be fast on a linux (perhaps on a mac too) and not on a windows (something weird with the 
#index-load on some machines)
#it is working though, note bam-files appearing in the subset-folder

#make a folder for the bam-files
dir.create("bam_out")

#then the command for alignment
alignInfo=align(index = "indx/NC_000913_WG",   #the index
                readfile1 = reads1,            #read pairs 1
                readfile2 = reads2,            #read pairs 2
                type = "rna",                  #we are doing rna
                input_format = "FASTQ",        #files are fastq-files, could also be fastq.gz
                output_format = "BAM",         #output will be BAM files
                PE_orientation = "fr",         #reads are [f]orward and [r]everse
                nthreads = 6,                  #number of threads
                output_file = gsub("_1.fastq",".bam",gsub("fastqs_subsets","bam_out",reads1)) 
)


#####
#5.  COUNT HOW THE SEQUENCES MAPPED TO THE GENES
#####

#lets get the list of BAM files
bam.files <- list.files("bam_out/", pattern = "bam$", full.names = T)

#Run the actual counting.
fc=featureCounts(files = bam.files,      #the bam-files we are counting 
                 annot.ext = saf,   #the external object we are annotating with
                 isPairedEnd = T)   #and our files are obviously paired


stats=fc$stat

features=fc$counts


