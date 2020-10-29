###bash script for downloading sequence data from the NCBI sra
##performed on the MSU HPCC cluster
##10/25/20

module load SRA-Toolkit
module load mash

#download fastqs for Capsella from Cornille et al. 2016
for i in {262..522}; do fastq-dump SRR2918${i}; done

#check for low-quality reads
ll

#remove low outliers
rm SRR2918268.fastq
rm SRR2918499.fastq
rm SRR2918375.fastq
rm SRR2918503.fastq

#clip reads to even length to avoid mash biases. This dataset has identical heads on several files, so tailcrop is better than head.
for i in *.fastq; do tail -n 1000000 ${i} > ${i}.tailcrop.fastq; done

#sketch reads with minimum depth of 2, k of 21, sampling 100k kmers
for i in *crop.fastq; do mash sketch -m 2 -k 21 -s 100000 ${i}; done

#list all files, then use a list sketch to combine them

ls *.fastq.msh > list
mash sketch -l list

#find all pairwise mash distances and write to a tsv
mash dist list.msh list.msh > tbl1_capsellatrim_tailcrop.tab