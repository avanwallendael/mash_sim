###bash script for aligning simulated polyploid reads for evaluation of the effect of missing data
##performed on the MSU HPCC cluster
##06/05/20

#load modules
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load SAMtools
module load BCFtools
module load BWA
module load VCFtools/0.1.15-Perl-5.28.0

#copy in the reference (r49.fa with no missing data)
cp ../ref.fa ./ref.fa

#index the ref
bwa index ref.fa

#align seqs to the ref
for i in *.fa; do bwa mem ref.fa ${i} > aln_${i}.sam; done

#generate pileup bcfs
for i in *.sam; do samtools mpileup -uf ref.fa ${i} > ${i}.bcf; done

#call snps
for i in *.bcf; do bcftools call ${i} -o ${i}.vcf -O v -c; done

#zip and index vcfs
for i in *.vcf; do bgzip ${i};done
for i in *.gz; do tabix -p vcf ${i}; done

#combine vcfs, very slow step
vcf-merge aln_r0_miss10.fa.sam.bcf.vcf.gz aln_r21_miss10.fa.sam.bcf.vcf.gz aln_r34_miss10.fa.sam.bcf.vcf.gz aln_r47_miss10.fa.sam.bcf.vcf.gz aln_r1_miss10.fa.sam.bcf.vcf.gz aln_r22_miss10.fa.sam.bcf.vcf.gz aln_r35_miss10.fa.sam.bcf.vcf.gz aln_r48_miss10.fa.sam.bcf.vcf.gz aln_r10_miss10.fa.sam.bcf.vcf.gz aln_r23_miss10.fa.sam.bcf.vcf.gz aln_r36_miss10.fa.sam.bcf.vcf.gz aln_r49_miss10.fa.sam.bcf.vcf.gz aln_r11_miss10.fa.sam.bcf.vcf.gz aln_r24_miss10.fa.sam.bcf.vcf.gz aln_r37_miss10.fa.sam.bcf.vcf.gz aln_r5_miss10.fa.sam.bcf.vcf.gz aln_r12_miss10.fa.sam.bcf.vcf.gz aln_r25_miss10.fa.sam.bcf.vcf.gz aln_r38_miss10.fa.sam.bcf.vcf.gz aln_r6_miss10.fa.sam.bcf.vcf.gz aln_r13_miss10.fa.sam.bcf.vcf.gz aln_r26_miss10.fa.sam.bcf.vcf.gz aln_r39_miss10.fa.sam.bcf.vcf.gz aln_r7_miss10.fa.sam.bcf.vcf.gz aln_r14_miss10.fa.sam.bcf.vcf.gz aln_r27_miss10.fa.sam.bcf.vcf.gz aln_r4_miss10.fa.sam.bcf.vcf.gz aln_r8_miss10.fa.sam.bcf.vcf.gz aln_r15_miss10.fa.sam.bcf.vcf.gz aln_r28_miss10.fa.sam.bcf.vcf.gz aln_r40_miss10.fa.sam.bcf.vcf.gz aln_r9_miss10.fa.sam.bcf.vcf.gz aln_r16_miss10.fa.sam.bcf.vcf.gz aln_r29_miss10.fa.sam.bcf.vcf.gz aln_r41_miss10.fa.sam.bcf.vcf.gz aln_r17_miss10.fa.sam.bcf.vcf.gz aln_r3_miss10.fa.sam.bcf.vcf.gz aln_r42_miss10.fa.sam.bcf.vcf.gz aln_r18_miss10.fa.sam.bcf.vcf.gz aln_r30_miss10.fa.sam.bcf.vcf.gz aln_r43_miss10.fa.sam.bcf.vcf.gz aln_r19_miss10.fa.sam.bcf.vcf.gz aln_r31_miss10.fa.sam.bcf.vcf.gz aln_r44_miss10.fa.sam.bcf.vcf.gz aln_r2_miss10.fa.sam.bcf.vcf.gz aln_r32_miss10.fa.sam.bcf.vcf.gz aln_r45_miss10.fa.sam.bcf.vcf.gz aln_r20_miss10.fa.sam.bcf.vcf.gz aln_r33_miss10.fa.sam.bcf.vcf.gz aln_r46_miss10.fa.sam.bcf.vcf.gz > outall.vcf.gz

#output is processed in R with polymiss_dist_github.R