#Merge realigned sequencing files of seven populations into one mpileup file using Samtools
samtools mpileup -B S1828Nr1_sort_rmdup_rg_realigned.bam \ 
                                  S1828Nr2_sort_rmdup_rg_realigned.bam \   
                                  S1828Nr3_sort_rmdup_rg_realigned.bam \ 
                                  S1828Nr4_sort_rmdup_rg_realigned.bam \         
               		  S1828Nr5_sort_rmdup_rg_realigned.bam \ 
                                  S1828Nr6_sort_rmdup_rg_realigned.bam \
                                  S1828Nr7_sort_rmdup_rg_realigned.bam > <working directory>/ S1828Nr1-7_mp.mpileup

#Convert mpileup file into sync file using mpileup2sync command in PoPoolation2
perl <popoolation2-path>/mpileup2sync.pl --fastq-type sanger --min-qual 20 --input S1828Nr1-7_mp.mpileup --output S1828Nr1-7_java_filtered.sync  --threads 8

#Calculate allele frequency for every SNP found in the genome using snp-frequency-diff command in PoPoolation2
perl <popoolation2-path>/snp-frequency-diff.pl --input S1828Nr1-7_java_filtered.sync --output-prefix pop-diff --min-count 14 --min-coverage 50 --max-coverage 500

#Calculate pairwise FST using fst-sliding command in PoPoolation2
perl <popoolation2-path>/fst-sliding.pl --input S1828Nr1-7_java_filtered.sync --output S1828Nr1-7.fst --min-count 14 --min-coverage 50 --max-coverage 500 --window-size 500000 --step-size 500000 --pool-size 80

#Get reference alleles for every SNP position from Aedes aegypti reference genome using getfasta command in bedtools
bedtools getfasta -fi Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa -bed common_snp.bed -tab -fo common_snp_ref.txt   

#Retrieve chromosome, position, and allele status from SNP frequency files 
cut -f1,2,5 S1828Nr1-7-diff_rc.txt > alt_snp.txt

#Download Aedes aegypti reference genome using SNPEff
java -jar <SNPEff path>/snpEff.jar download -v Aedes_aegypti_lvpagwg
                                       
#Add genome annotation to common evolved alleles using SNPEff
java -Xmx8g -jar <SNPEff path>/snpEff.jar Aedes_aegypti_lvpagwg common_snp_ref_cap.vcf > common_snp.ann.vcf

#Gene ontology enrichment analysis using Gowinda
java -Xmx4g -jar <Gowinda path>/Gowinda-1.12.jar --snp-file all_snp_pos.txt --candidate-snp-file candidate_snp.txt --gene-set-file Ae.aegypti_gene_set_noquotes.txt --annotation-file Aedes_aegypti_lvpagqg.AaegL5.47.CDS_EXON3.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 8 --output-file result_high_resolution.txt --mode snp --min-genes 1
