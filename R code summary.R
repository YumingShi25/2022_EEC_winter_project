##popoolation2 snp frequency result processing and allele frequency calculation
setwd("~")
snp_freq <- read.delim("S1828Nr1-7-diff_rc.txt")
library(dplyr)
library(tidyverse)

#select diallelic
snp_freq <- snp_freq[snp_freq$allele_count == 2,]

#show major/minor allele counts and total counts
snp_freq <- snp_freq %>% separate(maa_1,c("major_count1","total_count1"),"/")
snp_freq <- snp_freq %>% separate(maa_2,c("major_count2","total_count2"),"/")
snp_freq <- snp_freq %>% separate(maa_3,c("major_count3","total_count3"),"/")
snp_freq <- snp_freq %>% separate(maa_4,c("major_count4","total_count4"),"/")
snp_freq <- snp_freq %>% separate(maa_5,c("major_count5","total_count5"),"/")
snp_freq <- snp_freq %>% separate(maa_6,c("major_count6","total_count6"),"/")
snp_freq <- snp_freq %>% separate(maa_7,c("major_count7","total_count7"),"/")
snp_freq <- snp_freq %>% separate(mia_1,c("minor_count1",NA),"/")
snp_freq <- snp_freq %>% separate(mia_2,c("minor_count2",NA),"/")
snp_freq <- snp_freq %>% separate(mia_3,c("minor_count3",NA),"/")
snp_freq <- snp_freq %>% separate(mia_4,c("minor_count4",NA),"/")
snp_freq <- snp_freq %>% separate(mia_5,c("minor_count5",NA),"/")
snp_freq <- snp_freq %>% separate(mia_6,c("minor_count6",NA),"/")
snp_freq <- snp_freq %>% separate(mia_7,c("minor_count7",NA),"/")
snp_freq$major_count1<-as.numeric(snp_freq$major_count1)
snp_freq$major_count2<-as.numeric(snp_freq$major_count2)
snp_freq$major_count3<-as.numeric(snp_freq$major_count3)
snp_freq$major_count4<-as.numeric(snp_freq$major_count4)
snp_freq$major_count5<-as.numeric(snp_freq$major_count5)
snp_freq$major_count6<-as.numeric(snp_freq$major_count6)
snp_freq$major_count7<-as.numeric(snp_freq$major_count7)
snp_freq$minor_count1<-as.numeric(snp_freq$minor_count1)
snp_freq$minor_count2<-as.numeric(snp_freq$minor_count2)
snp_freq$minor_count3<-as.numeric(snp_freq$minor_count3)
snp_freq$minor_count4<-as.numeric(snp_freq$minor_count4)
snp_freq$minor_count5<-as.numeric(snp_freq$minor_count5)
snp_freq$minor_count6<-as.numeric(snp_freq$minor_count6)
snp_freq$minor_count7<-as.numeric(snp_freq$minor_count7)
snp_freq$total_count1<-as.numeric(snp_freq$total_count1)
snp_freq$total_count2<-as.numeric(snp_freq$total_count2)
snp_freq$total_count3<-as.numeric(snp_freq$total_count3)
snp_freq$total_count4<-as.numeric(snp_freq$total_count4)
snp_freq$total_count5<-as.numeric(snp_freq$total_count5)
snp_freq$total_count6<-as.numeric(snp_freq$total_count6)
snp_freq$total_count7<-as.numeric(snp_freq$total_count7)

#show reference allele and alternative
snp_freq <- snp_freq %>% separate(allele_states,c("reference","SNP"),"/")

#substrate allele used for popoolation2 calculation in each population
snp_freq$major_allele1 <- substr(snp_freq$major_alleles.maa.,1,1)
snp_freq$major_allele2 <- substr(snp_freq$major_alleles.maa.,2,2)
snp_freq$major_allele3 <- substr(snp_freq$major_alleles.maa.,3,3)
snp_freq$major_allele4 <- substr(snp_freq$major_alleles.maa.,4,4)
snp_freq$major_allele5 <- substr(snp_freq$major_alleles.maa.,5,5)
snp_freq$major_allele6 <- substr(snp_freq$major_alleles.maa.,6,6)
snp_freq$major_allele7 <- substr(snp_freq$major_alleles.maa.,7,7)
snp_freq$minor_allele1 <- substr(snp_freq$minor_alleles.mia.,1,1)
snp_freq$minor_allele2 <- substr(snp_freq$minor_alleles.mia.,2,2)
snp_freq$minor_allele3 <- substr(snp_freq$minor_alleles.mia.,3,3)
snp_freq$minor_allele4 <- substr(snp_freq$minor_alleles.mia.,4,4)
snp_freq$minor_allele5 <- substr(snp_freq$minor_alleles.mia.,5,5)
snp_freq$minor_allele6 <- substr(snp_freq$minor_alleles.mia.,6,6)
snp_freq$minor_allele7 <- substr(snp_freq$minor_alleles.mia.,7,7)

#fix polarity
# If the new reference character (most common we had across pops) MATCHES the major character, then we set our new majorC1 variable
# to the majorcounts1 variable we already had (i.e. nothing changes)
snp_freq$majorC6 <- snp_freq$majorC5 <- snp_freq$majorC4 <- snp_freq$majorC3 <- snp_freq$majorC2 <- snp_freq$majorC2 <- 0
ifelse(snp_freq$reference == snp_freq$major_allele1, snp_freq$majorC1 <- snp_freq$major_count1,snp_freq$majorC1 <- snp_freq$minor_count1)
ifelse(snp_freq$reference == snp_freq$major_allele2, snp_freq$majorC2 <- snp_freq$major_count2,snp_freq$majorC2 <- snp_freq$minor_count2)
ifelse(snp_freq$reference == snp_freq$major_allele3, snp_freq$majorC3 <- snp_freq$major_count3,snp_freq$majorC3 <- snp_freq$minor_count3)
ifelse(snp_freq$reference == snp_freq$major_allele4, snp_freq$majorC4 <- snp_freq$major_count4,snp_freq$majorC4 <- snp_freq$minor_count4)
ifelse(snp_freq$reference == snp_freq$major_allele5, snp_freq$majorC5 <- snp_freq$major_count5,snp_freq$majorC5 <- snp_freq$minor_count5)
ifelse(snp_freq$reference == snp_freq$major_allele6, snp_freq$majorC6 <- snp_freq$major_count6,snp_freq$majorC6 <- snp_freq$minor_count6)
ifelse(snp_freq$reference == snp_freq$major_allele7, snp_freq$majorC7 <- snp_freq$major_count7,snp_freq$majorC7 <- snp_freq$minor_count7)

#calculate minor count
snp_freq$minorC1 <- snp_freq$total_count1 - snp_freq$majorC1
snp_freq$minorC2 <- snp_freq$total_count2 - snp_freq$majorC2
snp_freq$minorC3 <- snp_freq$total_count3 - snp_freq$majorC3
snp_freq$minorC4 <- snp_freq$total_count4 - snp_freq$majorC4
snp_freq$minorC5 <- snp_freq$total_count5 - snp_freq$majorC5
snp_freq$minorC6 <- snp_freq$total_count6 - snp_freq$majorC6
snp_freq$minorC7 <- snp_freq$total_count7 - snp_freq$majorC7

#calculate allele frequency
snp_freq$major_frequency1 <- snp_freq$majorC1/snp_freq$total_count1
snp_freq$major_frequency2 <- snp_freq$majorC2/snp_freq$total_count2
snp_freq$major_frequency3 <- snp_freq$majorC3/snp_freq$total_count3
snp_freq$major_frequency4 <- snp_freq$majorC4/snp_freq$total_count4
snp_freq$major_frequency5 <- snp_freq$majorC5/snp_freq$total_count5
snp_freq$major_frequency6 <- snp_freq$majorC6/snp_freq$total_count6
snp_freq$major_frequency7 <- snp_freq$majorC7/snp_freq$total_count7

#filtering threshold = 5%
five <-which(snp_freq$major_frequency1 >= .05 & snp_freq$major_frequency1 <= .95 |
                      snp_freq$major_frequency2 >= .05 & snp_freq$major_frequency2 <= .95 | 
                      snp_freq$major_frequency3 >= .05 & snp_freq$major_frequency3 <= .95 |
                      snp_freq$major_frequency4 >= .05 & snp_freq$major_frequency4 <= .95 |
                      snp_freq$major_frequency5 >= .05 & snp_freq$major_frequency5 <= .95 |
                      snp_freq$major_frequency6 >= .05 & snp_freq$major_frequency6 <= .95)
snp_freq <- snp_freq[five,]

#output
write.table(snp_freq, file = "snp_freq_filtered.txt")


##fst value filtering
setwd("~")
fst <- read.table("S1828Nr1-7fst.txt")

library(stringr)
#substrate fst values
for (i in 6:26) {
  colnames(fst)[i] <- substr(fst[1,i],1,3)
  fst[,i] <- as.numeric(str_split(fst[,i],"=",simplify = TRUE)[,2])
}

colnames(fst)[1:5] <- c("chr","pos","SNPs number","fraction of window","ave.min.coverage") 
#subset chr, pos, fst values
fst <- fst[,c(1:2,5:26)]

#fst between ancestral and evolved should set a threshold on 0.05
#filtering fst between ancestral and evolved population, fst >= 0.05 as the symbol of evolution 
fst.evo.sig <- fst[which(fst$'1:2' >= 0.05 & fst$'1:3' >= 0.05 & fst$`1:4` >= 0.05
                         & fst$'1:5' >= 0.05 & fst$'1:6' >= 0.05 & fst$`1:7` >= 0.05),]
#filtering fst between evolved populations, fst <= 0.05 as the symbol of common
fst.evo.com <- fst.evo.sig[which(fst.evo.sig$'2:3' < 0.05 & fst.evo.sig$'2:4' < 0.05 & fst.evo.sig$`2:5` < 0.05
                                 & fst.evo.sig$'2:6' < 0.05 & fst.evo.sig$'2:7' < 0.05 & fst.evo.sig$`3:4` < 0.05
                                 & fst.evo.sig$'3:5' < 0.05 & fst.evo.sig$'3:6' < 0.05 & fst.evo.sig$`3:7` < 0.05
                                 & fst.evo.sig$'4:5' < 0.05 & fst.evo.sig$'4:6' < 0.05 & fst.evo.sig$`4:7` < 0.05
                                 & fst.evo.sig$'5:6' < 0.05 & fst.evo.sig$'5:7' < 0.05 & fst.evo.sig$`6:7` < 0.05),]

write.table(fst.evo.com, file = "fst_common_snp.txt")


##Manhattan plot with fst values
setwd("~")
library(qqman)
#four cols: snp(SNP name), chr, bp(position), p(p-value) are needed
fst_common_snp <- read.table("fst_common_snp.txt", header=TRUE, quote="", comment.char="")
#rename columns
colnames(fst_common_snp) <- c("chr","pos","cov","fst12","fst13","fst14","fst15","fst16","fst17",
                              "fst23","fst24","fst25","fst26","fst27","fst34","fst35","fst36","fst37",
                              "fst45","fst46","fst47","fst56","fst57","fst67")
#replace character chromomsome numbers with numeric numbers
fst_common_snp$chr <- gsub('"2"',2,fst_common_snp$chr)
fst_common_snp$chr <- gsub('"1"',1,fst_common_snp$chr)
fst_common_snp$chr <- gsub('"3"',3,fst_common_snp$chr)
fst_common_snp <- fst_common_snp[which(fst_common_snp <= 3),]
fst_common_snp$chr <- as.numeric(fst_common_snp$chr)
#remove NA
fst_common_snp <- na.omit(fst_common_snp)

library(tidyverse)
#combine fst on the same position into one column
fst_common_snp <- fst_common_snp %>% pivot_longer(fst12:fst67,names_to = "populations",names_prefix = "fst",values_to = "fst")

#manhattan plot
man_test <- manhattan(fst_common_snp,chr = "chr", bp = "pos", p = "fst", snp = "pos",
                      logp = FALSE, col = c("blue","red","green"), ylab = "Fst")
#show the fst = 0.05 threshold
abline(h = 0.05,col="black",lwd=5)


##create vcf file for SNPEff analysis
setwd("~")

library(tidyr)
library(dplyr)

common_snp <- read.table("common_snp_ref.txt", header=FALSE, stringsAsFactors = FALSE)
common_snp <- common_snp %>% separate(V1,c("chrom","chromStart","chromEnd"))
common_snp <- common_snp[,c(1,2,4)]
#all reference alleles in the genome
write.table(common_snp, "common_snp_ref_mod.txt", col.names = FALSE, row.names = FALSE, sep="\t")

#alternative alleles(evolved alleles)
alt <- read.table("alt_snp.txt", header=TRUE, stringsAsFactors = FALSE)
alt <- alt[which(alt$chr == 1 | alt$chr == 2 | alt$chr == 3),]
alt$allele_states <- gsub('/', ',', alt$allele_states)
write.table(alt, "alt_snp_mod.txt", col.names = TRUE, row.names = FALSE, sep = "\t")

#find coresponding reference allele for every alternative allele
common_snp$com <- paste(common_snp$chrom,common_snp$chromStart,sep = ".")
alt$com <- paste(alt$chr,alt$pos,sep = ".")
overlap <- as.data.frame(alt[alt$com %in% common_snp$com,])
common_snp$allele_states <- overlap$allele_states

#modify as required for vcf file
common_snp$ID <- rep(".",nrow(common_snp))
common_snp$QUAL <- rep(".",nrow(common_snp))
common_snp$FILTER <- rep(".",nrow(common_snp))
common_snp$INFO <- rep(".",nrow(common_snp))
common_snp <- common_snp[,c(1,2,6,3,5,7,8,9)]
colnames(common_snp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
write.table(common_snp, "common_snp_ref.vcf", row.names=FALSE, sep="\t", quote = FALSE)


##generating gene ontology term list for Gowinda
setwd("~")

GO_TERM <- read.delim("Aedes_aegypti_gene_ontology_term.txt")
GO_TERM <- GO_TERM[which(GO_TERM$Computed.GO.Component.IDs!="N/A"),]
GO_TERM <- GO_TERM[,c(6,7,1)]

colnames(GO_TERM) <- c("GOID","GOTERM","GENEID")
rownames(GO_TERM) <- c(1:nrow(GO_TERM))
GO_TERM$GOID <- as.character(GO_TERM$GOID)
GO_TERM$GOTERM <- as.character(GO_TERM$GOTERM)
GO_TERM$GENEID <- as.character(GO_TERM$GENEID)
#separate multiple GO IDs for one gene into several rows
for (i in 1:nrow(GO_TERM)) {
  if (grepl(";",GO_TERM$GOID[i])==TRUE) {
    temp_GOID <- data.frame(strsplit(GO_TERM$GOID[i],";"))
    temp_GOTERM <- data.frame(strsplit(GO_TERM$GOTERM[i],";"))
    temp_GENEID <- data.frame(Gene.ID <- rep(GO_TERM$GENEID[i],times= nrow(temp_GOID)))
    temp <- cbind(temp_GOID,temp_GOTERM,temp_GENEID)
    colnames(temp) <- c("GOID","GOTERM","GENEID")
    rownames(temp) <- paste(i,rownames(temp),sep = "-")
    GO_TERM <- rbind(GO_TERM,temp)
    rm(temp_GOID)
    rm(temp_GENEID)
    rm(temp_GOTERM)
  }
}
#remove blank space
GO_TERM <- GO_TERM[-c(which(grepl(";",GO_TERM$GOID[1:nrow(GO_TERM)])==TRUE)),]

#combine genes for the same GO term
for (i in 1:nrow(GO_TERM)) {
  for (j in (i+1):nrow(GO_TERM)) {
    if (isTRUE(GO_TERM$GOID[i]==GO_TERM$GOID[j])) {
      GO_TERM$GENEID[i] <- paste(GO_TERM$GENEID[i],GO_TERM$GENEID[j],sep = " ")
      GO_TERM <- GO_TERM[-j,]
    }
  }
}

write.table(GO_TERM,"Ae.aegypti_gene_set_noquotes.txt",row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)


##analysis of Gowinda result
setwd("~")
GO_result <- read.delim("result_high_resolution.txt", header=FALSE)
colnames(GO_result) <- c("ID","count_per_sim","count","pvalue","FDR","gene_count","gene_count_cat","total_gene","Description","geneID")

#Aedes aegypti reference genome
Ae.aegypti_gene_set <- read.delim("Ae.aegypti_gene_set_noquotes.txt")
Aedes_aegypti_lvpagqg.AaegL5.47.CDS_EXON3 <- read.delim("Aedes_aegypti_lvpagqg.AaegL5.47.CDS_EXON3.gtf", header=FALSE)

#create .db reference file for Aedes aegypti, two lists are needed: chromosome, go
library(AnnotationForge)
#chromosome = geneID + chromosome
aachr <- Aedes_aegypti_lvpagqg.AaegL5.47.CDS_EXON3[,c(1,9)]
aachr$geneID <- substr(aachr$V9,9,18)
aachr <- aachr[,c(3,1)]
colnames(aachr) <- c("GID","CHROMOSOME")
aachr$GID <- unique(aachr$GID)
aachr <- unique(aachr)

#go = geneID + GOID + Evidence
aaGO <- Ae.aegypti_gene_set[,c(3,1)]
#separate multiple GO ID for one Gene ID into several rows
for (i in 1:nrow(aaGO)) {
  if (grepl(" ",aaGO$GENE_ID[i])==TRUE) {
    temp_GENEID <- data.frame(strsplit(aaGO$GENE_ID[i]," "))
    temp_GOID <- data.frame(rep(aaGO$GO_ID[i],times= nrow(temp_GENEID)))
    temp <- cbind(temp_GENEID,temp_GOID)
    colnames(temp) <- c("GENE_ID","GO_ID")
    rownames(temp) <- paste(i,rownames(temp),sep = "-")
    aaGO <- rbind(aaGO,temp)
    rm(temp_GOID)
    rm(temp_GENEID)
  }
}
#remove blank space
aaGO <- aaGO[-c(which(grepl(" ",aaGO$GENE_ID[1:nrow(aaGO)])==TRUE)),]
colnames(aaGO) <- c("GID","GO")
aaGO$EVIDENCE <- rep("IEA",nrow(aaGO))
aaGO <- unique(aaGO)

#create the Aedes aegypti package using AnnotationForge
makeOrgPackage(chromosome=aachr, go=aaGO,
               version="0.1",
               maintainer="Yuming Shi <yuming.shi.25@gmail.com>",
               author="Yuming Shi <yuming.shi.25@gmail.com>",
               outputDir = "~",
               tax_id="7159",
               genus="Aedes",
               species="aegypti",
               goTable="go")

#GO enrichment analysis result visualize
GO_result <- read.delim("result_high_resolution.txt", header=FALSE)
colnames(GO_result) <- c("ID","count_per_sim","count","pvalue","qvalue","gene_count","gene_count_cat","total_gene","Description","geneID")
#replace gene ID with capital letters to match reference genome
library(stringr)
GO_result$geneID <- str_replace_all(GO_result$geneID,",","/")
GO_result$geneID <- str_replace_all(GO_result$geneID,"aael","AAEL")

install.packages("~/org.Aaegypti.eg.db", repos=NULL, type = "source")
library(org.Aaegypti.eg.db)
library(rrvgo)

#ont = gene ontlogy category. One of c("BP", "MF", "CC")
#keytype = the key word using to find matches in the reference genome
#filtering with FDR <= 0.01
GO_sig <- GO_result[which(GO_result$qvalue <= 0.01),]
#create a matrix for similarity calculation
GO_sig_simMatrix_BP <- calculateSimMatrix(GO_sig$ID,
                                          orgdb = "org.Aaegypti.eg.db",
                                          keytype = "GID",
                                          ont = "BP",
                                          method = "Rel")
GO_sig_simMatrix_MF <- calculateSimMatrix(GO_sig$ID,
                                          orgdb = "org.Aaegypti.eg.db",
                                          keytype = "GID",
                                          ont = "MF",
                                          method = "Rel")
GO_sig_simMatrix_CC <- calculateSimMatrix(GO_sig$ID,
                                          orgdb = "org.Aaegypti.eg.db",
                                          keytype = "GID",
                                          ont = "CC",
                                          method = "Rel")
#similarity score
GO_sig_scores <- setNames(-log10(GO_sig$qvalue), GO_sig$ID)

#reduced GO term list
#threshold = similarity threshold,  0.4 tiny, 0.5 small, 0.7 medium, large 0.9
GO_sig_reducedTerms_BP <- reduceSimMatrix(GO_sig_simMatrix_BP,
                                          GO_sig_scores,
                                          threshold = 0.5,
                                          orgdb="org.Aaegypti.eg.db",
                                          keytype = "GID")
GO_sig_reducedTerms_MF <- reduceSimMatrix(GO_sig_simMatrix_MF,
                                          GO_sig_scores,
                                          threshold = 0.5,
                                          orgdb="org.Aaegypti.eg.db",
                                          keytype = "GID")
GO_sig_reducedTerms_CC <- reduceSimMatrix(GO_sig_simMatrix_CC,
                                          GO_sig_scores,
                                          threshold = 0.5,
                                          orgdb="org.Aaegypti.eg.db",
                                          keytype = "GID")

#combine reduced term lists into one list
GO_sig_term_BP <- GO_sig_reducedTerms_BP
GO_sig_term_BP$category <- c("BP")
GO_sig_term_BP<- GO_sig_term_BP[order(-GO_sig_term_BP$size),]
GO_sig_term_MF <- GO_sig_reducedTerms_MF
GO_sig_term_MF$category <- c("MF")
GO_sig_term_MF<- GO_sig_term_MF[order(-GO_sig_term_MF$size),]
GO_sig_term_CC <- GO_sig_reducedTerms_CC
GO_sig_term_CC$category <- c("CC")
GO_sig_term_CC<- GO_sig_term_CC[order(-GO_sig_term_CC$size),]
GO_sig_term <- rbind(GO_sig_term_BP,GO_sig_term_MF,GO_sig_term_CC)
#subset useful info
GO_sig_term_s <- GO_sig_term[,c(1,5,6,7,9)] 

#GO category bar plot
library(ggplot2)
GO_sig_term_s$term <- factor(GO_sig_term_s$term, levels = GO_sig_term_s$term)
GO_bar <- ggplot(data = GO_sig_term_s, aes(x = term, y = log10(size), fill = category)) + 
  geom_bar(stat = "identity", colour = "black") + 
  coord_flip() + 
  labs(x = "Gene Ontology Terms") +
  ylab(bquote(log[10]("number of genes"))) + #type in subscript
  scale_fill_discrete(name = "Gene Ontology Category", 
                      labels = c("Biological Process", "Cell Component", "Molecular Function")) +
  theme(panel.background = element_rect(fill = "transparent"), # remove background of the panel
        plot.background = element_rect(fill = "transparent", color = NA),# remove background of the plot
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(size=10, face="bold"))

#plots provided in rrvgo
options(ggrepel.max.overlaps = Inf) #avoid label overlapping
#scatter plot(bubble plot)
scatterPlot(GO_sig_simMatrix_BP, GO_sig_reducedTerms_BP, size = "score",labelSize = 4)
#tree map
treemapPlot(GO_result_reducedTerms, size = "score")
#word clouds
wordcloudPlot(GO_result_reducedTerms, onlyParents = FALSE, min.freq=1, colors="black")


##SNPEff result pie chart
setwd("~")
SNPEff_result_region <- read.delim("SNPEff_result_region.txt")

library(ggplot2)
library(ggrepel)
library(dplyr)
#calculate position for labels
SNPEff_result_region <- SNPEff_result_region %>% arrange(desc(Type)) %>%
  mutate(prop = Count/sum(Count) *100) %>% mutate(ypos = cumsum(prop)- 0.5*prop)
#round the proportion with 2 digits
SNPEff_result_region$prop <- format(round(SNPEff_result_region$prop, 2))
SNPEff_result_region$prop <- as.numeric(SNPEff_result_region$prop)

#pie chart = curving bar chart
region_pie <- ggplot(SNPEff_result_region, aes(x = "", y = prop, fill = Type)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) + 
  theme_void() + 
  geom_label_repel(aes(y = ypos, label = paste0(Type, " ", prop, "%")), 
                   size = 5, nudge_x = 0.7, show.legend = FALSE, force = 10, force_pull = 5, xlim = 1000)



