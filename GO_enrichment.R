##Get the enriched GO terms for each of the DE gene sets

library(sleuth)
library(dplyr)
library(tidyr)
library(goseq)
library(gridExtra)
library(ggplot2)


##read in and process the data for each timepoint. Note-- please replace the file names with the appropriate names for your sleuth output files!
##this script is based on the RDS files generated from sleuth as in the file "Sleuth_KO_12DAA_WT-ref.R" which are is available on https://github.com/Uauy-Lab/GENIE3_scripts


##First, identify the universe of genes from which your subset is drawn. In this case, we use all genes which were included in the sleuth DE analysis at each timepoint

so_12DAA <- readRDS("Z://Anna/Genie3_codes_AB/sleuth output 07_05_2020/12DAA.rds")
filtered_genes_12DAA <- as.data.frame(unique(so_12DAA$obs_norm_filt$target_id))
colnames(filtered_genes_12DAA) <- c("Gene")
filtered_genes_12DAA$Gene <- gsub("\\..*","",filtered_genes_12DAA$Gene)
filtered_genes_12DAA <- unique(filtered_genes_12DAA)

so_22DAA <- readRDS("Z://Anna/Genie3_codes_AB/sleuth output 07_05_2020/22DAA.rds")
filtered_genes_22DAA <- as.data.frame(unique(so_22DAA$obs_norm_filt$target_id))
colnames(filtered_genes_22DAA) <- c("Gene")
filtered_genes_22DAA$Gene <- gsub("\\..*","",filtered_genes_22DAA$Gene)
filtered_genes_22DAA <- unique(filtered_genes_22DAA)

#extract transcript median length. This step allows us to take transcript length into account when carrying out the GO enrichment analyses.
gene_medians <- unique(so_12DAA$obs_norm %>% select(target_id, eff_len))

##now get median of eff_len over all of the transcripts
gene_medians$Gene <- gsub("\\..*","",gene_medians$target_id)
gene_medians <- gene_medians %>% group_by(Gene) %>% summarise(Median = median(eff_len))


#read in the DE genes. Here we are using files that contain only genes which are considered DE, in our case with a Wald test q-value < 0.05.
#These tables are extracted from the sleuth object, again see relevant files on Github (i.e. Sleuth_KO_12DAA_WT-ref.R)

NAMA1_vsWT_12DAANF <- read.csv("Y://Sophie/Networks/Genie3_NAMcomparison/Anna_RNAseqoutput/12DAA_gpcA1_ref_WT_NOFILTERING.txt", sep=" ", header=TRUE)
NAMA1_vsWT_22DAA <- read.csv("Y://Sophie/Networks/Genie3_NAMcomparison/Anna_RNAseqoutput/22DAA_gpcA1_ref_WT.txt", sep=" ", header=TRUE)


# read in GO terms
# GO terms can be downloaded from the Gene Ontology website, using the most recent wheat GO release. See methods for generating GO terms in Ramirez-Gonzalez et al. 2018 https://science.sciencemag.org/content/361/6403/eaar6089

all_go <- read.csv("Y:\\expression_browser\\WGA\\data_tables\\Formatted_Triticum_aestivum_V1_PGSB_GOA_by_orthology.tsv", sep="\t")
all_go <- all_go[,c(1,2)]

#if necessary, ensure gene IDs correspond between the GO terms and the DE genes, as the GO terms may still be based on RefSeqv1.0 annotation, while DE genes are using the v1.1 annotation
all_go$Gene <- (gsub("01G", "02G", all_go$Gene))


wd <- "Y://Sophie/Networks/Genie3_NAMcomparison/RNASeq_DE_GOEnrichment"


##function for GO terms enrichment, and saving the output. Edit file names as desired.

GO <- function(gene_list, name, universe){
  ##remove transcript data
  gene_list <- gsub("\\..*","",gene_list)
  print(length(gene_list))
  
  ##prepare the universe etc based on the category of the DE genes
  all_go <- subset(all_go, Gene %in% universe$Gene)
  t1 <- subset(gene_medians, Gene %in% universe$Gene)
  gene.lens <- as.numeric(t1$Median)
  print(length(gene.lens))
  all_genes <- as.vector(universe$Gene)
  
  ##get vector of all genes with 0/1 based on presence or absence in gene subset
  gene.vector <- as.integer(all_genes%in%gene_list)
  print(length(gene.vector))
  names(gene.vector)=all_genes
  #now carry out the GOseq analysis, save the pwf graph
  pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
  png(file.path(wd,paste0(name,"median_pwf.png")))
  plotPWF(pwf)
  dev.off()
  GO.wall = goseq(pwf, gene2cat = all_go, test.cats ="GO:BP")
  # add new column with over represented GO terms padj
  GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
  GO.wall <- GO.wall[GO.wall$over_rep_padj <0.05,]
  write.table(GO.wall, file = paste0(wd,"/",name,"_GOseq_enrichment_median.txt", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE)
}

##run the GO enrichment analysis on the desired subset of genes, and against the appropriate universe of genes.

GO(NAMA1_vsWT_12DAANF$target_id, "NAMA1_vsWT_12DAANF", filtered_genes_12DAA)
GO(NAMA1_vsWT_22DAA$target_id, "NAMA1_vsWT_22DAA", filtered_genes_22DAA)


