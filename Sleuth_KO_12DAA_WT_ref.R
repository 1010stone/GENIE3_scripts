##In this script, we are carrying out Differential Expression analysis on RNA-Seq data which had previously been pseudo-aligned against the wheat transcriptome using Kallisto.
##Here we show how we carried out the DE analysis for the 12 DAA timepoint. The same procedure was also used for 22 DAA.

##Please see Harrington et al. 2020 "The wheat GENIE3 network provides biologically-relevant information in polyploid wheat" for details on the transcriptome and Kallisto settings used.

##We recommend reading the primer on Sleuth before carrying out analysis, as that explains well the requirements for metadata formatting, etc. See https://pachterlab.github.io/sleuth/

library("sleuth")
library("dplyr")

wd <- file.path('//jic-hpc-data/Research-Groups/Cristobal-Uauy/Sophie/NAM-B1_DEAnalysis/TILLING_KO')

sample_id <- dir(file.path(wd,"Kallisto"))
kal_dirs <- file.path(wd,"Kallisto",sample_id,"30")

##filter the files to be included in the sleuth analysis by editing the overall metadata table (in this example, we exclude timepoints that we are not interested in using dplyr::filter)
s2c <- read.table(file.path(wd,"Sample_Metadata_A.txt"), header=TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::select(s2c, sample, genotype, stage)
s2c <- dplyr::mutate(s2c, path=kal_dirs)
s2c <- dplyr::filter(s2c, stage != 'HD')
s2c <- dplyr::filter(s2c, stage != '22DAA')

##here we run the LRT model of sleuth, comparing the expression of genes by the "genotype" category in the metadata
so <- sleuth_prep(s2c, ~genotype, num_cores=1, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~genotype, 'full')
so <- sleuth_fit(so, ~genotype, 'reduced')
so <- sleuth_lrt(so, 'reduced','full')

##here we run the Wald test sleuth model, saving this separately to the LRT model
so_a1 <- sleuth_wt(so, 'genotypegpc-a1')
head(so_a1)

##now we can filter the results to extract only those genes which are Differentially expressed, based on our q-value cutoff.
test_table_a1 <- sleuth_results(so_a1, 'genotypegpc-a1')
sleuth_significant_a1 <- dplyr::filter(test_table_a1, qval <= 0.05)
write.table(sleuth_significant_a1, "c:/Users/backhaua/Documents/R/R_outputs/0.05_2/12DAA_WT_ref_a1.txt", sep=" ")

##or we can save all of the genes, without filtering based on q-value
table_a1 <- sleuth_results(so_a1, 'genotypegpc-a1')
write.table(table_a1, "c:/Users/backhaua/Desktop/12DAA_WT_vs_a1.txt", sep=" ")
head(table_a1)
ggplot(table_a1, aes(x=qval)) + geom_histogram() + geom_vline(xintercept=0.05)

##here we can save both the LRT and the WT models
saveRDS(so, file = "12DAA.rds")
saveRDS(so_a1, file = "12DAA_WT.rds")


