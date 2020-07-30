# GENIE3_scripts
Scripts used in the manuscript "Validation and characterisation of the wheat GENIE3 network using an independent RNA-Seq dataset." from Harrington et al. 2020


ExtractCDS_fromPearceDEGs.py

  This script extracts the coding sequence for the differentially expressed (DE) genes identified in the Pearce et al. 2014 analysis (https://doi.org/10.1186/s12870-014-0368-2)
  In brief, it extracts the sequence identified as DE from the URGI contig provided.
  
BLAST_Pearce_cds_against_HC_and_LC.sh

  This script uses BLASTn to identify the gene model in the RefSeqv1.1 annotation that corresponds with the URGI-derived sequence obtained using the previous script.

Genie3_Statistics_SharedRatios_RNASeqDEGs_Figure1.Rmd

  This script was used to produce the distributions seen in Figure 1.
  In brief, this identifies the level of overlap ("shared ratio") between the targets of randomly-selected transcription factors in the GENIE3 network and the DE genes identified in the reanalysis of the Pearce et al. 2014 data.
  This also demonstrates how to pull out the targets of a given transcription factor from the GENIE3 network.
  
Genie3_Statistics_SharedRatios_AllTFs_AllTFFamilies_Fig3_SuppFig2.Rmd

  This script was used to produce the graphs in Figure 3 and Supplementary Figure 2.
  In brief, it calculates the "shared ratio" between randomly selected transcription factors, and transcription factors in a particular family.
  
Genie3_Statistics_SharedRatios_Homoeologs_MovementCategories_Figure4_SuppFig3.Rmd

  This script was used to produce the graphs in Figure 4 and Supplementary Figure 2, as well as the homoeolog distribution in Figure 3.
  In brief, it calculates the "shared ratio" between homoeologs, and then investigates how triad movement category (stable/mid80/dynamic) affects the shared ratio between homoeologs.
  It also looks at whether the genome of a transcription factor influences the distribution of targets across the three genomes (A/B/D).

GO_enrichment.R
  
   This script carried out GO term enrichment on the DE genes, and provides an R function which can be utilised to carry out GO enrichment on any desired subset of genes, given that gene set and the corresponding gene universe.
   
Sleuth_KO_12DAA_WT_ref.R

  This script exemplifies our use of the Sleuth differential expression software in re-analysing the RNA-Seq data. We provide this script for completeness, however encourage users to first visit https://pachterlab.github.io/sleuth/ for information regarding the use of Sleuth.
