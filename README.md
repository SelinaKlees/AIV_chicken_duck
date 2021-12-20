# AIV_chicken_duck
Data and code for reproducing the analyses presented in Klees et al. "Comparative Investigation of Gene Regulatory Processes Underlying Avian Influenza Viruses in Chicken and Duck". The results are provided in the corresponding folders, but may be reproduced by following the instructions below. 

# Requirements 
R programming language (4.1.0) Packages: - DESeq2 (1.34.0) - gprofiler2 (0.2.1) - ashr (2.2-47) - biomaRt (2.50.1) - tidyverse (1.3.1) - data.table (1.14.2) 

# Reproducing the results
1. Download or clone the github repository to retain the folder structure. 
    * The experimental data may be found in Duck/Experiment and Chicken/Experiment.
    * Please refer to the publicly available resources on Expression Atlas for [Duck](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2909/Results) and [Chicken](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2908/Results).
2. Run the R script analysis.R. 
    * Once finished, data will be written to both folders Duck/ and Chicken/
    * DEG/ contains the significant DESeq results
    * DEG_orth/ contains the orthologous genes to the DEGs.
    * GO/ contains the GO terms for the DEGs (divided into up-and downregulated)
    * Promoter/Info/ contains a more detailed gene list for the DEGs
3. For each file in Promoter/Info/, the promoter regions of the genes (-1000bp to +100bp) are obtained in fasta format from the [Ensembl](https://www.ensembl.org/index.html) BioMart and put into Promoter/Fasta/FG/. 
    * For chicken, the reference genome GRCg6a was used, for duck CAU_duck1.0.
4. For the creation of the background set, the promoter regions of all genes that are not differentially
expressed in any condition are obtained from [Ensembl](https://www.ensembl.org/index.html) BioMart and put into the web interface of
[BiasAway](https://biasaway.uio.no/biasaway/g/) using all of the sequences obtained in the step above as foreground. For each individual set
of promoters of the foreground, a random sample is then taken from the result of BiasAway. 
    * Please also refer to the section “Identification of Enriched TFs and TF-TF Cooperations” from the paper for more information
5. Transcription factor enrichment is carried out by [CiiiDER](http://ciiider.org/) using the previously described foreground and background sequences.
6. Transcription factor collaboration analysis on the promoter sequences in Promoter/Fasta/FG is done by [PC-TraFF](http://pctraffpro.bioinf.med.uni-goettingen.de/index.html).
    * The extension of this algorithm, described in [PMID: 29896218](https://www.frontiersin.org/articles/10.3389/fgene.2018.00189/full)<!--**PMID: 31095599 ??**-->
    , was used where the results are further filtered by running the algorithm 100 times on permutations of the original sequences.
7. Upstream analysis is carried out on the genexplain platform.
    * The TF pairs obtained from the previous step are combined for lung and ileum in the conditions regarding the H5N1 infection for both duck and chicken.
    * The resulting list of TFs is analysed by genexplain’s upstream analysis.
    * Please also refer to the section “Identification of Master Regulators” in the methods of the paper.



