# odna-coevolution
Coevolution of mtDNA and ptDNA genome content

Requires R with libraries phytools, ape, phangorn, phylolm, nlme, lme4, pheatmap, Oncotree, igraph, ggplot2, ggpubr, ggrepel, ggtree, ggtreeExtra, ggVennDiagram.

Phylogenetic information, gene counts, and gene contents are pulled from a previous project https://github.com/StochasticBiology/odna-loss -- these are the `[mt/pt]-*.csv` datafiles. 

`species-list-info.csv` is a combination of the gene count information from this pipeline and manual annotation of various ecological features of the species involved. The `-pruned` version has some awkward outliers (misannotated records etc) removed prior to analysis. The ecological features are summarised from a set of sources including The Encyclopedia of Life, Wikipedia, World Register of Marine Species, the USDA PLANTS Database, World Flora Online, and AlgaeBase.

`tree-for-traits-clean-mt.phy` is a phylogeny from the NCBI Taxonomy Common Tree tool, connecting the observed species.

`coevol-loop.R`, the main analysis code, loops through several different subsets of the complete dataset and works through the analysis: Naive and corrected linear models, mixed models, phylogenetic linear models are used to explore mt-pt gene count correlations. Specific ecological features are also explored. At the gene-specific level, clustering, Oncotrees, and minimum-leaf arborescences are used to explore mt-pt coevolution. The code produces and labels plots as it goes, outputting some to files for use in an associated manuscript.
