# odna-coevolution
Coevolution of mtDNA and ptDNA genome content

Requires R with libraries phytools, ape, phangorn, phylolm, nlme, lme4, pheatmap, Oncotree, igraph, ggplot2, ggpubr, ggrepel, ggtree, ggtreeExtra, ggVennDiagram.

Phylogenetic information, gene counts, and gene contents are pulled from a previous project https://github.com/StochasticBiology/odna-loss -- these are the `.csv` datafiles. `coevol-loop.R`, the main analysis code, loops through several different subsets of the complete dataset and works through the analysis: Naive and corrected linear models, mixed models, phylogenetic linear models are used to explore mt-pt gene count correlations. Specific ecological features are also explored. At the gene-specific level, clustering, Oncotrees, and minimum-leaf arborescences are used to explore mt-pt coevolution. The code produces and labels plots as it goes, outputting some to files for use in an associated manuscript.
