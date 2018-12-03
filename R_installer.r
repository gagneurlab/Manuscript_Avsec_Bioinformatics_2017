install.packages(readLines("Manuscript_Avsec_Bioinformatics_2017/r_packages.txt"))
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(readLines("Manuscript_Avsec_Bioinformatics_2017/r_bioc_packages.txt"))
