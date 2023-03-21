# evoDNAm_paper

This is the code associated to the following article:  
*DNA methylation patterns of transcription factor binding regions characterize their functional and evolutionary contexts*
Martina Rimoldi, Ning Wang, Jilin Zhang, Diego Villar, Duncan T. Odom, Jussi Taipale, Paul Flicek, Maša Roller
bioRxiv 2022.07.21.500978; doi: https://doi.org/10.1101/2022.07.21.500978 

There are 4 bash pipelines, they each run with their config file:
- BSvsChip.sh : performs the integration of ChIPseq and WGBS data  
Usage: BSvsChip.sh [-c config file] [-s stage]  
config file: BSvsChip.config  
  
- Annotations.sh : generate the annotation of TFBRs with methylome fragmentation information  
Usage: Annotations.sh [-c config file ]  
[-c config file]: provide config file  
[-s stage]: if specified, it must be either "CGI", "MethFrag", "cluster"    
config file: annotations.config    
  
- BPR.sh : identification of clustered methylation profiles
Usage: BPR.sh [-c config file] [-s stage]  
[-c config file]: provide config file  
[-s stage]: if specified, it must be either "input", "parameters", "clusters" or "features"  
config file: BPR.config  
  
- TFbindingDivergence.sh: performs cross species projection and annotation of TFBRs
Usage: TFbindingDivergence.sh [-c configFile] [-s stage]  
[-c config]: provide config file as input  
[-s stage]: must be one of "EPO", "intersect", "CrossSpecies", "LostPeak" or "finalTables"  
-h: print this help and exit  
config file: TFbindingDivergence.config

The figures and statistics of the article are within 4 Rmarkdown files:  
- MethAtTFBRs_Fig1-2.Rmd  
- Cluster_features_fig3.Rmd
- TFbindingDivergence_Fig4.Rmd

