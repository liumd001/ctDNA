This project aims at 

    1) developing an algorithm to detect allele frequency changes of SNPs with SNP NGS data when sample DNA is from different sources;
    
    2) developing an algorithm to estimate the mixing fractions of different sources; 
    
    3) wraping up the algorithm into an application for use.
    
 Files:
 
    1) ctDNA_assignment.Rmd: RMD file documenting the algorithm developing process and R code. For easy read without R code, please use ctDNA_assignment.html.
    2) callinformativeSNP.R: wrap-up function, its input include: 
        control (files of control sample);
        sample: sample files; 
        outdir: the directory to output;
        p_thresh_for_SNP: threshold of p value to call informative SNP;
        missing: missingness handling: "remove" or "add1".
    3) InformativeSNPs.csv/MixingFractions.csv : two output files
    4) snp_genotyping.csv: csv file to show how mixing change the allele frequency given major and minor genotypes.
