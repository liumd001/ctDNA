This project aims at 
    1) develop algorithm to detect allele frequency changes of SNPs with SNP NGS data when a sample is mixed with another sample; 
     2) develop algorithm to estimate the mixing fractions; 3) wrap up the algorithm into an application for use.

Invitae_assignment.Rmd: RMD file documenting the algorithm developing process and R code. For easy read without R code, please use Invitae_assignment.html.

callinformativeSNP.R: wrap-up function, its input include: 
      control (files of control sample);
      sample: sample files; 
      outdir: the directory to output;
      p_thresh_for_SNP: threshold of p value to call informative SNP;
      missing: missingness handling: "remove" or "add1".
      
InformativeSNPs.csv/MixingFractions.csv : two output files

snp_genotyping.csv: csv file to show how mixing change the allele frequency given major and minor genotypes.
