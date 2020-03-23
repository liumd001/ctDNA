################################################################
#function to call informative SNP and get the mixing fraction
################################################################

@control: data file of control sample;
@sample: data files of samples;
@p_thresh_for_SNP: threshold of p value to call a SNP informative
@p_thresh_for_fraction: threshold of p value to claim mixing fraction

callInformativeSNP =  function(control,
                               sample,
                               p_thresh_for_SNP = 0.01, 
                               p_thresh_for_fraction = 0.001,
                               missing = c("remove","add1"),
                               outdir = NULL){
  
  ## Function to calculate R
  calculateR = function(rf){
    tmp = data.frame(h1 = rf -0,
                     h2 = rf - 0.5,
                     h3 = rf - 1)
    return(as.vector(apply(tmp, 1, function(x) return(x[abs(x) == min(abs(x))][1]))))
  }
  
  #create a list to save file
  dta = vector(mode = "list", length = length(control) + length(sample))
  
  for (idx in 1:length(dta)){
    thisfile = c(control,sample)[idx]
    tmp = read.csv(thisfile, stringsAsFactors = FALSE)
    tmp$sample = thisfile
    
    #reshape to wide
    dta[[idx]] = reshape(tmp, v.names = "count", idvar = "dbSNP.ID", timevar = "allele", direction = "wide")
    dta[[idx]] = dta[[idx]][,c("sample","dbSNP.ID","count.R", "count.A")]
  }
  
  #combining wide data
  dta = do.call(rbind, dta)
  
  dta$sample = factor(dta$sample, levels = c(control, sample))
  
  #handle missing
  if (missing == "remove"){
    dta <- na.omit(dta)
  }else if(missing == "add1"){
    dta[is.na(dta[,4]),4] <- 1
    dta[is.na(dta[,3]),3] <- 1
  }
  
  #data transforming
  dta$pctR = dta$count.R/rowSums(dta[,c("count.R", "count.A")])
  
  #calculate minimum difference between pctR and 0,0.5,1
  tmp1 = data.frame(h1 = dta$pctR - 0,
                   h2 =  dta$pctR - 0.5,
                   h3 = dta$pctR - 1
  )
  
  dta$R = as.vector(apply(tmp1, 1, function(x) return(x[abs(x) == min(abs(x))])))
  dta$absR = abs(dta$R)
  
  #identify informative SNP
  #ctr_mean = mean(dta$R[dta$sample == control])
  ctr_sd = sd(dta$R[dta$sample == control])
  ctr_mean = 0
  
  #calculate p values 
  dta$pvalue = NA
  dta$pvalue = pnorm(dta$R, mean = ctr_mean, sd = ctr_sd)

  dta$pcat = FALSE
  dta$pcat[dta$sample %in% sample & (dta$pvalue <  p_thresh_for_SNP| dta$pvalue > (1 - p_thresh_for_SNP))] = TRUE

  #mixing fraction
  op2 = as.data.frame(matrix(NA, ncol= 3, nrow = length(sample)))
  
  for (i in 1:length(sample)){
    samp = sample[i]
    dta_s = dta[dta$sample %in% c(control, samp) & dta$dbSNP.ID %in% dta$dbSNP.ID[dta$pcat[dta$sample == samp]],]
    myfit = lm(absR ~ sample + dbSNP.ID, data = dta_s)
    op2[i,] = summary(myfit)$coefficients[2, c(1,2,4)]
  }
  op2$sample = sample
  colnames(op2) = c("Fraction","SE","p_value","Sample")
  
  if(!is.null(outdir)){
    outputfile1 = file.path(outdir,"InforamtiveSNPs.csv")
    outputfile2 = file.path(outdir,"MixingFractions.csv")
  }else{
    outputfile1 = "InforamtiveSNPs.csv"
    outputfile2 = "MixingFractions.csv"
  }
  write.csv(dta, outputfile1)
  write.csv(op2, outputfile2)
  
 return(list(output1 = dta,
             output2 = op2))
}


##Example
if(FALSE){
 dtafolder = "/home/lmd/Documents/job_application/invitae_assignment"
 myout = callInformativeSNP(control = file.path(dtafolder,"hw_sample_0.csv"),
                            sample = file.path(dtafolder,c("hw_sample_1.csv","hw_sample_2.csv")),
                            missing = "add1")
}