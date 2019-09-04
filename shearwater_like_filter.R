#-------------------------------------------------
# Shearwater-like filter for WGS (post-CaVEMan)
# Tim Coorens - November 2018
#-------------------------------------------------
options(stringsAsFactors=F)

#-------------------------------------------------
# Libraries
#-------------------------------------------------

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

#-------------------------------------------------
# Functions
#-------------------------------------------------

estimateRho_gridml = function(x, mu) {
  # Estimate rho by MLE grid approach
  rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
  mm = x[,2]
  #cov = c(x[,3:4])+c(x[,1:2])
  cov = c(x[,1])
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
  rhovec[ll==max(ll)][1]
}

shearwater_probability=function(patient, save=NULL, path_prefix='', rho=10^-3){
  #Function to calculcate probability of presence of mutation based on Shearwarer
  #'patient' is the name of the patient-specific subfolder
  #'path_prefix' is any prefix to the path necessary to find that subfolder
  #'rho' is the constant for the overdispersion parameter. If rho=NULL, calculate
  # it from the data (much slower)
  #'save' is path for output. If NULL, returns matrix
  # output of this function will be a matrix (muts by samples) of p values
  
  
  #A file of mutations in patient subdirectory (format: Chr_Ref_Pos_Alt)
  Muts_patient = read.table(paste0(path_prefix,patient,"/All_mutations_filtered.txt"))[,1] 
  #A file of with the sample names belonging to this patient
  case_samples=read.table(paste0(path_prefix,patient,"/samples.txt"))[,1]
  
  #Select all normal samples not belonging to this patient
  normal_panel = normal_samples[normal_samples%in%samples&!normal_samples%in%case_samples]
  norm_all_counts = all_counts[normal_panel,,]
  
  #Set up pval matrix
  pval_mat = matrix(1,ncol=length(case_samples),nrow=length(Muts_patient))
  rownames(pval_mat)=Muts_patient
  colnames(pval_mat)=case_samples
  
  coords_proj = substr(Muts_patient,1,nchar(Muts_patient)-4)
  Alt=substr(Muts_patient,nchar(Muts_patient),nchar(Muts_patient))
  Ref=substr(Muts_patient,nchar(Muts_patient)-2,nchar(Muts_patient)-2)
  
  for (s in case_samples){
    rho_est=rep(NA,length(coords_proj))
    test_counts = all_counts[s,coords_proj,]
    for (k in 1:length(coords_proj)) {
      n = sum(test_counts[coords_proj[k],])
      x = test_counts[coords_proj[k],Alt[k]]
      
      N_indiv = rowSums(norm_all_counts[,coords_proj[k],])
      X_indiv = norm_all_counts[,coords_proj[k],c("A","C","G","T")!=Ref[k]]
      pseudo = .Machine$double.eps    
      N=sum(N_indiv)
      X=sum(X_indiv)
      
      mu = max(X,pseudo)/max(N,pseudo)
      counts = cbind(N,X)
      if(is.null(rho)) rho = estimateRho_gridml(counts,mu)
      rdisp = (1 - rho)/rho
      
      prob0 = (X + x)/(N + n); prob0[prob0==0] = pseudo
      prob1s = x/(n+pseudo); prob1s[prob1s==0] = pseudo
      prob1c = X/(N+pseudo); prob1c[prob1c==0] = pseudo
      
      prob1s = pmax(prob1s,prob1c) # Min error rate is that of the population (one-sided test)
      nu0 = prob0 * rdisp; nu1s = prob1s * rdisp; nu1c = prob1c * rdisp; 
      
      # Likelihood-Ratio Tests
      LL = logbb(x, n, nu0, rdisp) + logbb(X, N, nu0, rdisp) - logbb(x, n, nu1s, rdisp) - logbb(X, N, nu1c, rdisp)
      pvals = pchisq(-2*LL, df=1, lower.tail=F)/2 # We divide by 2 as we are performing a 1-sided test
      # Saving the result
      pval_mat[k,s] = pvals
    } 
  }
  if(is.null(save)){
    return(pval_mat)
  }else{
    write.table(pval_mat,save)
  }
}



#-------------------------------------------------
# Input
#-------------------------------------------------

# Vector of all sample names
samples=read.table("samples_all_final.txt")[,1] 

# Vector of samples excluding tumours/aneuploid samples
normal_samples=read.table("normal_samples_all.txt")[,1] 

# "Bed" file of all mutations to be considered (across all patients)
# Format: Chr Ref Pos Alt
muts=read.table("AlleleCounts_final_June/All_mutations_filtered.bed") 
coords=paste(muts$V1,muts$V2,sep="_")

# Read in data from AlleleCounter
all_counts = array(0,dim=c(length(samples),length(coords),4),
                   dimnames=list(samples,coords,c("A","C","G","T")))
print(length(samples))
for (k in 1:length(samples)){
  #Read in allele counts per sample
  if(file.exists(paste0("AlleleCounts/",samples[k],".allelecounts.txt"))){
    data=read.table(paste0("AlleleCounts/",samples[k],".allelecounts.txt"),comment.char = '',header=T)
    #to be safe, avoid dulicated mutations in bed file:
    muts_data=paste(data$X.CHR,data$POS,sep="_")
    data=data[!duplicated(muts_data),]
    muts_data=muts_data[!duplicated(muts_data)]
    rownames(data)=muts_data
    all_counts[k,,]=as.matrix(data[coords,3:6])
  }
  print(k)
}

# Vector of all different projects/patients
patients=read.table("patients_all.txt")[,1]
for (patient in patients){
  shearwater_probability(patient=patient,save=paste0(patient,"/shearwater_pval_mat.txt"))
}
#If return mat:
pval_mat=shearwater_probability(patient=patient)
#Multiple testing correction
qval_mat=apply(pval_mat,2,function(x) p.adjust(x,method="BH",n = length(pval_mat)))
