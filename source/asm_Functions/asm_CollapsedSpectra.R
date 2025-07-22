asm_CollapsedSpectra <- function(Spectra){
     
     nSpec = length(Spectra[[1]])
     mz = NULL;
     ab = NULL;
     
     for(i in 1:nSpec){
          mz_i = as.numeric(Spectra[[1]][i][[1]][,1])
          ab_i = as.numeric(Spectra[[1]][i][[1]][,2])
          mz = c(mz,mz_i)
          ab = c(ab,ab_i)
     }
     
     # CollapsedSpec = as.data.table(cbind(mz,ab))
     # CollapsedSpec = CollapsedSpec[order(mz),]
     
     unique_mz = sort(unique(round(mz,3)))
     
     mz_f = unique_mz
     ab_f = numeric(length(unique_mz))
     for(i in 1:length(unique_mz)){
        a = which(round(mz,3) == unique_mz[i])
        ab_f[i] = sum(ab[a])
     }
     
     return(CollapsedSpec=cbind(mz=mz_f,ab=ab_f))
     
}