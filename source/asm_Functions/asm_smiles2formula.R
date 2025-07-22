library(stringr)
asm_struc2formula = function(smiles){
     
     smiles = toupper(smiles)
     
     #Cl
     #Br
     #
     str_count
     
     a = grepl("CL",smiles)
     tempData = strsplit(smiles,"")[[1]]
     
     a = which(is.na(suppressWarnings(as.numeric(tempData))))
     if(length(a)>0)
     atoms = tempData[]
     blanks = which(tempData=="")
     tempData = tempData[-blanks]
     natoms = as.numeric(tempData[1])
     nconns = as.numeric(tempData[2])
     
     atomblock = character(natoms)
     connblock = array(0,dim=c(nconns,3))
     atomS = 5
     atomE = atomS + natoms - 1
     connS = atomE + 1
     connE = connS + nconns - 1
     
     DeuteratedInfoLine = grep("ISO",sdf_struc)
     if(length(DeuteratedInfoLine)>0){
       DeuteratedInfo = strsplit(sdf_struc[DeuteratedInfoLine]," ")[[1]]
       spaces = which(DeuteratedInfo=="")
       DeuteratedInfo = DeuteratedInfo[-spaces]
       DeuteratedSet = DeuteratedInfo[seq(4,length(DeuteratedInfo),2)]
     } else {
       DeuteratedSet = NULL
     }
     
     
     p = seq(atomS,atomE)
     for(i in 1:natoms){
          t = strsplit(sdf_struc[p[i]], " ")[[1]]
          spaces = which(t=="")
          t = t[-spaces]
          if(i %in% DeuteratedSet){
            if(t[4]=="H"){
              atomblock[i] = "D"
            }
          } else {
            atomblock[i] = t[4]  
          }
     }
     
     p = seq(connS,connE)
     for(i in 1:nconns){
          t = strsplit(sdf_struc[p[i]], " ")[[1]]
          spaces = which(t=="")
          t = t[-spaces]
          connblock[i,] = c(t[1],t[2], t[3]) 
     }
     

     FullMol = table(atomblock)
     FullMol_formula = NULL
     for(j in 1:dim(FullMol)){
          FullMol_formula = paste0(FullMol_formula,names(FullMol)[j],FullMol[j])
     }
     
     return(FullMol_formula)
}