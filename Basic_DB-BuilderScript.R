## CRAFTS LAB Database Builder (DB) Script - Basic DB
##
## Generates R data table (.RDS) for use with NIST/NIJ DART-MS DIT (v3.22).
##
## Developed by: A.S.M; arunmoorthy@trentu.ca
##
## Adapted from code originally developed by:
##.              A.S.M; arunmoorthy@trentu.ca
##               E.S; edward.sisco@nist.gov
##               S.S.T; stephen.tennyson@nist.gov
##               M.G.A; meghan.appley@nist.gov
##
## Revised: July 20, 2025
## =============================================================================

# Setup
rm(list=ls())

header = "CRAFTS LAB Basic Databse Builder\nRevised July 20th, 2025.\n\nA basic builder with no quality checks."
cat(header)
cat("\n\n")

parent_path = getwd()
parent_path = paste0(parent_path,"/source")
child_path = paste0(parent_path,"/asm_Preliminaries.R")
source(child_path)

user_OS = .Platform$OS.type   # for OS specific system commands
DEThreshold = 0.30;           # Dimer Error Threshold
# NoiseThreshold = 0.45;        # Unexplained Peak Intensity Threshold
# IsotopeRatioThreshold = 0.90; # similarity between  measured
#                               # and calculated molecular ion envelope


## MAIN Loop
operation = TRUE
while (operation){

  ## Step 0: USER Data
  potential_master_files = list.files(".",pattern=".xlsx")
  if(length(potential_master_files)<1){
    operation = asm_ErrorFlagFatal("Step 0. No potential master_files in current directory.")
    break
  } else {

    isError = TRUE
    while (isError){
      cat(paste0("There are ", length(potential_master_files), " potential master files currently in the directory.\n\n"))
      for(i in 1:length(potential_master_files)){
        cat(paste0(i,": ",potential_master_files[i],"\n"))
      }
      cat("\n")

      a <- readline(prompt="Choose Master File (integer value) for creating database: ")
      if (a %in% as.character(seq(1,length(potential_master_files)))){
        master_file = potential_master_files[as.numeric(a)]
        isError = FALSE
        cat("\n")
      } else {
        cat("\n")
        cat("INVALID INPUT. Try again.\n")
      }
    }
  }

  potential_main_folders = list.dirs(".",recursive=FALSE)
  a = which(potential_main_folders == "./source")
  if(length(a)>0) potential_main_folders = potential_main_folders[-a]
  a = which(potential_main_folders == "./.git")
  if(length(a)>0) potential_main_folders = potential_main_folders[-a]
  
  if(length(potential_main_folders)<1){
    operation = asm_ErrorFlagFatal("Step 0. No potential main_folders in current directory.")
    break
  } else {

    isError = TRUE
    while (isError){
      cat(paste0("There are ", length(potential_main_folders), " potential main folders currently in the directory.\n\n"))
      for(i in 1:length(potential_main_folders)){
        folder_name = strsplit(potential_main_folders[i],"/")[[1]][2]
        cat(paste0(i,": ",folder_name,"\n"))
      }
      cat("\n")

      a <- readline(prompt="Choose Main Folder (integer value) for creating database: ")
      if (a %in% as.character(seq(1,length(potential_main_folders)))){
        main_folder = strsplit(potential_main_folders[as.integer(a)],"/")[[1]][2]
        isError = FALSE
        cat("\n")
      } else {
        cat("\n")
        cat("INVALID INPUT. Try again.\n")
      }
    }
  }

  potential_RequireFolderStructure = list(c("LOW","MID","HIGH"),c("1","2","3"), c("+30 V","+60 V","+90 V"))
  isError = TRUE
  while (isError){
    cat(paste0("Subfolder structures currently allowed by the DB builder:\n\n"))
    for(i in 1:length(potential_RequireFolderStructure)){
      cat(paste0(i,": "))
      for(j in 1:length(potential_RequireFolderStructure[[i]])){
        cat(paste0(potential_RequireFolderStructure[[i]][j]," "))
      }
      cat("\n")
    }
    cat("\n")

    a <- readline(prompt="Choose subfolder structure (integer value) for creating database: ")
    if (a %in% as.character(seq(1,length(potential_RequireFolderStructure)))){
      RequireFolderStructure = potential_RequireFolderStructure[[as.numeric(a)]]
      isError = FALSE
      cat("\n")
    } else {
      cat("\n")
      cat("INVALID INPUT. Try again.\n")
    }
  }

  mz_res = 0.005 # default value
  isError = TRUE
  while (isError){

    a <- readline(prompt="Enter m/z tolerance in Daltons (e.g. 0.005): ")
    b <- suppressWarnings(as.numeric(a))
    if (!is.na(b)){
      if(b>0){
        mz_res = b
        isError = FALSE
      } else {
        cat("INVALID INPUT. The m/z tolerance must be greater than 0 Daltons. Try again.\n")
      }
    } else {
      cat("INVALID INPUT. Try again.\n")
    }
  }

  # MAY Revisit for future iterations
  # potential_build_style = c("Traditional","Collapsed")
  # isError = TRUE
  # while (isError){
  #   cat(paste0("Build a traditional (low/med/high) or collapsed spectral library?\n\n"))
  #   for(i in 1:length(potential_build_style)){
  #     cat(paste0(i,": ",potential_build_style[i],"\n"))
  #   }
  #   cat("\n")
  # 
  #   a <- readline(prompt="Choose build type (integer value) for creating database: ")
  #   if (a %in% as.character(seq(1,length(potential_build_style)))){
  #     build_style = potential_build_style[as.numeric(a)]
  #     isError = FALSE
  #     cat("\n")
  #   } else {
  #     cat("\n")
  #     cat("INVALID INPUT. Try again.\n")
  #   }
  # }
  build_style = "Traditional"

  output_name = "Output_DARTMS_Database"
  isError = TRUE
  while (isError){

    a <- readline(prompt="Enter name for generated database files (e.g. output_database): ")
    output_name = a
    isError = FALSE
  }

## MAIN CODE
  # Preliminary Data Validation (make sure the spectra are txt files)
  asm_jsp2txt(main_folder)

  ## Definition of monoisotpic mass for atoms of interest
  child_path = paste0(parent_path,"/asm_MIM_Definitions.R")
  source(child_path)

  cat("Step 1. Initiating database metadata by reading master_file.\n")
  LibMaster = as.data.table(read_excel(master_file))
  nCompounds = dim(LibMaster)[1]

  cat('Step 2. Confirming correct number of spectra exist in each "energy" subdirectory of the main_folder.\n\n')
  EnergyFolders = list.dirs(main_folder,recursive = FALSE)
  nEnergyFolders = length(EnergyFolders)
  ReplicateList = NULL;

    if(nEnergyFolders<1) {
    operation = asm_ErrorFlagFatal("Step 2. Improper main_folder format. No subdirectories with spectra.")
    break
    }
    d = 0;
    for(i in 1:nEnergyFolders){
    InFolderSpectra = list.files(EnergyFolders[i],pattern = ".txt")
    nInFolderSpectra = length(InFolderSpectra)

          a = unlist(strsplit(InFolderSpectra,".txt"))
          a2 = grep("-",a)
          if(length(a2)>0){
            ReplicateList = c(ReplicateList,strsplit(a[a2],"-")[[1]][1])
            a = a[-a2]
          }

          b = a %in% LibMaster[,Code]
          c = which(b==FALSE);
          if(length(c)>0){
          message = paste("The following spectra in ",
                          EnergyFolders[i],
                          " do not have metadata in the  master file:\n",
                          sep="");
          d = d + 1;
          for(j in 1:length(c)){
            message = paste(message,a[c[j]],"\n",sep="")
          }
          cat(message);
          cat("\n")
          }

          b = LibMaster[,Code] %in% a;
          c = which(b==FALSE);
          if(length(c)>0){
            message = paste("The following codes with metadata in the master file",
                            " do not have spectra in ",
                          EnergyFolders[i],
                          ":\n",
                          sep="");
            d = d + 1;
            for(j in 1:length(c)){
             message = paste(message,LibMaster[c[j],Code],"\n",sep="")
          }
          cat(message);
          cat("\n")
          }

    }

    if (d>0){
      operation = asm_ErrorFlagFatal("Step 2: Folder and metadata errors.")
      break
    }
    

  cat('\nStep 3. Collecting spectra from folders.\n')

  eFixed = character(nEnergyFolders)
    for(i in 1:nEnergyFolders){
      a = strsplit(EnergyFolders[i],"/")[[1]]
      eFixed[i] = a[2]
    }

    if(RequireFolderStructure[1] %in% eFixed == FALSE){
      operation = asm_ErrorFlagFatal("Step C2: No LOW fragmentation spectra folder found (required).")
      break
    }

    if(sum(eFixed %in% RequireFolderStructure) != nEnergyFolders){
      operation = asm_ErrorFlagFatal("Step C2: Spectra subfolder labeling inconsistent with specifications.")
      break
    }

  eFixed = RequireFolderStructure[RequireFolderStructure %in% eFixed]

  peaks = character(nEnergyFolders)       # individual peak list
  numPeaks = numeric(nEnergyFolders)      # number of peaks per list
  Energies = character(nCompounds)        # likely fixed at nEnergyFolders
  PeakLists = character(nCompounds)       # combined peak lists for a given entry
  NumSpectra = numeric(1);                # likely fixed at nEnergyFolders
  NumPeaksList = character(nCompounds)    # likely fixed at nEnergyFolders
  DimerProb = numeric(nCompounds)         # check if dimer
  DimerErrorProb = numeric(nCompounds)    # abundance ratio of dimer : base peak
  MassCaliError = numeric(nCompounds)     # mass error
  PotentialBPs = numeric(nCompounds)      # number of "potential" base peaks
  BP = numeric(nCompounds)                # m/z value of base peak in +30 V spectrum
  PotentialErrors = numeric(nCompounds)   # potential error spectra (mass calibration)
  formulasubscript = character(nCompounds)# formula with subscript # added by Stephen
  polarity = character(nCompounds)        # the library's polarity (positive or negative)
  sourceGas = character(nCompounds)       # the library's source gas (He, N, etc.)
  IsotopeRatioSim = numeric(nCompounds)   # the similarity between measured and computed protonated molecule envelope
  CollapsedSpectra = character(nCompounds) # a placeholder for collapsed mass spectra

  FragmentationMetrics = numeric(nCompounds) # a set of internal metrics to measure fragmentation consistency
  PotentialErrorsFM1 = numeric(nCompounds) # potential errors based on fragmentation inconsistency

    pb <- txtProgressBar(min=0,max=(nInFolderSpectra*nEnergyFolders),initial=0,char="=",style=3)
    counter = 0
    for (i in 1:nInFolderSpectra) {
    h = strsplit(InFolderSpectra[i],"\\.")[[1]][1] # split at the period (get just the code)
    k = which(LibMaster[,Code]==h)

    preMZ = LibMaster[k,3][[1]];
    for(j in 1:nEnergyFolders) {

      txtDirectory = paste(main_folder,"/",eFixed[j],"/",InFolderSpectra[i],sep="")
      # cat(paste(file.info(txtDirectory)$size, "\n",sep="")); # diagnostic for over-sized spectra files

      data = readLines(txtDirectory)
      b = asm_ListCreator_v1.02(data,preMZ) ## updated to acquire base peak information
      peaks[j] = b[1];
      if(j==1){
        MassCaliError[k] = as.numeric(b[[2]])
        DimerProb[k] = as.numeric(b[[3]])
        DimerErrorProb[k] = as.numeric(b[[5]])
        PotentialBPs[k] = as.numeric(b[[4]])
        BP[k] = as.numeric(b[[6]])
        if(abs(as.numeric(b[[2]]))>mz_res){
          PotentialErrors[k]=1
        }
      }
      .npeaks = length(data)
      numPeaks[j] = .npeaks;

      counter = counter + 1;
      setTxtProgressBar(pb,counter)

    }

    Energies[k] = list(eFixed)
    PeakLists[k] = list(peaks)
    NumPeaksList[k] = list(numPeaks)
    NumSpectra[k] = length(EnergyFolders);
  }
  
  cat('Step 4. Creating library for the DIT.\n')
  
  Structure_gen = character(nCompounds)
  for(i in 1:nCompounds){
    Structure_gen[i] = list(c("",                                                                     
                              "Arun's Random Structure Block",                                               
                              "",                                                                     
                              " 4 1  0  0  0  0  0  0  0  0999 V2000",                              
                              "    0  1    0.0000 No   0  0  0  0  0  0  0  0  0  0  0  0",
                              "    0  -1    0.0000 Structure   0  0  0  0  0  0  0  0  0  0  0  0",
                              "    0   2    0.0000 .   0  0  0  0  0  0  0  0  0  0  0  0",
                              "    0  -2    0.0000 .   0  0  0  0  0  0  0  0  0  0  0  0",
                              "  1  2  1  0  0  0  0",                                                
                              "M  END"))
  }
  
  Energies_306090 = rep(list(c("+30 V","+60 V","+90 V")),nCompounds)
  Library_RDT = as.data.table(cbind(LibMaster[,1],
                                    LibMaster[,2],
                                    'Cas #' = rep("NA",nCompounds), 
                                    'Synonymns' = rep("NA",nCompounds),
                                    'IUPAC/ Formal Name' = rep("NA",nCompounds),
                                    'Formula' = rep("NA",nCompounds),
                                    'AccurateMolecularMass' = c(LibMaster[,3][[1]]),
                                    'Class' = rep("NA",nCompounds), 
                                    'Canonical SMILES' = rep("NA",nCompounds), 
                                    'InChi Code' = rep("NA",nCompounds),
                                    'InChIKey' = rep("NA",nCompounds),
                                    'MW_gen' = c(LibMaster[,3][[1]]),
                                    'AccurateMass_gen' = c(LibMaster[,3][[1]]),
                                    'PrecursorMZ_gen' = c(LibMaster[,3][[1]]),
                                    'Energies' = Energies_306090,
                                    NumPeaksList,
                                    NumSpectra,
                                    PeakLists,
                                    'SMILES' = rep("NA",nCompounds),
                                    'InChIKey_gen' = rep("NA",nCompounds),
                                    MassCaliError,
                                    DimerProb,
                                    DimerErrorProb,
                                    BP,
                                    'theoBP' = rep(NaN,nCompounds),
                                    'theoBP_MolForm' = rep("NA",nCompounds),
                                    PotentialBPs,
                                    PotentialErrors,
                                    FragmentationMetrics,
                                    PotentialErrorsFM1,
                                    Structure_gen,
                                    'RefinedAnnotations' = rep("NA",nCompounds),
                                    'NoiseMetric' = rep("NA",nCompounds),
                                    'formulasubscript' = rep("NA",nCompounds),
                                    'polarity' = rep("NA",nCompounds),
                                    'sourceGas' = rep("NA",nCompounds),
                                    'IsotopeRatioSim' = rep(NaN,nCompounds),
                                    'pmMajIsoMZ30V' = rep(NaN,nCompounds),
                                    'pmMajIsoAb30V' = rep(NaN,nCompounds),
                                    'bpMajIsoMZ30V' = rep(NaN,nCompounds),
                                    'bpMajIsoAb30V' = rep(NaN,nCompounds),
                                    'mfMajIsoMZ30V' = rep(NaN,nCompounds),
                                    'mfMajIsoAb30V' = rep(NaN,nCompounds),
                                    'mfMZ30V' = rep(NaN,nCompounds)
                                    ))

  RDSfilename = paste0(output_name,".RDS")
  saveRDS(Library_RDT,RDSfilename)

  ## End code (check if another library is to be built?)
  potential_continue = c("Yes","yes","Y","y")
  potential_exit = c("No","no","N","n")
  isError = TRUE
  cat("\n")
  while (isError){
    a <- readline(prompt="Do you want to create another database? (yes/no) ")
    if (a %in% potential_continue){
      isError = FALSE
    } else if (a %in% potential_exit){
      cat("\nExiting NIST DART-MS Database Builder program.\n\n")
      operation = FALSE
      isError =  FALSE
    } else {
      cat("INVALID INPUT. Do you want to create another database? (yes/no) ")
    }
  }




}
