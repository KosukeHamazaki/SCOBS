########################################################################################################################################################
######  Title: 3.4_SimulatedCrop_SCOBS_simulate_breeding_scheme_via_genomic_selection_with_simply_optimized_selection_index_copy_files_for_paper  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                                                          ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                                                                 ######
######  Date: 2022/05/16 (Created), 2022/05/16 (Last Updated)                                                                                     ######
########################################################################################################################################################





###### 1. Settings ######
##### 1.0. Reset workspace ######
# rm(list=ls())



##### 1.1. Setting working directory to the "SCOBS" directory #####
cropName <- "SimulatedCrop"
project <- "SCOBS"
os <- osVersion

isRproject <- function(path = getwd()) {
  files <- list.files(path)

  if (length(grep(".Rproj", files)) >= 1) {
    out <- TRUE
  } else {
    out <-  FALSE
  }
  return(out)
}

if (!isRproject()) {
  if (os == "macOS Mojave 10.14.5") {
    dirResearchBase <- "/Users/hamazaki/research/"  ### for mac OS
    dirScriptBase <- "~/research_secret/"
  } else if (os == "Ubuntu 16.04.6 LTS") {
    dirResearchBase <- "/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/"   ### for Ubunutu 1
    dirScriptBase <- "~/research_secret/"
  } else if (os == "Ubuntu 18.04.5 LTS") {
    dirResearchBase <- "/media/hamazaki/HDD2/research/"     ### for Ubuntu 2
    dirScriptBase <- "~/GitHub/research_secret/"
  } else if (os == "Ubuntu 20.04.2 LTS") {
    dirResearchBase <- "/media/hamazaki/HDD1/research/"     ### for Ubuntu 2
    dirScriptBase <- "~/GitHub/research_secret/"
  } else {
    stop("Which type of work-station do you use?")
  }

  setwd(paste0(dirResearchBase, cropName, "/Project/", project))
}


scriptID <- "3.4"





dirRes <- "results/"
dirResFp <- paste0(dirRes, scriptID, "_",
                   "Scenario_1_Trial_2_Geno_1_model_true_VerySimpleGS_results_for_paper/")
if (!dir.exists(dirResFp)) {
  dir.create(dirResFp)
}


dirGeno <- "midstream/0.1_Scenario_1/0.1_SCOBS_Trial_2/0.1_Geno_1/"
filesForTenTraits <- list.files(dirGeno)
# whichMoveFormat <- unlist(sapply(X = c("*.csv", "*.png", "*.rds"),
#                                  FUN = function(pattern) {
#                                    grep(pattern = pattern,
#                                         x = filesForTenTraits)
#                                  }))
whichMoveScriptIDForTenTraits <- grep(pattern = "3.3_*",
                                      x = filesForTenTraits)
filesForTenTraitsMove <- filesForTenTraits[whichMoveScriptIDForTenTraits]

file.copy(from = paste0(dirGeno, filesForTenTraitsMove),
          to = paste0(dirResFp, filesForTenTraitsMove))


dirModelTrue <- paste0(dirGeno, "0.1_Pheno_1/0.1_Model_with_true/")
filesForOneTrait <- list.files(dirModelTrue)
whichMoveScriptIDForOneTrait <- grep(pattern = "3.3_*",
                                     x = filesForOneTrait)
filesForOneTraitMove <- filesForOneTrait[whichMoveScriptIDForOneTrait]

file.copy(from = paste0(dirModelTrue, filesForOneTraitMove),
          to = paste0(dirResFp, filesForOneTraitMove))


dirStrategies <- filesForOneTrait[grep(pattern = "2.2_",
                                       x = filesForOneTrait)]
sapply(X = dirStrategies,
       FUN = function(dirStrategyEach) {
         dirStrategyEachFull <- paste0(dirModelTrue, dirStrategyEach, "/")
         filesForConvergence <- list.files(dirStrategyEachFull)
         whichMoveScriptIDForConvergence <- grep(pattern = "3.3_*",
                                              x = filesForConvergence)
         filesForConvergenceMove <- filesForConvergence[whichMoveScriptIDForConvergence]

         file.copy(from = paste0(dirStrategyEachFull, filesForConvergenceMove),
                   to = paste0(dirResFp, filesForConvergenceMove))
       })
