########################################################################################################
######  Title: 2.1.1_SimulatedCrop_SCOBS_simulate_breeding_scheme_via_genomic_selection_until_G6  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                          ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                 ######
######  Date: 2020/01/22 (Created), 2021/01/22 (Last Updated)                                     ######
########################################################################################################





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
  } else {
    stop("Which type of work-station do you use?")
  }

  setwd(paste0(dirResearchBase, cropName, "/Project/", project))
}

scriptIDData <- "0.1"
scriptID <- "2.1"
setwd("/home/hamazaki/research/SimulatedCrop/Project/SCOBS")


##### 1.2. Setting some parameters #####
dirMidSCOBSBase <- "midstream/"


#### 1.2.1. Setting some parameters related to names of R6 class ####
scenarioNo <- 1
simName <- paste0("Scenario_", scenarioNo)
breederName <- "K. Hamazaki"

dirMidSCOBSScenario <- paste0(dirMidSCOBSBase, scriptIDData, "_", simName, "/")
if (!dir.exists(dirMidSCOBSScenario)) {
  dir.create(dirMidSCOBSScenario)
}


#### 1.2.2. Setting some parameters related to simulation ####
trialNo <- 1
trialName <- paste0("Trial_", trialNo)

dirMidSCOBSTrial <- paste0(dirMidSCOBSScenario, scriptIDData,
                           "_", project, "_", trialName, "/")
if (!dir.exists(dirMidSCOBSTrial)) {
  dir.create(dirMidSCOBSTrial)
}

nSimGeno <- 1
nSimPheno <- 10

nIterSimulation <- 1000
nGenerationProceed <- 4
nRefreshMemoryEvery <- 2
updateBreederInfo <- rep(TRUE, nGenerationProceed)
phenotypingInds <- rep(FALSE, nGenerationProceed)
nRepForPhenoInit <- 3
nRepForPheno <- rep(1, nGenerationProceed)
updateModels <- rep(FALSE, nGenerationProceed)

simBsNameBase <- "VerySimpleGS-Until-G6-for-Trait_1"
strategyNames <- c(
  # "selectBV", "selectWBV",
  # "selectOHV", "selectOPV",
  "selectEMBV"
)

nTrainingPopInits <- list("true", 2000,  1500,
                          1000, 750, 500, 250)

lociEffMethods <- c("true", "estimated")
trainingPopType <- "latest"
trainingPopInit <- 1
trainingIndNamesInit <- NULL


#### 1.2.3. Setting some parameters related to estimation of marker effects ####
methodMLRInit <- "BayesB"
multiTraitInit <- FALSE
methodMLR <- "LASSO"
multiTrait <- FALSE

alpha <- 0.5
nIter <- 20000
burnIn <- 5000
thin <- 5
bayesian <- TRUE


#### 1.2.4. Setting some parameters related to selection of parent candidates ####
nSelectionWaysVec <- rep(1, nGenerationProceed)
traitNoSelList <- rep(list(list(1)), nGenerationProceed)

blockSplitMethod <- "minimumSegmentLength"
nMrkInBlock <- 10
minimumSegmentLength <- 20
nSelInitOPVList <- rep(list(50), nGenerationProceed)
nIterOPV <- 5e04
nProgeniesEMBVVec <- rep(40, nGenerationProceed)
nIterEMBV <- 3
nCoresEMBV <- 1

clusteringForSelList <- rep(list(FALSE), nGenerationProceed)
nClusterList <- rep(list(1), nGenerationProceed)
nTopClusterList <- rep(list(1), nGenerationProceed)
nTopEachList <- rep(list(15), nGenerationProceed)
nSelList <- rep(list(15), nGenerationProceed)

multiTraitsEvalMethodList <- rep(list("sum"), nGenerationProceed)
hSelList <- rep(list(list(1)), nGenerationProceed)


#### 1.2.5. Setting some parameters related to mating and resource allocation ####
matingMethodVec <- rep("diallelWithSelfing", nGenerationProceed)
allocateMethodVec <- rep("equalAllocation", nGenerationProceed)
weightedAllocationMethodList <- rep(list("selectBV"), nGenerationProceed)
traitNoRAList <- rep(list(1), nGenerationProceed)
hList <- rep(list(1), nGenerationProceed)
includeGVPVec <- rep(FALSE, nGenerationProceed)
nNextPopVec <- rep(250, nGenerationProceed)


#### 1.2.6. Setting some parameters related to mating and resource allocation ####
nCores <-  60
nCoresWhenOPV <-  60
nCoresWhenEMBV <- 50
nCoresSummary <- 10


nameMethod <- "pairBase"
overWriteRes <- FALSE
showProgress <- TRUE
returnMethod <- "summary"
evaluateGVMethod <- "true"
nTopEval <- 5
traitNoEval <- 1
hEval <- 1
verbose <- TRUE

saveAllResAtBase <- "all_results"
summaryAllResAtBase <- "all_results"


#### 1.2.7. Save parameters ####
fileParamsSCOBS <- paste0(dirMidSCOBSTrial, scriptID,
                          "_", project, "_", trialName,
                          "_", simBsNameBase,
                          "_all_parameters.RData")
save.image(fileParamsSCOBS)



##### 1.3. Import packages #####
require(myBreedSimulatR)
require(data.table)
require(MASS)
require(rrBLUP)
require(BGLR)
require(RAINBOWR)
require(ggplot2)
require(plotly)





##### 1.4. Project options #####
options(stringAsFactors = FALSE)


# simGenoNo <- 1
# simPhenoNo <- 1
# strategyName <- strategyNames[1]
# lociEffMethod <- lociEffMethods[1]
# trainingPopInitNo <- 1
for (simGenoNo in 1:nSimGeno) {
  genoName <- paste0("Geno_", simGenoNo)

  dirMidSCOBSGeno <- paste0(dirMidSCOBSTrial, scriptIDData,
                            "_", genoName, "/")


  for (simPhenoNo in 1:nSimPheno) {
    phenoName <- paste0("Pheno_", simPhenoNo)

    dirMidSCOBSPheno <- paste0(dirMidSCOBSGeno, scriptIDData,
                               "_", phenoName, "/")
    for (trainingPopInitNo in 1:length(nTrainingPopInits)) {
      nTrainingPopInit <- nTrainingPopInits[[trainingPopInitNo]]
      simBsNameTrainingPop <- paste0(simBsNameBase, "-", nTrainingPopInit)

      dirMidSCOBSTrainingPop <- paste0(dirMidSCOBSPheno, scriptIDData,
                                       "_Model_with_", nTrainingPopInit, "/")
      if (!dir.exists(dirMidSCOBSTrainingPop)) {
        dir.create(dirMidSCOBSTrainingPop)
      }


      ###### 2. Load bsInfoInit & breederInfoInit & estimate marker effects ######
      if (nTrainingPopInit == "true") {
        dirMidSCOBSTrainingPopFalse <- paste0(dirMidSCOBSPheno, scriptIDData,
                                              "_Model_with_", nTrainingPopInits[[2]], "/")
        fileNameBsInfoInitRmInit <- paste0(dirMidSCOBSTrainingPopFalse, scriptIDData,
                                           "_bsInfoInit_remove_initial_population.rds")
        fileNameBreederInfoInitWithMrkEffRmInit <- paste0(dirMidSCOBSTrainingPopFalse, scriptIDData,
                                                          "_breederInfoInit_with_estimated_marker_effects_",
                                                          methodMLRInit, "_", nTrainingPopInits[[2]],
                                                          "_individuals_remove_initial_population.rds")

        lociEffMethod <- lociEffMethods[1]
      } else {
        fileNameBsInfoInitRmInit <- paste0(dirMidSCOBSTrainingPop, scriptIDData,
                                           "_bsInfoInit_remove_initial_population.rds")
        fileNameBreederInfoInitWithMrkEffRmInit <- paste0(dirMidSCOBSTrainingPop, scriptIDData,
                                                          "_breederInfoInit_with_estimated_marker_effects_",
                                                          methodMLRInit, "_", nTrainingPopInit,
                                                          "_individuals_remove_initial_population.rds")

        lociEffMethod <- lociEffMethods[2]
      }
      myBsInfoInit <- readRDS(file = fileNameBsInfoInitRmInit)

      # if (file.exists(fileNameBreederInfoInitWithMrkEff)) {
      myBreederInfoInit <- readRDS(file = fileNameBreederInfoInitWithMrkEffRmInit)
      # } else {
      #   fileNameBreederInfoInit <- paste0(dirMidSCOBSPheno, scriptIDData,
      #                                     "_breederInfoInit.rds")
      #   myBreederInfoInit <- readRDS(file = fileNameBreederInfoInit)
      #
      #   myBreederInfoInit$estimateMrkEff(trainingPop = NULL,
      #                                    methodMLR = methodMLRInit,
      #                                    multiTrait = multiTraitInit,
      #                                    alpha = alpha,
      #                                    nIter = nIter,
      #                                    burnIn = burnIn,
      #                                    thin = thin,
      #                                    bayesian = bayesian)
      #
      #
      #   saveRDS(object = myBreederInfoInit, file = fileNameBreederInfoInitWithMrkEff)
      # }


      # mySimBsBVList <- list()
      for (strategyName in strategyNames) {
        print(paste0(genoName, "  ", phenoName, "  ",
                     nTrainingPopInit, "  ",
                     strategyName, "  ", lociEffMethod))
        st <- Sys.time()
        simBsName <- paste0(simBsNameTrainingPop, "-",
                            strategyName, "-", lociEffMethod)
        selectionMethodList <- rep(list(strategyName), nGenerationProceed)




        dirMidSCOBSBV <- paste0(dirMidSCOBSTrainingPop, scriptID, "_", simBsName, "/")
        if (!dir.exists(dirMidSCOBSBV)) {
          dir.create(dirMidSCOBSBV)
        }

        # if ((simGenoNo == 1) & (simPhenoNo == 1)) {
        #   saveAllResAt <- paste0(dirMidSCOBSBV, scriptID, "_", saveAllResAtBase)
        # } else {
        #   saveAllResAt <- NULL
        # }
        saveAllResAt <- NULL
        # summaryAllResAt <- paste0(dirMidSCOBSBV, scriptID, "_", summaryAllResAtBase)
        summaryAllResAt <- NULL


        fileNameSimBsBVNow <- paste0(dirMidSCOBSBV, scriptID, "_", simBsName, "_simBs.rds")

        if (file.exists(fileNameSimBsBVNow)) {
          mySimBsBVNow <- readRDS(file = fileNameSimBsBVNow)

          # mySimBsBVList[[simBsName]] <- mySimBsBVNow

        } else {

          ###### 3. Simulate breeding scheme using simple GS with each BV ######
          ##### 3.1. Set simBs class object #####
          mySimBsBVNow <- myBreedSimulatR::simBs$new(simBsName = simBsName,
                                                     bsInfoInit = myBsInfoInit,
                                                     breederInfoInit = myBreederInfoInit,
                                                     lociEffMethod = lociEffMethod,
                                                     trainingPopType = trainingPopType,
                                                     trainingPopInit = trainingPopInit,
                                                     trainingIndNamesInit = trainingIndNamesInit,
                                                     methodMLRInit = methodMLRInit,
                                                     multiTraitInit = multiTraitInit,
                                                     nIterSimulation = nIterSimulation,
                                                     nGenerationProceed = nGenerationProceed,
                                                     nRefreshMemoryEvery = nRefreshMemoryEvery,
                                                     updateBreederInfo = updateBreederInfo,
                                                     phenotypingInds = phenotypingInds,
                                                     nRepForPhenoInit = nRepForPhenoInit,
                                                     nRepForPheno = nRepForPheno,
                                                     updateModels = updateModels,
                                                     methodMLR = methodMLR,
                                                     multiTrait = multiTrait,
                                                     nSelectionWaysVec = nSelectionWaysVec,
                                                     selectionMethodList = selectionMethodList,
                                                     traitNoSelList = traitNoSelList,
                                                     blockSplitMethod = blockSplitMethod,
                                                     nMrkInBlock = nMrkInBlock,
                                                     minimumSegmentLength = minimumSegmentLength,
                                                     nSelInitOPVList = nSelInitOPVList,
                                                     nIterOPV = nIterOPV,
                                                     nProgeniesEMBVVec = nProgeniesEMBVVec,
                                                     nIterEMBV = nIterEMBV,
                                                     nCoresEMBV = nCoresEMBV,
                                                     clusteringForSelList = clusteringForSelList,
                                                     nClusterList = nClusterList,
                                                     nTopClusterList = nTopClusterList,
                                                     nTopEachList = nTopEachList,
                                                     nSelList = nSelList,
                                                     multiTraitsEvalMethodList = multiTraitsEvalMethodList,
                                                     hSelList = hSelList,
                                                     matingMethodVec = matingMethodVec,
                                                     allocateMethodVec = allocateMethodVec,
                                                     weightedAllocationMethodList = weightedAllocationMethodList,
                                                     traitNoRAList = traitNoRAList,
                                                     hList = hList,
                                                     includeGVPVec = includeGVPVec,
                                                     nNextPopVec = nNextPopVec,
                                                     nameMethod = nameMethod,
                                                     nCores =  nCores,
                                                     overWriteRes = overWriteRes,
                                                     showProgress = showProgress,
                                                     returnMethod = returnMethod,
                                                     saveAllResAt = saveAllResAt,
                                                     evaluateGVMethod = evaluateGVMethod,
                                                     nTopEval = nTopEval,
                                                     traitNoEval = traitNoEval,
                                                     hEval = hEval,
                                                     summaryAllResAt = summaryAllResAt,
                                                     verbose = verbose)


          ##### 3.2. Start simulation, summary results, and save simBs class object #####
          print("Start simulation....")

          checkFinishSimulation <- FALSE
          while (!checkFinishSimulation) {
            if (strategyName == "selectOPV") {
              mySimBsBVNow$nCores <- nCoresWhenOPV
            }

            if (strategyName == "selectEMBV") {
              mySimBsBVNow$nCores <- nCoresWhenEMBV
            }
            system.time(
              mySimBsBVNow$startSimulation()
            )

            if (!is.null(saveAllResAt)) {
              saveAllResAtSplit <- stringr::str_split(string = list.files(saveAllResAt),
                                                      pattern = "_")
              saveAllResAtSplitLast <- lapply(X = saveAllResAtSplit,
                                              FUN = function (saveAllResAtSplitVec) {
                                                return(saveAllResAtSplitVec[length(saveAllResAtSplitVec)])
                                              })

              saveAllNumeric <- unique(sort(as.numeric(stringr::str_remove(saveAllResAtSplitLast, ".rds"))))

              checkFinishSimulation <- all(unlist(lapply(mySimBsBVNow$trueGVMatList, length)) == (nGenerationProceed + 1)) |
                all((1:nIterSimulation) %in% saveAllNumeric)
            } else {
              checkFinishSimulation <- all(unlist(lapply(mySimBsBVNow$trueGVMatList, length)) == (nGenerationProceed + 1))
            }
          }

          if (strategyName %in% c("selectOPV", "selectEMBV")) {
            mySimBsBVNow$nCores <- nCores
          }

          if (all(unlist(lapply(mySimBsBVNow$trueGVMatList, length)) == (nGenerationProceed + 1))) {
            mySimBsBVNow$summaryResults()
          } else {
            print("Start summarization....")
            mySimBsBVNow$summaryAllResAt <- saveAllResAt
            mySimBsBVNow$nCores <- nCoresSummary
            mySimBsBVNow$summaryResults()
            mySimBsBVNow$nCores <- nCores
          }



          saveRDS(object = mySimBsBVNow, file = fileNameSimBsBVNow)
          # mySimBsBVList[[simBsName]] <- mySimBsBVNow


          rm(mySimBsBVNow)
          gc(reset = TRUE)
          ed <- Sys.time()

          print(paste0("Computation time required: ", round(ed - st, 2), " ",
                       attr(ed - st, "units")))

          # mySimEval <- myBreedSimulatR::simEval$new(simEvalName = simBsNameTrainingPop,
          #                                           simBsList = mySimBsBVList,
          #                                           verbose = verbose)
          #
          # fileNameSimEval <- paste0(dirMidSCOBSTrainingPop, scriptID, "_", simBsNameTrainingPop, "_simEval.rds")
          # saveRDS(object = mySimEval, file = fileNameSimEval)
          #
          # rm(mySimEval)
          gc(reset = TRUE)
        }
      }
    }
  }
}
