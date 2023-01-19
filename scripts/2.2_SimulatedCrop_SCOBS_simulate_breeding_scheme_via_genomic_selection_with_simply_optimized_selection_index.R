###################################################################################################################################
######  Title: 2.2_SimulatedCrop_SCOBS_simulate_breeding_scheme_via_genomic_selection_with_simply_optimized_selection_index  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                                     ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                                            ######
######  Date: 2020/12/26 (Created), 2021/01/15 (Last Updated)                                                                ######
###################################################################################################################################





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

scriptIDData <- "0.1"
scriptID <- "2.2"


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
trialNo <- 2
trialName <- paste0("Trial_", trialNo)

dirMidSCOBSTrial <- paste0(dirMidSCOBSScenario, scriptIDData,
                           "_", project, "_", trialName, "/")
if (!dir.exists(dirMidSCOBSTrial)) {
  dir.create(dirMidSCOBSTrial)
}

nSimGeno <- 1
nSimPheno <- 1

nIterSimulation <- 1000
nGenerationProceed <- 4
nGenerationProceedSimulation <- 4

nRefreshMemoryEvery <- 2
updateBreederInfo <- rep(TRUE, nGenerationProceed)
phenotypingInds <- rep(FALSE, nGenerationProceed)
nRepForPhenoInit <- 3
nRepForPheno <- rep(1, nGenerationProceed)
updateModels <- rep(FALSE, nGenerationProceed)

updateBreederInfoSimulation <- rep(TRUE, nGenerationProceedSimulation)
phenotypingIndsSimulation <- rep(FALSE, nGenerationProceedSimulation)
nRepForPhenoSimulation <- rep(1, nGenerationProceedSimulation)
updateModelsSimulation <- rep(FALSE, nGenerationProceedSimulation)

# strategyNames <- c("selectBV", "selectWBV",
#                    "selectOHV", "selectOPV")
strategyNames <- c(
  # "selectBV"
  # ,
  "selectWBV"
  # ,
  # "selectOHV"
)

nTrainingPopInit <- "true"

lociEffMethods <- c("true", "estimated")
trainingPopType <- "latest"
trainingPopInit <- 1
trainingIndNamesInit <- NULL


#### 1.2.3. Setting some parameters related to optimization of hyperparameters ####
setGoalAsFinalGeneration <- TRUE
performOptimization <- c(TRUE, rep(FALSE, nGenerationProceed - 1))
useFirstOptimizedValue <- TRUE
sameAcrossGeneration <- FALSE
hMin <- 0
hMax <- 2
nTotalIterForOneOptimization <- 1e06
nIterSimulationPerEvaluation <- 50
nIterOptimization <- 20000

simBsNameBase <- paste0("VerySimpleGS-SimplestOpt-",
                        nIterSimulationPerEvaluation, "_BS-",
                        nIterOptimization, "_Eval-for-Trait_1")

# simBsNameBase <- paste0("VerySimpleGS-AllocateBVOpt-",
#                         nIterSimulationPerEvaluation, "_BS-",
#                         nIterOptimization, "_Eval-for-Trait_1")

# simBsNameBase <- paste0("VerySimpleGS-IncludeGVPOpt-",
#                         nIterSimulationPerEvaluation, "_BS-",
#                         nIterOptimization, "_Eval-for-Trait_1")

# simBsNameBase <- paste0("VerySimpleGS-AllocateBVIncludeGVPOpt-",
#                         nIterSimulationPerEvaluation, "_BS-",
#                         nIterOptimization, "_Eval-for-Trait_1")

# simBsNameBase <- paste0("VerySimpleGS-IncludeTwoBVOpt-",
#                         nIterSimulationPerEvaluation, "_BS-",
#                         nIterOptimization, "_Eval-for-Trait_1")

# simBsNameBase <- paste0("VerySimpleGS-IncludeTwoBVGVPOpt-",
#                         nIterSimulationPerEvaluation, "_BS-",
#                         nIterOptimization, "_Eval-for-Trait_1")

nChildrenPerExpansion <- 3
nMaxEvalPerNode <- ceiling(x = nIterOptimization / (log(x = nIterOptimization) ^ 3))
maxDepth <- min(ceiling(x = sqrt(x = nIterOptimization / nMaxEvalPerNode)),
                floor(x = logb(x = 9e15, base = nChildrenPerExpansion)))
confidenceParam <- 1 / sqrt(x = nIterOptimization)
returnOptimalNodes <- seq(from = nIterOptimization / 100, to = nIterOptimization,
                          by = nIterOptimization / 100)
saveTreesAtBase <- "current_trees_backup"
whenToSaveTrees <- seq(from = nIterOptimization / 100, to = nIterOptimization,
                       by = nIterOptimization / 100)

nTopEvalForOpt <- 5
rewardWeightVec <- c(rep(0, nGenerationProceedSimulation - 1), 1)
# discountedRate <- 0.8
# rewardWeightVec <- sapply(1:nGenerationProceedSimulation,
#                           function(genProceedNo) discountedRate ^ (genProceedNo - 1))
digitsEval <- 3



#### 1.2.4. Setting some parameters related to estimation of marker effects ####
methodMLRInit <- "BayesB"
multiTraitInit <- FALSE
methodMLR <- "LASSO"
multiTrait <- FALSE

alpha <- 0.5
nIter <- 20000
burnIn <- 5000
thin <- 5
bayesian <- TRUE


#### 1.2.5. Setting some parameters related to selection of parent candidates ####
nSelectionWaysVec <- rep(1, nGenerationProceedSimulation)
traitNoSelList <- rep(list(list(1)), nGenerationProceedSimulation)

blockSplitMethod <- "minimumSegmentLength"
nMrkInBlock <- 10
minimumSegmentLength <- 100
nSelInitOPVList <- rep(list(50), nGenerationProceedSimulation)
nIterOPV <- 2e04
nProgeniesEMBVVec <- rep(40, nGenerationProceedSimulation)
nIterEMBV <- 5
nCoresEMBV <- 1

clusteringForSelList <- rep(list(FALSE), nGenerationProceedSimulation)
nClusterList <- rep(list(1), nGenerationProceedSimulation)
nTopClusterList <- rep(list(1), nGenerationProceedSimulation)
nTopEachList <- rep(list(15), nGenerationProceedSimulation)
nSelList <- rep(list(15), nGenerationProceedSimulation)

multiTraitsEvalMethodList <- rep(list("sum"), nGenerationProceedSimulation)
hSelList <- rep(list(list(1)), nGenerationProceedSimulation)


#### 1.2.6. Setting some parameters related to mating and resource allocation ####
matingMethodVec <- rep("diallelWithSelfing", nGenerationProceedSimulation)
allocateMethodVec <- rep("weightedAllocation", nGenerationProceedSimulation)
traitNoRAList <- rep(list(1), nGenerationProceedSimulation)
hList <- rep(list(1), nGenerationProceedSimulation)
if (!stringr::str_detect(string = simBsNameBase, pattern = "GVP")) {
  includeGVPVec <- rep(FALSE, nGenerationProceedSimulation)
} else {
  includeGVPVec <- rep(TRUE, nGenerationProceedSimulation)
}
nNextPopVec <- rep(250, nGenerationProceedSimulation)


#### 1.2.7. Setting some parameters related to mating and resource allocation ####
nCores <- 1
nCoresPerOptimization <- 50
nCoresPerOptimizationWhenOPV <- 50
nCoresPerOptimizationWhenEMBV <- 34
nCoresSummary <- 50

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


#### 1.2.8. Save parameters ####
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
for (simGenoNo in 1:nSimGeno) {
  genoName <- paste0("Geno_", simGenoNo)

  dirMidSCOBSGeno <- paste0(dirMidSCOBSTrial, scriptIDData,
                            "_", genoName, "/")


  for (simPhenoNo in 1:nSimPheno) {
    phenoName <- paste0("Pheno_", simPhenoNo)

    dirMidSCOBSPheno <- paste0(dirMidSCOBSGeno, scriptIDData,
                               "_", phenoName, "/")


    ###### 2. Load bsInfoInit & breederInfoInit & estimate marker effects ######
    simBsNameTrainingPop <- paste0(simBsNameBase, "-", nTrainingPopInit)

    dirMidSCOBSTrainingPop <- paste0(dirMidSCOBSPheno, scriptIDData,
                                     "_Model_with_", nTrainingPopInit, "/")
    if (!dir.exists(dirMidSCOBSTrainingPop)) {
      dir.create(dirMidSCOBSTrainingPop)
    }



    if (nTrainingPopInit == "true") {
      dirMidSCOBSTrainingPopFalse <- paste0(dirMidSCOBSPheno, scriptIDData,
                                            "_Model_with_", 2000, "/")
      fileNameBsInfoInitRmInit <- paste0(dirMidSCOBSTrainingPopFalse, scriptIDData,
                                         "_bsInfoInit_remove_initial_population.rds")
      fileNameBreederInfoInitWithMrkEffRmInit <- paste0(dirMidSCOBSTrainingPopFalse, scriptIDData,
                                                        "_breederInfoInit_with_estimated_marker_effects_",
                                                        methodMLRInit, "_", 2000,
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

    # if (file.exists(fileNameBreederInfoInitWithMrkEffRmInit)) {
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
    #   saveRDS(object = myBreederInfoInit, file = fileNameBreederInfoInitWithMrkEffRmInit)
    # }


    if (any(includeGVPVec)) {
      if (is.null(myBsInfoInit$lociInfo$recombBetweenMarkersList)) {
        myBsInfoInit$lociInfo$computeRecombBetweenMarkers()
        myBsInfoInit$traitInfo$lociInfo$computeRecombBetweenMarkers()
        myBsInfoInit$populations[[myBsInfoInit$generation]]$traitInfo$lociInfo$computeRecombBetweenMarkers()
      }
    }

    # mySimBsBVList <- list()
    for (strategyName in strategyNames) {
      print(paste0(genoName, "  ", phenoName, "  ",
                   strategyName, "  ", lociEffMethod))
      st <- Sys.time()
      simBsName <- paste0(simBsNameTrainingPop, "-", strategyName, "-", lociEffMethod)
      selectionMethodList <- rep(list(strategyName), nGenerationProceedSimulation)
      if (!any(stringr::str_detect(string = simBsNameBase, pattern = c("TwoBV", "AllocateBV", "AllocateWBV")))) {
        weightedAllocationMethodList <- rep(list(strategyName), nGenerationProceedSimulation)
      } else if (stringr::str_detect(string = simBsNameBase, pattern = "AllocateBV")) {
        weightedAllocationMethodList <- rep(list("selectBV"), nGenerationProceedSimulation)
      } else if (stringr::str_detect(string = simBsNameBase, pattern = "AllocateWBV")) {
        weightedAllocationMethodList <- rep(list("selectWBV"), nGenerationProceedSimulation)
      } else {
        weightedAllocationMethodList <- rep(list(strategyNames[1:2]), nGenerationProceedSimulation)
      }
      dirMidSCOBSBV <- paste0(dirMidSCOBSTrainingPop, scriptID, "_", simBsName, "/")
      if (!dir.exists(dirMidSCOBSBV)) {
        dir.create(dirMidSCOBSBV)
      }

      # if ((simGenoNo == 1) & (simPhenoNo == 1)) {
      #   saveAllResAt <- paste0(dirMidSCOBSBV, scriptID, "_", saveAllResAtBase)
      # } else {
      saveAllResAt <- NULL
      # }
      # summaryAllResAt <- paste0(dirMidSCOBSBV, scriptID, "_", summaryAllResAtBase)
      summaryAllResAt <- NULL
      saveTreesAt <- paste0(dirMidSCOBSBV, scriptID, "_", saveTreesAtBase, "/")
      if (!dir.exists(saveTreesAt)) {
        dir.create(saveTreesAt)
      }
      saveTreeNameBase <-  paste0(saveTreesAt, scriptID, "_", simBsName)


      fileNameSimBsBVNow <- paste0(dirMidSCOBSBV, scriptID, "_", simBsName, "_simBs.rds")

      if (!file.exists(fileNameSimBsBVNow)) {

        ###### 3. Simulate breeding scheme using simple GS with each BV ######
        ##### 3.1. Set simBs class object #####

        mySimBsBVNow <- myBreedSimulatR::simBsOpt$new(simBsName = simBsName,
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
                                                      nGenerationProceedSimulation = nGenerationProceedSimulation,
                                                      setGoalAsFinalGeneration = setGoalAsFinalGeneration,
                                                      performOptimization = performOptimization,
                                                      useFirstOptimizedValue = useFirstOptimizedValue,
                                                      sameAcrossGeneration = sameAcrossGeneration,
                                                      hMin = hMin,
                                                      hMax = hMax,
                                                      nTotalIterForOneOptimization = nTotalIterForOneOptimization,
                                                      nIterSimulationPerEvaluation = nIterSimulationPerEvaluation,
                                                      nIterOptimization = nIterOptimization,
                                                      nMaxEvalPerNode = nMaxEvalPerNode,
                                                      maxDepth = maxDepth,
                                                      nChildrenPerExpansion = nChildrenPerExpansion,
                                                      confidenceParam = confidenceParam,
                                                      returnOptimalNodes = returnOptimalNodes,
                                                      saveTreeNameBase = saveTreeNameBase,
                                                      whenToSaveTrees = whenToSaveTrees,
                                                      nTopEvalForOpt = nTopEvalForOpt,
                                                      rewardWeightVec = rewardWeightVec,
                                                      digitsEval = digitsEval,
                                                      nRefreshMemoryEvery = nRefreshMemoryEvery,
                                                      updateBreederInfo = updateBreederInfo,
                                                      phenotypingInds = phenotypingInds,
                                                      nRepForPhenoInit = nRepForPhenoInit,
                                                      nRepForPheno = nRepForPheno,
                                                      updateModels = updateModels,
                                                      updateBreederInfoSimulation = updateBreederInfoSimulation,
                                                      phenotypingIndsSimulation = phenotypingIndsSimulation,
                                                      nRepForPhenoSimulation = nRepForPhenoSimulation,
                                                      updateModelsSimulation = updateModelsSimulation,
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
                                                      includeGVPVec = includeGVPVec,
                                                      nNextPopVec = nNextPopVec,
                                                      nameMethod = nameMethod,
                                                      nCores = nCores,
                                                      nCoresPerOptimization = nCoresPerOptimization,
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
            mySimBsBVNow$nCoresPerOptimization <- nCoresPerOptimizationWhenOPV
          }

          if (strategyName == "selectEMBV") {
            mySimBsBVNow$nCoresPerOptimization <- nCoresPerOptimizationWhenEMBV
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
          mySimBsBVNow$nCoresPerOptimization <- nCoresPerOptimization
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

      }
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
