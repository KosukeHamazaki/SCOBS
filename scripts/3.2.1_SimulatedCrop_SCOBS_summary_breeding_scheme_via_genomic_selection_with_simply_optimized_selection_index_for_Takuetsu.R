#####################################################################################################################################################
######  Title: 3.2.1_SimulatedCrop_SCOBS_summary_breeding_scheme_via_genomic_selection_with_simply_optimized_selection_index_for_Takuetsu      ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                                                       ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                                                              ######
######  Date: 2021/02/26 (Created), 2021/02/26 (Last Updated)                                                                                  ######
#####################################################################################################################################################





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
scriptIDMainOpt <- "2.2"
scriptIDSummaryOpt <- "3.2"
scriptID <- "3.2.1"

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

# nSimGeno <- 1
# nSimPheno <- 1
#
# nIterSimulation <- 1000
# nGenerationProceed <- 4
# nGenerationProceedSimulation <- 4
#
# nRefreshMemoryEvery <- 2
# updateBreederInfo <- rep(TRUE, nGenerationProceed)
# phenotypingInds <- rep(FALSE, nGenerationProceed)
# nRepForPhenoInit <- 3
# nRepForPheno <- rep(1, nGenerationProceed)
# updateModels <- rep(FALSE, nGenerationProceed)
#
# updateBreederInfoSimulation <- rep(TRUE, nGenerationProceedSimulation)
# phenotypingIndsSimulation <- rep(FALSE, nGenerationProceedSimulation)
# nRepForPhenoSimulation <- rep(1, nGenerationProceedSimulation)
# updateModelsSimulation <- rep(FALSE, nGenerationProceedSimulation)
#
# # strategyNames <- c("selectBV", "selectWBV",
# #                    "selectOHV", "selectOPV")
# strategyNames <- c("selectBV", "selectWBV",
#                    "selectOHV")
#
# nTrainingPopInit <- "true"
#
# lociEffMethods <- c("true", "estimated")
# trainingPopType <- "latest"
# trainingPopInit <- 1
# trainingIndNamesInit <- NULL
#
#
# #### 1.2.3. Setting some parameters related to optimization of hyperparameters ####
# setGoalAsFinalGeneration <- TRUE
# performOptimization <- c(TRUE, rep(FALSE, nGenerationProceed - 1))
# useFirstOptimizedValue <- TRUE
# sameAcrossGeneration <- FALSE
# hMin <- 0
# hMax <- 2
# nTotalIterForOneOptimization <- 1e06
nIterSimulationPerEvaluation <- 100
nIterOptimization <- 10000

simBsNameBase <- paste0("VerySimpleGS-SimplestOpt-",
                        nIterSimulationPerEvaluation, "_BS-",
                        nIterOptimization, "_Eval-for-Trait_1")

# nChildrenPerExpansion <- 3
# nMaxEvalPerNode <- ceiling(x = nIterOptimization / (log(x = nIterOptimization) ^ 3))
# maxDepth <- min(ceiling(x = sqrt(x = nIterOptimization / nMaxEvalPerNode)),
#                 floor(x = logb(x = 9e15, base = nChildrenPerExpansion)))
# confidenceParam <- 1 / sqrt(x = nIterOptimization)
# returnOptimalNodes <- seq(from = 100, to = nIterOptimization,
#                           by = 100)
#
# nTopEvalForOpt <- 5
# rewardWeightVec <- c(rep(0, nGenerationProceedSimulation - 1), 1)
# # discountedRate <- 0.8
# # rewardWeightVec <- sapply(1:nGenerationProceedSimulation,
# #                           function(genProceedNo) discountedRate ^ (genProceedNo - 1))
# digitsEval <- 3
#
#
#
# #### 1.2.4. Setting some parameters related to estimation of marker effects ####
# methodMLRInit <- "BayesB"
# multiTraitInit <- FALSE
# methodMLR <- "LASSO"
# multiTrait <- FALSE
#
# alpha <- 0.5
# nIter <- 20000
# burnIn <- 5000
# thin <- 5
# bayesian <- TRUE
#
#
# #### 1.2.5. Setting some parameters related to selection of parent candidates ####
# nSelectionWaysVec <- rep(1, nGenerationProceedSimulation)
# traitNoSelList <- rep(list(list(1)), nGenerationProceedSimulation)
#
# blockSplitMethod <- "minimumSegmentLength"
# nMrkInBlock <- 10
# minimumSegmentLength <- 100
# nSelInitOPVList <- rep(list(50), nGenerationProceedSimulation)
# nIterOPV <- 2e04
# nProgeniesEMBVVec <- rep(40, nGenerationProceedSimulation)
# nIterEMBV <- 5
# nCoresEMBV <- 1
#
# clusteringForSelList <- rep(list(FALSE), nGenerationProceedSimulation)
# nClusterList <- rep(list(1), nGenerationProceedSimulation)
# nTopClusterList <- rep(list(1), nGenerationProceedSimulation)
# nTopEachList <- rep(list(15), nGenerationProceedSimulation)
# nSelList <- rep(list(15), nGenerationProceedSimulation)
#
# multiTraitsEvalMethodList <- rep(list("sum"), nGenerationProceedSimulation)
# hSelList <- rep(list(list(1)), nGenerationProceedSimulation)
#
#
# #### 1.2.6. Setting some parameters related to mating and resource allocation ####
# matingMethodVec <- rep("diallelWithSelfing", nGenerationProceedSimulation)
# allocateMethodVec <- rep("weightedAllocation", nGenerationProceedSimulation)
# traitNoRAList <- rep(list(1), nGenerationProceedSimulation)
# hList <- rep(list(1), nGenerationProceedSimulation)
# includeGVPVec <- rep(FALSE, nGenerationProceedSimulation)
# nNextPopVec <- rep(250, nGenerationProceedSimulation)
#
#
# #### 1.2.7. Setting some parameters related to mating and resource allocation ####
# nCores <- 1
# nCoresPerOptimization <- 50
# nCoresPerOptimizationWhenOPV <- 50
# nCoresPerOptimizationWhenEMBV <- 34
# nCoresSummary <- 50
#
# nameMethod <- "pairBase"
# overWriteRes <- FALSE
# showProgress <- TRUE
# returnMethod <- "summary"
# evaluateGVMethod <- "true"
# nTopEval <- 5
# traitNoEval <- 1
# hEval <- 1
# verbose <- TRUE
#
# saveAllResAtBase <- "all_results"
# summaryAllResAtBase <- "all_results"


#### 1.2.8. Save parameters ####
fileParamsSCOBS <- paste0(dirMidSCOBSTrial, scriptIDMainOpt,
                          "_", project, "_", trialName,
                          "_", simBsNameBase,
                          "_all_parameters.RData")
load(fileParamsSCOBS)

nTrainingPopInit <- 2000

strategyNames <- c(
  "selectBV", "selectWBV"
  # , "selectOHV"
  # , "selectOPV"
  # , "selectEMBV"
)
scriptIDData <- "0.1"
scriptIDMainOpt <- "2.2"
scriptIDSummaryOpt <- "3.2"
scriptID <- "3.2.1"

nCores <-  10
nCoresWhenOPV <-  10
nCoresWhenEMBV <- 10
nCoresSummary <- 10

nIterSimulation <- 10000

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


simGenoNo <- 1
simPhenoNo <- 1
strategyName <- strategyNames[1]
lociEffMethod <- lociEffMethods[1]
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
      lociEffMethod <- lociEffMethods[1]
    } else {
      lociEffMethod <- lociEffMethods[2]
    }

    mySimBsBVConvListAll <- NULL
    for (strategyName in strategyNames) {
      print(paste0(genoName, "  ", phenoName, "  ",
                   strategyName, "  ", lociEffMethod))
      simBsName <- paste0(simBsNameTrainingPop, "-", strategyName, "-", lociEffMethod)
      selectionMethodList <- rep(list(strategyName), nGenerationProceedSimulation)
      weightedAllocationMethodList <- rep(list(strategyName), nGenerationProceedSimulation)

      dirMidSCOBSBV <- paste0(dirMidSCOBSTrainingPop, scriptIDMainOpt, "_", simBsName, "/")
      if (!dir.exists(dirMidSCOBSBV)) {
        dir.create(dirMidSCOBSBV)
      }

      # if ((simGenoNo == 1) & (simPhenoNo == 1)) {
      #   saveAllResAt <- paste0(dirMidSCOBSBV, scriptIDMainOpt, "_", saveAllResAtBase)
      # } else {
      saveAllResAt <- NULL
      # }
      # summaryAllResAt <- paste0(dirMidSCOBSBV, scriptIDMainOpt, "_", summaryAllResAtBase)
      summaryAllResAt <- NULL

      fileNameSimEval <- paste0(dirMidSCOBSBV, scriptIDSummaryOpt, "_", simBsName,
                                "_convergence_simEval.rds")

      if (file.exists(fileNameSimEval)) {
        mySimEval <- readRDS(file = fileNameSimEval)
        mySimEval$plot(targetTrait = 1,
                       targetPopulation = 1:5,
                       plotType = "lines",
                       plotTarget = "max",
                       returnGain = TRUE)
        mySimBsBVConvList <- mySimEval$clone(deep = FALSE)$simBsList
        nUniq <- length(mySimBsBVConvList) - 1
        uniqIterationNames <- unlist(lapply(
          X = stringr::str_split(string = names(mySimBsBVConvList[1:nUniq]),
                                 pattern = "-"),
          FUN = function(x) x[length(x)]
        ))

        trueGVSummaryArrayForConv <- do.call(
          what = abind::abind,
          args = lapply(X = mySimBsBVConvList,
                        FUN = function(mySimBsBV) {
                          apply(X = mySimBsBV$trueGVSummaryArray,
                                MARGIN = c(1, 3, 2),
                                FUN = mean)
                        })
        )

        dimnames(trueGVSummaryArrayForConv)[[3]] <- c(uniqIterationNames, "Equal")

        whichMaxEst <- which.max(trueGVSummaryArrayForConv[1, 5, ])
        finalEst <- nUniq


        mySimBsBVConvListNow <- mySimBsBVConvList[unique(c(whichMaxEst, finalEst, nUniq + 1))]
        mySimBsBVConvListAll <- c(mySimBsBVConvListAll, mySimBsBVConvListNow)
      }
    }
    selectionMethodNamesNow <- unlist(x = lapply(X = mySimBsBVConvListAll,
                                                 FUN = function(x) {
                                                   x$selectionMethodList[[1]]
                                                 }), use.names = FALSE)
    nBV <- sum(selectionMethodNamesNow %in% "selectBV")
    nWBV <- sum(selectionMethodNamesNow %in% "selectWBV")
    if (nBV == 2) {
      methodNamesBVOptsBV <- paste0("BV-", c("Max=Opt", "Equal"))
    } else {
      methodNamesBVOptsBV <- paste0("BV-", c("Max", "Opt", "Equal"))
    }

    if (nWBV == 2) {
      methodNamesBVOptsWBV <- paste0("WBV-", c("Max=Opt", "Equal"))
    } else {
      methodNamesBVOptsWBV <- paste0("WBV-", c("Max", "Opt", "Equal"))

    }
    methodNamesBVOpts <- c(methodNamesBVOptsBV, methodNamesBVOptsWBV)
    mySimBsBVConvListAll <- lapply(X = 1:length(mySimBsBVConvListAll),
                                   FUN = function(x) {
                                     mySimBsBVConvListAll[[x]]$simBsName <- methodNamesBVOpts[x]

                                     return(mySimBsBVConvListAll[[x]])
                                   })
    names(mySimBsBVConvListAll) <- methodNamesBVOpts


    simEvalOptMethodName <- paste0(simBsNameBase, "-", nTrainingPopInit,
                                   "-BV_WBV-Opt_vs_Equal")
    mySimEvalOptMethod <- myBreedSimulatR::simEval$new(simEvalName = simEvalOptMethodName,
                                                       simBsList = mySimBsBVConvListAll)
    fileNameSimEvalOpt <- paste0(dirMidSCOBSTrainingPop, scriptID, "_", simEvalOptMethodName,
                                 "-BV_WBV", "_simEval.rds")
    saveRDS(object = mySimEvalOptMethod,
            file = fileNameSimEvalOpt)

    pltTwoOptOrigin <- mySimEvalOptMethod$plot(targetTrait = 1,
                                               targetPopulation = 1:5,
                                               plotType = "lines",
                                               plotTarget = "max",
                                               returnGain = TRUE)
    pltTwoOptOrigin <- pltTwoOptOrigin %>% plotly::layout(showlegend = TRUE,
                                                          xaxis = list(tickvals = as.list(0:4),
                                                                       ticktext = as.list(paste0("Population_", 1:5)),
                                                                       tickmode = "array"))


    fileNamePltTwoOptMaxPng <- paste0(dirMidSCOBSTrainingPop, scriptID, "_", simEvalOptMethodName,
                                      "-BV_WBV", "_plotly_lines_max.png")
    fileNamePltTwoOptMaxHtml <- paste0(dirMidSCOBSTrainingPop, scriptID, "_", simEvalOptMethodName,
                                       "-BV_WBV", "_plotly_lines_max.html")
    plotly::orca(p = pltTwoOptOrigin,
                 format = "png",
                 width = 1200,
                 height = 800,
                 file = fileNamePltTwoOptMaxPng)

    htmlwidgets::saveWidget(widget = plotly::partial_bundle(pltTwoOptOrigin),
                            file = here::here(fileNamePltTwoOptMaxHtml))
  }
}
