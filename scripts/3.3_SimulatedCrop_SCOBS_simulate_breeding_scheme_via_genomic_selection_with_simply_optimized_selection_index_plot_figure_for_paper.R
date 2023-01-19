#########################################################################################################################################################
######  Title: 3.3_SimulatedCrop_SCOBS_simulate_breeding_scheme_via_genomic_selection_with_simply_optimized_selection_index_plot_figure_for_paper  ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                                                           ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                                                                  ######
######  Date: 2022/05/06 (Created), 2022/05/16 (Last Updated)                                                                                      ######
#########################################################################################################################################################





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
scriptIDConv <- "3.2"
scriptID <- "3.3"



##### 1.2. Import packages #####
require(myBreedSimulatR)
require(data.table)
require(MASS)
require(rrBLUP)
require(BGLR)
require(RAINBOWR)
require(ggplot2)
require(plotly)





##### 1.3. Project options #####
options(stringAsFactors = FALSE)





##### 1.4. Setting some parameters #####
dirMidSCOBSBase <- "midstream/"
overWriteRes <- FALSE


#### 1.4.1. Setting some parameters related to names of R6 class ####
scenarioNo <- 1
simName <- paste0("Scenario_", scenarioNo)
breederName <- "K. Hamazaki"

dirMidSCOBSScenario <- paste0(dirMidSCOBSBase, scriptIDData, "_", simName, "/")
if (!dir.exists(dirMidSCOBSScenario)) {
  dir.create(dirMidSCOBSScenario)
}


#### 1.4.2. Setting some parameters related to simulation ####
trialNo <- 2
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
# #### 1.4.3. Setting some parameters related to optimization of hyperparameters ####
# setGoalAsFinalGeneration <- TRUE
# performOptimization <- c(TRUE, rep(FALSE, nGenerationProceed - 1))
# useFirstOptimizedValue <- TRUE
# sameAcrossGeneration <- FALSE
# hMin <- 0
# hMax <- 2
# nTotalIterForOneOptimization <- 1e06
nIterSimulationPerEvaluation <- 50
nIterOptimization <- 20000


allocateStrategyNames <- c("SimplestOpt",
                           "AllocateBVOpt",
                           "IncludeGVPOpt",
                           "AllocateBVIncludeGVPOpt",
                           "IncludeTwoBVOpt",
                           "IncludeTwoBVGVPOpt")
# allocateStrategyName <- allocateStrategyNames[1]

geno1Name <- "Geno_1"
dirMidSCOBSGeno1 <- paste0(dirMidSCOBSTrial, scriptIDData,
                           "_", geno1Name, "/")
pheno1Name <- "Pheno_1"

dirMidSCOBSPheno1 <- paste0(dirMidSCOBSGeno1, scriptIDData,
                            "_", pheno1Name, "/")

dirMidSCOBSTrainingPop1 <- paste0(dirMidSCOBSPheno1, scriptIDData,
                                  "_Model_with_true/")


fileNameGeneticGainList <- paste0(dirMidSCOBSTrainingPop1, scriptID, "_",
                                  "all_change_in_genetic_gain_for_7_strategies.rds")
fileNameGeneticVarList <- paste0(dirMidSCOBSTrainingPop1, scriptID, "_",
                                 "all_change_in_genetic_variance_for_7_strategies.rds")
fileNameGeneticGainIterList <- paste0(dirMidSCOBSTrainingPop1, scriptID, "_",
                                      "final_genetic_gain_of_10000_iterations_for_7_strategies.rds")
fileNameGeneticGainFinalTenList <- paste0(dirMidSCOBSGeno1, scriptID, "_",
                                          "genetic_gain_for_ten_traits_for_7_strategies.rds")

if (nIterOptimization == 4000) {
  nSimPheno <- 10
  fileExists <- file.exists(fileNameGeneticGainFinalTenList)
} else {
  nSimPheno <- 1
  fileExists <- file.exists(fileNameGeneticGainList)
}



if (overWriteRes | (!fileExists)) {
  geneticGainList <- list()
  geneticVarList <- list()
  geneticGainIterList <- list()
  geneticGainFinalTenList <- list()
  for (allocateStrategyName in allocateStrategyNames) {
    if ((nIterOptimization == 20000) & (allocateStrategyName == "IncludeTwoBVGVPOpt")) {
      nIterOptimization <- 50000
    } else if ((nIterOptimization == 50000) & (allocateStrategyName != "IncludeTwoBVGVPOpt")) {
      nIterOptimization <- 20000
    }


    simBsNameBase <- paste0("VerySimpleGS-", allocateStrategyName, "-",
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
    # #### 1.4.4. Setting some parameters related to estimation of marker effects ####
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
    # #### 1.4.5. Setting some parameters related to selection of parent candidates ####
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
    # #### 1.4.6. Setting some parameters related to mating and resource allocation ####
    # matingMethodVec <- rep("diallelWithSelfing", nGenerationProceedSimulation)
    # allocateMethodVec <- rep("weightedAllocation", nGenerationProceedSimulation)
    # traitNoRAList <- rep(list(1), nGenerationProceedSimulation)
    # hList <- rep(list(1), nGenerationProceedSimulation)
    # includeGVPVec <- rep(FALSE, nGenerationProceedSimulation)
    # nNextPopVec <- rep(250, nGenerationProceedSimulation)
    #
    #
    # #### 1.4.7. Setting some parameters related to mating and resource allocation ####
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


    #### 1.4.8. Load parameters ####
    fileParamsSCOBS <- paste0(dirMidSCOBSTrial, scriptIDMainOpt,
                              "_", project, "_", trialName,
                              "_", simBsNameBase,
                              "_all_parameters.RData")
    load(fileParamsSCOBS)

    strategyNames <- c(
      # "selectBV"
      # ,
      "selectWBV"
      # , "selectOHV"
      # , "selectOPV"
      # , "selectEMBV"
    )
    scriptIDData <- "0.1"
    scriptIDMainOpt <- "2.2"
    scriptIDConv <- "3.2"
    scriptID <- "3.3"


    nCores <-  50
    nCoresWhenOPV <-  50
    nCoresWhenEMBV <- 34
    nCoresSummary <- 34

    nIterSimulation <- 10000
    nTrainingPopInit <- "true"



    simGenoNo <- 1
    simPhenoNo <- 1
    lociEffMethod <- lociEffMethods[1]







    ###### 2. Start summarization ######
    ##### 2.1. Loop for genome #####
    for (simGenoNo in 1:nSimGeno) {
      genoName <- paste0("Geno_", simGenoNo)

      dirMidSCOBSGeno <- paste0(dirMidSCOBSTrial, scriptIDData,
                                "_", genoName, "/")


      ##### 2.2. Loop for phenotype #####
      geneticGainFinalTen <- list()
      for (simPhenoNo in 1:nSimPheno) {
        phenoName <- paste0("Pheno_", simPhenoNo)

        dirMidSCOBSPheno <- paste0(dirMidSCOBSGeno, scriptIDData,
                                   "_", phenoName, "/")


        simBsNameTrainingPop <- paste0(simBsNameBase, "-", nTrainingPopInit)

        dirMidSCOBSTrainingPop <- paste0(dirMidSCOBSPheno, scriptIDData,
                                         "_Model_with_", nTrainingPopInit, "/")
        if (!dir.exists(dirMidSCOBSTrainingPop)) {
          dir.create(dirMidSCOBSTrainingPop)
        }



        if (nTrainingPopInit == "true") {
          # dirMidSCOBSTrainingPopFalse <- paste0(dirMidSCOBSPheno, scriptIDData,
          #                                       "_Model_with_", 2000, "/")
          # fileNameBsInfoInitRmInit <- paste0(dirMidSCOBSTrainingPopFalse, scriptIDData,
          #                                    "_bsInfoInit_remove_initial_population.rds")
          # fileNameBreederInfoInitWithMrkEffRmInit <- paste0(dirMidSCOBSTrainingPopFalse, scriptIDData,
          #                                                   "_breederInfoInit_with_estimated_marker_effects_",
          #                                                   methodMLRInit, "_", 2000,
          #                                                   "_individuals_remove_initial_population.rds")
          #
          lociEffMethod <- lociEffMethods[1]
        } else {
          # fileNameBsInfoInitRmInit <- paste0(dirMidSCOBSTrainingPop, scriptIDData,
          #                                    "_bsInfoInit_remove_initial_population.rds")
          # fileNameBreederInfoInitWithMrkEffRmInit <- paste0(dirMidSCOBSTrainingPop, scriptIDData,
          #                                                   "_breederInfoInit_with_estimated_marker_effects_",
          #                                                   methodMLRInit, "_", nTrainingPopInit,
          #                                                   "_individuals_remove_initial_population.rds")
          #
          lociEffMethod <- lociEffMethods[2]
        }
        # myBsInfoInit <- readRDS(file = fileNameBsInfoInitRmInit)

        # if (file.exists(fileNameBreederInfoInitWithMrkEffRmInit)) {
        # myBreederInfoInit <- readRDS(file = fileNameBreederInfoInitWithMrkEffRmInit)
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


        # if (any(includeGVPVec)) {
        #   if (is.null(myBsInfoInit$lociInfo$recombBetweenMarkersList)) {
        #     myBsInfoInit$lociInfo$computeRecombBetweenMarkers()
        #     myBsInfoInit$traitInfo$lociInfo$computeRecombBetweenMarkers()
        #     myBsInfoInit$populations[[myBsInfoInit$generation]]$traitInfo$lociInfo$computeRecombBetweenMarkers()
        #   }
        # }



        ##### 2.3. Loop for Selection strategy #####
        strategyName <- strategyNames[1]
        # for (strategyName in strategyNames) {
        print(paste0(allocateStrategyName, "  ", genoName, "  ", phenoName, "  ",
                     strategyName, "  ", lociEffMethod))
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


        ##### 2.4. Load simulation results #####
        fileNameSimEval <- paste0(dirMidSCOBSBV, scriptIDConv, "_", simBsName,
                                  "_convergence_simEval.rds")
        if (file.exists(fileNameSimEval)) {
          mySimEval <- readRDS(file = fileNameSimEval)

          trueGVSummaryArrayList <- lapply(X = mySimEval$simBsList,
                                           FUN = function(simBsEach) {
                                             trueGVSummaryArray <- simBsEach$trueGVSummaryArray
                                             trueGVSummaryArrayNeeded <- trueGVSummaryArray[c(1, 5), , , ]

                                             return(trueGVSummaryArrayNeeded)
                                           })




          ##### 2.5. Convergence check #####
          #### 2.5.1. Convergence plot #####
          trueGVSummaryConv <- do.call(what = rbind,
                                       args = lapply(X = trueGVSummaryArrayList,
                                                     FUN = function(trueGVSummaryArrayNeeded) {
                                                       apply(X = trueGVSummaryArrayNeeded[1, , ],
                                                             MARGIN = 1, mean, na.rm = TRUE)
                                                     }))

          nIterConv <- unlist(lapply(
            X = stringr::str_split(string = rownames(trueGVSummaryConv)[-nrow(trueGVSummaryConv)],
                                   pattern = "_"),
            FUN = function(x) {
              as.numeric(x[length(x)])
            }
          ))

          nIterConvDiff <- diff(c(nIterConv, nIterOptimization + 1))
          trueGVFinalGenConv <- rep(trueGVSummaryConv[-nrow(trueGVSummaryConv),
                                                      ncol(trueGVSummaryConv)] -
                                      trueGVSummaryConv[-nrow(trueGVSummaryConv), 1],
                                    nIterConvDiff)


          dfForGGPlot <- data.frame(Iteration = seq(from = 0, to = nIterOptimization, by = 1),
                                    Value = trueGVFinalGenConv)
          rownames(dfForGGPlot) <- 1:nrow(dfForGGPlot)


          pltConv <- ggplot2::ggplot(data = dfForGGPlot) +
            ggplot2::geom_line(mapping = aes(x = Iteration,
                                             y = Value)) +
            ggplot2::geom_hline(yintercept = trueGVSummaryConv[nrow(trueGVSummaryConv), ncol(trueGVSummaryConv)] -
                                  trueGVSummaryConv[nrow(trueGVSummaryConv), 1],
                                color = "red", linetype = "dashed") +
            ggplot2::xlim(0, min(nIterOptimization, 20000)) +
            ggplot2::theme(text = element_text(size = 48)) +
            ggplot2::ylab(label = "Genetic gain")
          fileNamePltConv <- paste0(dirMidSCOBSBV, scriptID, "_", simBsName,
                                    "_true_convergence_evaluation.png")

          png(filename = fileNamePltConv, width = 1200, height = 1000)
          print(pltConv)
          dev.off()



          #### 2.5.2. Optimized parameters #####
          optimizedParamsList <- lapply(X = mySimEval$simBsList,
                                        FUN = function(simBsEach) {
                                          hList <- simBsEach$hList
                                          optimizedParams <- do.call(what = rbind,
                                                                     args = hList)
                                          return(optimizedParams)
                                        })
          optimizedParamsList[[length(optimizedParamsList)]][] <- 0

          fileNameOptParams <- paste0(dirMidSCOBSBV, scriptID, "_", names(optimizedParamsList),
                                      "_optimized_parameteres.csv")

          saveOptimizedParams <- sapply(X = 1:length(optimizedParamsList),
                                        FUN = function(optimizedParamNo) {
                                          write.csv(x = optimizedParamsList[[optimizedParamNo]],
                                                    file = fileNameOptParams[optimizedParamNo])
                                        })



          ##### 2.6. Change over four generations #####
          #### 2.6.1. Change in genetic gain ####
          trueGVSummaryMeanList <- lapply(X = trueGVSummaryArrayList,
                                          FUN = function(trueGVSummaryArrayNeeded) {
                                            trueGVSummaryMean <- apply(X = trueGVSummaryArrayNeeded,
                                                                       MARGIN = c(1, 2), FUN = mean)

                                            return(trueGVSummaryMean)
                                          })
          geneticGainMat <- do.call(
            what = rbind,
            args = lapply(X = trueGVSummaryMeanList,
                          FUN = function(trueGVSummaryMean) {
                            trueGVSummaryMeanMax <- trueGVSummaryMean[1, ]
                            trueGVSummaryMeanGenGain <- trueGVSummaryMeanMax - trueGVSummaryMeanMax[1]

                            return(trueGVSummaryMeanGenGain)
                          })
          )


          if (nSimPheno > 1) {
            geneticGainFinalTen[[phenoName]] <- geneticGainMat[, ncol(geneticGainMat)]
          }




          if (nSimPheno == 1) {
            #### 2.6.2. Change in genetic variance ####
            geneticVarMat <- do.call(
              what = rbind,
              args = lapply(X = trueGVSummaryMeanList,
                            FUN = function(trueGVSummaryMean) {
                              trueGVSummaryMeanVar <- trueGVSummaryMean[2, ]

                              return(trueGVSummaryMeanVar)
                            })
            )





            ##### 2.7. Evaluation in final generation #####
            #### 2.7.1. Evaluation by density ####
            geneticGainIter <- do.call(what = rbind,
                                       args = lapply(X = trueGVSummaryArrayList,
                                                     FUN = function(trueGVSummaryArrayNeeded) {
                                                       trueGVSummaryArrayNeeded[1, ncol(trueGVSummaryArrayNeeded), ] -
                                                         trueGVSummaryArrayNeeded[1, 1, ]
                                                     }))
          }
        }





        # }
      }
    }


    if (nSimPheno == 1) {
      geneticGainList[[allocateStrategyName]] <- geneticGainMat
      geneticVarList[[allocateStrategyName]] <- geneticVarMat
      geneticGainIterList[[allocateStrategyName]] <- geneticGainIter
    } else {
      geneticGainFinalTenList[[allocateStrategyName]] <- geneticGainFinalTen
    }
  }
} else {
  if (nSimPheno == 1) {
    geneticGainList <- readRDS(file = fileNameGeneticGainList)
    geneticVarList <- readRDS(file = fileNameGeneticVarList)
    geneticGainIterList <- readRDS(file = fileNameGeneticGainIterList)
  } else {
    geneticGainFinalTenList <- readRDS(file = fileNameGeneticGainFinalTenList)
  }
}







##### 3. Draw plots for figures #####
##### 3.1. Some definitions #####
allocateStrategyNamesAbb <- c("WBV", "BV", "WBVGVP",
                              "BVGVP", "TBV", "TBVGVP")
allocateStrategyNamesAbbAllOrd <- c("EQ", "BV", "WBV", "TBV",
                                    "BVGVP", "WBVGVP", "TBVGVP")
if (nSimPheno == 1) {
  # allocateStrategyCols <- c("black", "yellow", "blue", "green",
  #                           "orange", "purple", "brown")
  allocateStrategyCols <- c("grey40", "goldenrod1", "blue", "violetred",
                            "goldenrod1", "blue", "violetred")
  allocateStrategyLtys <- c(1, 2, 2, 2, 1, 1, 1)
  allocateStrategyLwds <- 1.75
} else {
  # allocateStrategyCols <- c("grey80", "yellow", "lightblue", "lightgreen",
  #                           "orange1", "mediumpurple1", "rosybrown1")
  allocateStrategyCols <- c("grey80", "yellow", "lightblue", "rosybrown1",
                            "yellow", "lightblue", "rosybrown1")
}




##### 3.2. Results for one trait #####
if (nSimPheno == 1) {
  names(geneticGainList) <- names(geneticGainIterList) <-
    names(geneticVarList) <- allocateStrategyNamesAbb


  if (overWriteRes | (!fileExists)) {
    saveRDS(object = geneticGainList,
            file = fileNameGeneticGainList)
    saveRDS(object = geneticVarList,
            file = fileNameGeneticVarList)
    saveRDS(object = geneticGainIterList,
            file = fileNameGeneticGainIterList)
  }



  #### 3.2.1. Change in genetic gain ####
  geneticGainOptMat <- do.call(
    what = rbind,
    args = lapply(X = geneticGainList,
                  FUN = function(geneticGainStr) {
                    geneticGainStr[(nrow(geneticGainStr) - 1), ]
                  })
  )
  geneticGainEqMat <- do.call(
    what = rbind,
    args = lapply(X = geneticGainList,
                  FUN = function(geneticGainStr) {
                    geneticGainStr[nrow(geneticGainStr), ]
                  })
  )
  geneticGainEq <- apply(X = geneticGainEqMat,
                         MARGIN = 2, FUN = mean)
  geneticGainMatAll <- rbind(geneticGainOptMat, EQ = geneticGainEq)
  geneticGainMatAllOrd <- geneticGainMatAll[allocateStrategyNamesAbbAllOrd, ]

  geneticGainDf <- data.frame(Strategy = rep(allocateStrategyNamesAbbAllOrd,
                                             ncol(geneticGainMatAllOrd)),
                              Population = rep(colnames(geneticGainMatAllOrd),
                                               each = nrow(geneticGainMatAllOrd)),
                              Value = c(geneticGainMatAllOrd))
  geneticGainDf$Strategy <- factor(geneticGainDf$Strategy,
                                   levels = allocateStrategyNamesAbbAllOrd)
  # geneticGainDf$Population <- factor(geneticGainDf$Population,
  #                                    labels = paste0("Population-", 1:5))
  geneticGainDf$Population <- factor(geneticGainDf$Population,
                                     labels = paste0(0:4))

  pltChangeGain <- ggplot2::ggplot(data = geneticGainDf) +
    ggplot2::geom_line(mapping = aes(x = Population, y = Value,
                                     group = Strategy,
                                     color = Strategy),
                       lty = rep(allocateStrategyLtys, each = 5),
                       lwd = allocateStrategyLwds) +
    ggplot2::scale_color_manual(values = allocateStrategyCols) +
    ggplot2::theme(text = element_text(size = 96)) +
    ggplot2::xlab(label = "Generation number") +
    ggplot2::ylab(label = "Genetic gain")

  fileNamePltChangeGain <- paste0(dirMidSCOBSTrainingPop1, scriptID, "_",
                                  "change_in_genetic_gain_over_four_generations_for_7_strategies_fp.png")
  # png(filename = fileNamePltChangeGain, width = 1800, height = 1100)
  png(filename = fileNamePltChangeGain, width = 1600, height = 1600)
  print(pltChangeGain)
  dev.off()




  #### 3.2.2. Change in genetic variance ####
  geneticVarOptMat <- do.call(
    what = rbind,
    args = lapply(X = geneticVarList,
                  FUN = function(geneticVarStr) {
                    geneticVarStr[(nrow(geneticVarStr) - 1), ]
                  })
  )
  geneticVarEqMat <- do.call(
    what = rbind,
    args = lapply(X = geneticVarList,
                  FUN = function(geneticVarStr) {
                    geneticVarStr[nrow(geneticVarStr), ]
                  })
  )
  geneticVarEq <- apply(X = geneticVarEqMat,
                        MARGIN = 2, FUN = mean)
  geneticVarMatAll <- rbind(geneticVarOptMat, EQ = geneticVarEq)
  geneticVarMatAllOrd <- geneticVarMatAll[allocateStrategyNamesAbbAllOrd, ]

  geneticVarDf <- data.frame(Strategy = rep(allocateStrategyNamesAbbAllOrd,
                                            ncol(geneticVarMatAllOrd)),
                             Population = rep(colnames(geneticVarMatAllOrd),
                                              each = nrow(geneticVarMatAllOrd)),
                             Value = c(geneticVarMatAllOrd))
  geneticVarDf$Strategy <- factor(geneticVarDf$Strategy,
                                  levels = allocateStrategyNamesAbbAllOrd)
  # geneticVarDf$Population <- factor(geneticVarDf$Population,
  #                                   labels = paste0("Population-", 1:5))
  geneticVarDf$Population <- factor(geneticVarDf$Population,
                                    labels = paste0(0:4))

  pltChangeVar <- ggplot2::ggplot(data = geneticVarDf) +
    ggplot2::geom_line(mapping = aes(x = Population, y = Value,
                                     group = Strategy,
                                     color = Strategy),
                       lty = rep(allocateStrategyLtys, each = 5),
                       lwd = allocateStrategyLwds) +
    ggplot2::scale_color_manual(values = allocateStrategyCols) +
    ggplot2::theme(text = element_text(size = 96)) +
    ggplot2::xlab(label = "Generation number") +
    ggplot2::ylab(label = "Genetic variance")

  fileNamePltChangeVar <- paste0(dirMidSCOBSTrainingPop1, scriptID, "_",
                                 "change_in_genetic_variance_over_four_generations_for_7_strategies_fp.png")
  # png(filename = fileNamePltChangeVar, width = 1800, height = 1100)
  png(filename = fileNamePltChangeVar, width = 1600, height = 1600)
  print(pltChangeVar)
  dev.off()


  geneticGainIterOptList <- lapply(
    X = geneticGainIterList,
    FUN = function(geneticGainIter) {
      geneticGainIter[(nrow(geneticGainIter) - 1), ]
    }
  )

  geneticGainIterEq <- unlist(lapply(
    X = geneticGainIterList,
    FUN = function(geneticGainIter) {
      geneticGainIter[nrow(geneticGainIter), ]
    }
  ))

  # geneticGainIterEq <- apply(
  #   X = do.call(
  #     what = rbind,
  #     args = lapply(
  #       X = geneticGainIterList,
  #       FUN = function(geneticGainIter) {
  #         geneticGainIter[nrow(geneticGainIter), ]
  #       }
  #     )),
  #   MARGIN = 2,
  #   FUN = mean
  # )




  #### 3.2.3. CDF for final genetic gain ####
  geneticGainIterAllList <- c(
    geneticGainIterOptList,
    list(EQ = geneticGainIterEq)
  )
  geneticGainIterAllOrdList <- geneticGainIterAllList[allocateStrategyNamesAbbAllOrd]

  geneticGainIterDf <- data.frame(
    Strategy = rep(allocateStrategyNamesAbbAllOrd,
                   lapply(X = geneticGainIterAllOrdList,
                          FUN = length)),
    Value = unlist(geneticGainIterAllOrdList)
  )
  geneticGainIterDf$Strategy <- factor(geneticGainIterDf$Strategy,
                                       levels = allocateStrategyNamesAbbAllOrd)
  pltCdf <- ggplot2::ggplot(data = geneticGainIterDf,
                            aes(x = Value, group = Strategy, color = Strategy)) +
    ggplot2::stat_ecdf(pad = TRUE,
                       lty = rep(allocateStrategyLtys,
                                 c(7414, 1240, 2433, 1618, 808, 1213, 1097)),
                       lwd = allocateStrategyLwds) +
    ggplot2::scale_color_manual(values = allocateStrategyCols) +
    ggplot2::theme(text = element_text(size = 96)) +
    ggplot2::xlab(label = "Genetic gain") +
    ggplot2::ylab(label = "CDF")

  fileNamePltCdf <- paste0(dirMidSCOBSTrainingPop1, scriptID, "_",
                           "density_plot_of_final_genetic_gain_for_7_strategies_fp.png")
  # png(filename = fileNamePltCdf, width = 1800, height = 1100)
  png(filename = fileNamePltCdf, width = 2400, height = 1000)
  print(pltCdf)
  dev.off()


  ##### 3.3. Results for ten traits #####
} else {
  names(geneticGainFinalTenList) <- allocateStrategyNamesAbb

  if (overWriteRes | (!fileExists)) {
    saveRDS(object = geneticGainFinalTenList,
            file = fileNameGeneticGainFinalTenList)
  }

  #### 3.3.1. Boxplot for ten traits ####
  geneticGainFinalTenOptList <- lapply(
    X = geneticGainFinalTenList,
    FUN = function(geneticGainFinalTenStr) {
      do.call(
        what = rbind,
        args = lapply(
          X = geneticGainFinalTenStr,
          FUN = function(geneticGainFinalTenPheno) {
            geneticGainFinalTenPheno[(length(geneticGainFinalTenPheno) - 1):length(geneticGainFinalTenPheno)]
          }
        )
      )
    }
  )
  names(geneticGainFinalTenOptList) <- allocateStrategyNamesAbb



  geneticGainFinalTenOptModiList <- c(
    lapply(
      X = geneticGainFinalTenOptList,
      FUN = function(geneticGainFinalTenOpt) {
        geneticGainFinalTenOpt[,1]
      }
    ),
    list(
      EQ = do.call(
        what = c,
        args = lapply(
          X = geneticGainFinalTenOptList,
          FUN = function(geneticGainFinalTenOpt) {
            geneticGainFinalTenOpt[, 2]
          }
        )
      )
    )
  )

  geneticGainFinalTenOptDf <- data.frame(Strategy = rep(names(geneticGainFinalTenOptModiList),
                                                        lapply(X = geneticGainFinalTenOptModiList, FUN = length)),
                                         Value = unlist(geneticGainFinalTenOptModiList))
  geneticGainFinalTenOptDf$Strategy <- factor(x = geneticGainFinalTenOptDf$Strategy,
                                              levels = allocateStrategyNamesAbbAllOrd)

  pltBox <- ggplot2::ggplot(data = geneticGainFinalTenOptDf) +
    ggplot2::geom_boxplot(mapping = aes(x = Strategy, y = Value, fill = Strategy)) +
    ggplot2::scale_fill_manual(values = allocateStrategyCols) +
    ggplot2::theme(text = element_text(size = 72)) +
    ggplot2::ylab(label = "Genetic gain")

  fileNamePltBox <- paste0(dirMidSCOBSGeno1, scriptID, "_",
                           "comparison_of_final_genetic_gain_between_7_strategies_boxplot_fp.png")
  # png(filename = fileNamePltBox, width = 1800, height = 1100)
  png(filename = fileNamePltBox, width = 2400, height = 1000)
  print(pltBox)
  dev.off()


  pltViolin <- ggplot2::ggplot(data = geneticGainFinalTenOptDf) +
    ggplot2::geom_violin(mapping = aes(x = Strategy, y = Value, fill = Strategy)) +
    ggplot2::scale_fill_manual(values = allocateStrategyCols) +
    ggplot2::theme(text = element_text(size = 72)) +
    ggplot2::ylab(label = "Genetic gain")

  fileNamePltViolin <- paste0(dirMidSCOBSGeno1, scriptID, "_",
                              "comparison_of_final_genetic_gain_between_7_strategies_violin_fp.png")
  # png(filename = fileNamePltViolin, width = 1800, height = 1100)
  png(filename = fileNamePltViolin, width = 2400, height = 1000)
  print(pltViolin)
  dev.off()




  geneticGainFinalTenOptModiList$TBVGVP / tapply(X = geneticGainFinalTenOptModiList$EQ, INDEX = rep(1:10, 6), mean)
  geneticGainFinalTenOptEqDiff <- do.call(what = cbind,
                                          args = geneticGainFinalTenOptModiList[1:6]) -
    matrix(rep(tapply(X = geneticGainFinalTenOptModiList$EQ, INDEX = rep(1:10, 6), mean), 6), ncol = 6)
  boxplot(geneticGainFinalTenOptEqDiff)

  geneticGainFinalTenOptEqImp <- round((do.call(what = cbind,
                                                args = geneticGainFinalTenOptModiList[1:6]) /
                                          matrix(rep(tapply(X = geneticGainFinalTenOptModiList$EQ, INDEX = rep(1:10, 6), mean), 6), ncol = 6) - 1) * 100, 3)
  geneticGainFinalTenOptEqImpOrd <- geneticGainFinalTenOptEqImp[, allocateStrategyNamesAbbAllOrd[-1]]


  fileNameGeneticGainFinalTenOptEqImpOrd <- paste0(dirMidSCOBSGeno, scriptID, "_",
                                                   "improvement_from_equal_allocation_for_6_strategies_fp.csv")
  fwrite(x = data.frame(geneticGainFinalTenOptEqImpOrd),
         file = fileNameGeneticGainFinalTenOptEqImpOrd)
}

