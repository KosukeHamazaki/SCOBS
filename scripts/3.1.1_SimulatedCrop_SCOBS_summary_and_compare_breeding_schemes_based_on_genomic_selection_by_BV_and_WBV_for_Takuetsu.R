#############################################################################################################################################
######  Title: 3.1.1_SimulatedCrop_SCOBS_summary_and_compare_breeding_schemes_based_on_genomic_selection_by_BV_and_WBV_for_Takuetsu    ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                                                                               ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo                                                      ######
######  Date: 2020/02/25 (Created), 2021/02/26 (Last Updated)                                                                          ######
#############################################################################################################################################





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
  if (os == "macOS  10.16") {
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
scriptIDMain <- "2.1"
scriptIDSummary <- "3.1"
scriptID <- "3.1.1"


##### 1.2. Setting some parameters #####
dirMidSCOBSBase <- "midstream/"


#### 1.2.1. Setting some parameters related to names of R6 class ####
scenarioNo <- 1
simName <- paste0("Scenario_", scenarioNo)
# breederName <- "K. Hamazaki"

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
# nSimPheno <- 10
#
# nIterSimulation <- 100
# nGenerationProceed <- 9
# nRefreshMemoryEvery <- 2
# updateBreederInfo <- rep(TRUE, nGenerationProceed)
# phenotypingInds <- rep(FALSE, nGenerationProceed)
# nRepForPhenoInit <- 3
# nRepForPheno <- rep(1, nGenerationProceed)
# updateModels <- rep(FALSE, nGenerationProceed)

simBsNameBase <- "VerySimpleGS-for-Trait_1"
# strategyNames <- c(
#   "selectBV", "selectWBV",
#   "selectOHV", "selectOPV"
#   # , "selectEMBV"
# )
#
# nTrainingPopInits <- list("true", 2000,  1500,
#                           1000, 750, 500, 250)
#
# lociEffMethods <- c("true", "estimated")
# trainingPopType <- "latest"
# trainingPopInit <- 1
# trainingIndNamesInit <- NULL


#### 1.2.3. Setting some parameters related to estimation of marker effects ####
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


#### 1.2.4. Setting some parameters related to selection of parent candidates ####
# nSelectionWaysVec <- rep(1, nGenerationProceed)
# traitNoSelList <- rep(list(list(1)), nGenerationProceed)
#
# blockSplitMethod <- "minimumSegmentLength"
# nMrkInBlock <- 10
# minimumSegmentLength <- 20
# nSelInitOPVList <- rep(list(50), nGenerationProceed)
# nIterOPV <- 5e04
# nProgeniesEMBVVec <- rep(40, nGenerationProceed)
# nIterEMBV <- 3
# nCoresEMBV <- 1
#
# clusteringForSelList <- rep(list(FALSE), nGenerationProceed)
# nClusterList <- rep(list(1), nGenerationProceed)
# nTopClusterList <- rep(list(1), nGenerationProceed)
# nTopEachList <- rep(list(15), nGenerationProceed)
# nSelList <- rep(list(15), nGenerationProceed)
#
# multiTraitsEvalMethodList <- rep(list("sum"), nGenerationProceed)
# hSelList <- rep(list(list(1)), nGenerationProceed)


#### 1.2.5. Setting some parameters related to mating and resource allocation ####
# matingMethodVec <- rep("diallelWithSelfing", nGenerationProceed)
# allocateMethodVec <- rep("equalAllocation", nGenerationProceed)
# weightedAllocationMethodList <- rep(list("selectBV"), nGenerationProceed)
# traitNoRAList <- rep(list(1), nGenerationProceed)
# hList <- rep(list(1), nGenerationProceed)
# includeGVPVec <- rep(FALSE, nGenerationProceed)
# nNextPopVec <- rep(250, nGenerationProceed)


#### 1.2.6. Setting some parameters related to mating and resource allocation ####
# nCores <-  100
# nCoresWhenOPV <-  100
# nCoresWhenEMBV <- 34
# nCoresSummary <- 10
#
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


#### 1.2.7. Save parameters ####
fileParamsSCOBS <- paste0(dirMidSCOBSTrial, scriptIDMain,
                          "_", project, "_", trialName,
                          "_", simBsNameBase,
                          "_all_parameters.RData")
load(fileParamsSCOBS)

strategyNames <- c(
  "selectBV", "selectWBV",
  "selectOHV", "selectOPV"
  , "selectEMBV"
)
scriptIDData <- "0.1"
scriptIDMain <- "2.1"
scriptIDSummary <- "3.1"
scriptID <- "3.1.1"

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
trainingPopInitNo <- 1
for (simGenoNo in 1:nSimGeno) {
  genoName <- paste0("Geno_", simGenoNo)

  dirMidSCOBSGeno <- paste0(dirMidSCOBSTrial, scriptIDData,
                            "_", genoName, "/")

  trueGVSummaryMeanArrayPheno <- array(data = NA,
                                       dim = c(5, 10, 14, nSimPheno))
  for (simPhenoNo in 1:nSimPheno) {
    phenoName <- paste0("Pheno_", simPhenoNo)
    print(phenoName)

    dirMidSCOBSPheno <- paste0(dirMidSCOBSGeno, scriptIDData,
                               "_", phenoName, "/")
    simBsNamePheno <- paste0(simBsNameBase, "_",
                             genoName, "_",
                             phenoName)

    fileNameSimEval <- paste0(dirMidSCOBSPheno, scriptIDSummary, "_", simBsNamePheno, "_simEval.rds")
    mySimEval <- readRDS(file = fileNameSimEval)
    print("Loaded")

    twoMethodsNos <- sort(unlist(x = sapply(X = c("selectBV", "selectWBV"),
                                            FUN = function (x) {
                                              grep(pattern = x,
                                                   x = names(mySimEval$simBsList))
                                            }, simplify = FALSE)))


    fileNameSimEvalTwo <- paste0(dirMidSCOBSPheno, scriptID, "_", simBsNamePheno,
                                 "-BV_WBV", "_simEval.rds")
    mySimEvalTwo <- myBreedSimulatR::simEval$new(simEvalName = paste0("simBsNamePheno", "-BV_WBV"),
                                                 simBsList = mySimEval$simBsList[twoMethodsNos],
                                                 verbose = TRUE)

    saveRDS(object = mySimEvalTwo, file = fileNameSimEvalTwo)
    print("Saved two methods")

    pltTwoOrigin <- mySimEvalTwo$plot(targetTrait = 1,
                                      targetPopulation = 1:10,
                                      plotType = "lines",
                                      plotTarget = "max",
                                      returnGain = FALSE,
                                      plotGVMethod = "true")
    # simBsEach <- mySimEvalTwo$simBsList[[1]]
    pltTwoOrigin <- pltTwoOrigin %>% plotly::layout(showlegend = TRUE,
                                                    xaxis = list(tickvals = as.list(0:9),
                                                                 ticktext = as.list(paste0("Population_", 1:10)),
                                                                 tickmode = "array"))
    pltTwo <- pltTwoOrigin %>% plotly::layout(showlegend = FALSE)


    fileNamePltTwoMaxPng <- paste0(dirMidSCOBSPheno, scriptID, "_", simBsNamePheno,
                                   "-BV_WBV", "_plotly_lines_max.png")
    fileNamePltTwoMaxHtml <- paste0(dirMidSCOBSPheno, scriptID, "_", simBsNamePheno,
                                    "-BV_WBV", "_plotly_lines_max.html")
    plotly::orca(p = pltTwo,
                 format = "png",
                 width = 1200,
                 height = 800,
                 file = fileNamePltTwoMaxPng)

    htmlwidgets::saveWidget(widget = plotly::partial_bundle(pltTwoOrigin),
                            file = here::here(fileNamePltTwoMaxHtml))

    print("Saved two images")

    trueGVSummaryMeanArray <- do.call(
      what = abind::abind,
      args = lapply(X = mySimEvalTwo$simBsList,
                    FUN = function (simBsEach) {
                      trueGVSummaryArrayMean <-
                        apply(X = simBsEach$trueGVSummaryArray,
                              MARGIN = c(1, 3, 2),
                              FUN = mean)
                      dimnames(trueGVSummaryArrayMean)[[3]] <- simBsEach$simBsName

                      return(trueGVSummaryArrayMean)
                    })
    )

    trueGVSummaryMeanArrayPheno[, , , simPhenoNo] <- trueGVSummaryMeanArray
  }

  dimnames(trueGVSummaryMeanArrayPheno) <-
    c(dimnames(trueGVSummaryMeanArray)[1],
      list(paste0("Population_", 1:10)),
      dimnames(trueGVSummaryMeanArray)[3],
      list(paste0("Pheno_", 1:nSimPheno)))


  fileTrueGVSummaryMeanArrayPheno <- paste0(dirMidSCOBSGeno, scriptID,
                                            "_",
                                            simBsNameBase,
                                            "_summary_of_true_GV_array.rds")
  saveRDS(object = trueGVSummaryMeanArrayPheno,
          file = fileTrueGVSummaryMeanArrayPheno)

  trueGVMaxMeanFinal <- t(trueGVSummaryMeanArrayPheno[1, 10, , ])

  trueGVMaxMeanFinalDf <- data.frame(Method = factor(rep(rep(c("BV", "WBV"), ncol(trueGVMaxMeanFinal) / 2),
                                                         each = nrow(trueGVMaxMeanFinal)),
                                                     levels = c("BV", "WBV")),
                                     Train = factor(rep(unlist(nTrainingPopInits),
                                                        each = 2 * nrow(trueGVMaxMeanFinal)),
                                                    levels = unlist(nTrainingPopInits)),
                                     Pheno = rep(rownames(trueGVMaxMeanFinal),
                                                 ncol(trueGVMaxMeanFinal)),
                                     Value = c(trueGVMaxMeanFinal))

  pltBoxMaxFinal <- ggplot2::ggplot(data = trueGVMaxMeanFinalDf) +
    ggplot2::geom_boxplot(mapping = aes(x = Train, y = Value, fill = Method)) +
    ggplot2::theme(text = element_text(size = 28))

  filePltBoxMaxFinal <- paste0(dirMidSCOBSGeno, scriptID,
                               "_",
                               simBsNameBase,
                               "_true_GV_max5_in_final_population_BV_WBV_boxplot.png")
  png(filename = filePltBoxMaxFinal, width = 1100, height = 900)
  print(pltBoxMaxFinal)
  dev.off()

}
