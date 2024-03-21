library(LiPAnalyzeR)
library(gginnards)
library(ggplot2)
library(patchwork)

# Loading data
LiPMS <- read.csv("../Documents/LiPAnalyzer/Data/aSyn/OriginalData/2019-04-29_Normalization_LiP_MSStats_Report.csv")
LiPMS$R.FileName <- paste0(LiPMS$R.Condition, "_", LiPMS$R.Replicate)
TrPMS <- read.csv("../Documents/LiPAnalyzer/Data/aSyn/OriginalData/2019-04-29_Normalization_Trp_MSStats_Report.csv")
TrPMS$R.FileName <- paste0(TrPMS$R.Condition, "_", TrPMS$R.Replicate)

annotPepProt <- readRDS("../Documents/LiPAnalyzer/Data/aSyn/PepProtInfo/PepProtInfo_131221.rds")
row.names(annotPepProt) <- annotPepProt$Pep
colnames(annotPepProt) <- c("Peptide", "Protein", "startPosition", "endPosition", "isTryptic")
annotPepProt$quantID <- row.names(annotPepProt)
annotPepProt <- annotPepProt[, c(6,1:5)]

annotSample <- readRDS("../Documents/LiPAnalyzer/Data/aSyn/SampleInfo/Info_221121.rds")
annotSampleAsynF <- annotSample[annotSample$aSyn == "F", ]
annotSampleAsynM <- annotSample[annotSample$aSyn == "M", ]

annotSampleAsyn5 <- annotSample[annotSample$Quantity == 5, ]
annotSampleAsyn20 <- annotSample[annotSample$Quantity == 20, ]

# Extracting peptide and protein quantities 
## Extracting LiPPep, TrpPep and TrpProt into a list
QuantityListRaw <- extractSpectroData(spectroLiP = LiPMS,
                                      spectroTrP = TrPMS)

## Filtering peptide and protein lists, removing peptides with NAs, keeping only FT peptides, ?PreprocessQuantityMatrix for details
QuantityList <- preprocessQuantityMatrix(quantityList = QuantityListRaw,
                                         annotPP = annotPepProt,
                                         annotS = annotSample,  
                                         nameFT = "full")

QuantityListFTHT <- preprocessQuantityMatrix(quantityList = QuantityListRaw,
                                             annotPP = annotPepProt,
                                             annotS = annotSample,
                                             mode = "FTHTjoin", )

# LiPAnalyzeR
## Wrapper function for extracting results
sumModelRes_Quantity <- function(mData){
    summarizeModelResults(resModel = mData, 
                          evalCovariable = "Quantity_Contrast", 
                          correctPval = "all")
}
sumModelRes_aSyn <- function(mData){
    summarizeModelResults(resModel = mData, 
                          evalCovariable = "aSyn_Contrast", 
                          correctPval = "all")
}
## Running LiPAnalyzeR (regressing TrPPep and TrPProt) on complete data set
### Contrast model quantity effect
modelAllQuant <- runModel(quantityList = QuantityList, 
                          annotS = annotSample,
                          formulaContrast = "Y ~ Quantity")
resModelAllQuant <- sumModelRes_Quantity(modelAllQuant)

### Contrast model aSyn effect
modelAllaSyn <- runModel(quantityList = QuantityList, 
                         annotS = annotSample, 
                         formulaContrast = "Y ~ aSyn")
resModelAllaSyn <- sumModelRes_aSyn(modelAllaSyn)

## Running LiPAnalyzeR (regressing TrPPep and TrPProt) on sub-data sets
### Contrast model quantity effect
modelFQuant <- runModel(quantityList = QuantityList, 
                        annotS = annotSampleAsynF, 
                        formulaContrast = "Y ~ Quantity")

resModelFQuant <- sumModelRes_Quantity(modelFQuant)

modelMQuant <- runModel(quantityList = QuantityList, 
                        annotS = annotSampleAsynM, 
                        formulaContrast = "Y ~ Quantity")
resModelMQuant <- sumModelRes_Quantity(modelMQuant)

### Contrast model aSyn effect
model5aSyn <- runModel(quantityList = QuantityList, 
                       annotS = annotSampleAsyn5, 
                       formulaContrast = "Y ~ aSyn")
resModel5aSyn <- sumModelRes_aSyn(model5aSyn)

model20aSyn <- runModel(quantityList = QuantityList, 
                        annotS = annotSampleAsyn20, 
                        formulaContrast = "Y ~ aSyn")
resModel20aSyn <- sumModelRes_aSyn(model20aSyn)


## Running LiPAnalyzeR (regressing only TrPProt) from HT and FT LiPPeps on complete data set
### Contrast model quantity effect
modelAllQuantFTHT <- runModel(quantityList = QuantityListFTHT, 
                             annotS = annotSample,
                             formulaContrast = "Y ~ Quantity", 
                             formulaRUV = "Y ~ XProt", 
                             lowRUV = c(-Inf, 0),
                             upRUV = c(Inf, Inf))
resModelAllQuantFTHT <- sumModelRes_Quantity(modelAllQuantFTHT)

### Contrast model aSyn effect
modelAllaSynFTHT <- runModel(quantityList = QuantityListFTHT, 
                            annotS = annotSample, 
                            formulaContrast = "Y ~ aSyn", 
                            formulaRUV = "Y ~ XProt", 
                            lowRUV = c(-Inf, 0),
                            upRUV = c(Inf, Inf))
resModelAllaSynFTHT <- sumModelRes_aSyn(modelAllaSynFTHT)

## Running LiPAnalyzeR (regressing only TrPProt) from HT and FT LiPPeps on sub-data set
### Contrast model quantity effect
modelFQuantFTHT <- runModel(quantityList = QuantityListFTHT, 
                            annotS = annotSampleAsynF, 
                            formulaContrast = "Y ~ Quantity", 
                            formulaRUV = "Y ~ XProt", 
                            lowRUV = c(-Inf, 0),
                            upRUV = c(Inf, Inf))
resModelFQuantFTHT <- sumModelRes_Quantity(modelFQuantFTHT)

modelMQuantFTHT <- runModel(quantityList = QuantityListFTHT, 
                            annotS = annotSampleAsynM, 
                            formulaContrast = "Y ~ Quantity", 
                            formulaRUV = "Y ~ XProt", 
                            lowRUV = c(-Inf, 0),
                            upRUV = c(Inf, Inf))
resModelMQuantFTHT <- sumModelRes_Quantity(modelMQuantFTHT)

### Contrast model aSyn effect
model5aSynFTHT <- runModel(quantityList = QuantityListFTHT, 
                           annotS = annotSampleAsyn5, 
                           formulaContrast = "Y ~ aSyn", 
                           formulaRUV = "Y ~ XProt", 
                           lowRUV = c(-Inf, 0),
                           upRUV = c(Inf, Inf))
resModel5aSynFTHT <- sumModelRes_aSyn(model5aSynFTHT)

model20aSynFTHT <- runModel(quantityList = QuantityListFTHT, 
                            annotS = annotSampleAsyn20, 
                            formulaContrast = "Y ~ aSyn", 
                            formulaRUV = "Y ~ XProt", 
                            lowRUV = c(-Inf, 0),
                            upRUV = c(Inf, Inf))
resModel20aSynFTHT <- sumModelRes_aSyn(model20aSynFTHT)

### Plotting results in woods plots
getPlot <- function(resModel, title){
    p <- makeWoodsPlotSingleProtein(sumDf = resModel,
                                    annotPP = annotPepProt, 
                                    protName = "P37840.1", 
                                    nameFT = "full", 
                                    nameHT = "half", 
                                    ylim = c(-5, 2.7),
                                    deltaColorIsTryptic = T) + ggtitle(title) + theme_bw (base_size = 11) + theme(legend.position = "none")
    ## Add coloring NAC region
    append_layers(p, geom_rect(aes(xmin = 61, xmax = 95, 
                                   ymin = -Inf, ymax = Inf), 
                               alpha = .3, color = "#eddddd", 
                               fill = "#eddddd"), position = "bottom")
}

#### Extended Data Figure 7b
getPlot(resModelAllaSynFTHT, "aSyn effect in all samples")
#### Extended Data Figure 7c
getPlot(resModel5aSynFTHT, "aSyn effect in F1 & M1")
#### Extended Data Figure 7d
getPlot(resModel20aSynFTHT, "aSyn effect in F2 & M2")

#### Extended Data Figure 7e
getPlot(resModelAllQuantFTHT, "Concentration effect in all samples")
#### Extended Data Figure 7f
getPlot(resModelMQuantFTHT, "Concentration effect in M1 & M2")
#### Extended Data Figure 7g
getPlot(resModelFQuantFTHT, "Concentration effect in F1 & F2")

