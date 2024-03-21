library(limma)
library(LiPAnalyzeR)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(patchwork)

# Loading data
LiPSpectro <- read.csv("../Documents/LiPAnalyzer/Data/PombeStrains/RawSpectroFiles/Export_210422_LiPAnalyzerSchema_AddMS2RawQ/LiP_SpectroOut_LipAnalyzerScheme_210422.csv")
TrPSpectro <- read.csv("../Documents/LiPAnalyzer/Data/PombeStrains/RawSpectroFiles/Export_210422_LiPAnalyzerSchema_AddMS2RawQ/Trp_SpectroOut_LipAnalyzerScheme_210422.csv")
annotSample <- readRDS("/cellfile/cellnet/Luise/LiPData/PombeStrains_LiP/SampleInfo/README_SampleInfo_250422.RDS")
annotSample <- annotSample[order(annotSample$Strain),]
annotSample <- annotSample[order(annotSample$Replicate),]
row.names(annotSample) <- paste0(annotSample$Strain, "_", annotSample$Replicate)

# Extracting peptide and protein quantities 
## Extracting LiPPep, TrPPep and TrPProt into a list
QuantityListRaw <- extractSpectroData(spectroLiP = LiPSpectro, 
                                      spectroTrP = TrPSpectro)

QuantityListLiPRaw <- extractSpectroData(spectroLiP = LiPSpectro,
                                         LiPonly = T)

### Adjusting colnames of data
colnames(QuantityListRaw[[1]]) <- unlist(lapply(strsplit(colnames(QuantityListRaw[[1]]), "_"), function(x) x[[4]]))
colnames(QuantityListRaw[[1]]) <- row.names(annotSample)[match(colnames(QuantityListRaw[[1]]), annotSample$RunNumber_LiP)]

colnames(QuantityListRaw[[2]]) <- unlist(lapply(strsplit(colnames(QuantityListRaw[[2]]), "_"), function(x) x[[4]]))
colnames(QuantityListRaw[[2]]) <- row.names(annotSample)[match(colnames(QuantityListRaw[[2]]), annotSample$RunNumber_Trp)]

colnames(QuantityListRaw[[3]]) <- unlist(lapply(strsplit(colnames(QuantityListRaw[[3]]), "_"), function(x) x[[4]]))
colnames(QuantityListRaw[[3]]) <- row.names(annotSample)[match(colnames(QuantityListRaw[[3]]), annotSample$RunNumber_Trp)]

colnames(QuantityListLiPRaw[[1]]) <- unlist(lapply(strsplit(colnames(QuantityListLiPRaw[[1]]), "_"), function(x) x[[4]]))
colnames(QuantityListLiPRaw[[1]]) <- row.names(annotSample)[match(colnames(QuantityListLiPRaw[[1]]), annotSample$RunNumber_LiP)]

colnames(QuantityListLiPRaw[[2]]) <- unlist(lapply(strsplit(colnames(QuantityListLiPRaw[[2]]), "_"), function(x) x[[4]]))
colnames(QuantityListLiPRaw[[2]]) <- row.names(annotSample)[match(colnames(QuantityListLiPRaw[[2]]), annotSample$RunNumber_LiP)]

# Further sorting of the data
## Creating data.frame with information to all peptides in data
annotPepProt <- getPepProtAnnotSpectro(spectroOut = TrPSpectro, 
                                       spectroOut2 = LiPSpectro)

## Sorting cols
QuantityListRaw <- lapply(QuantityListRaw, \(x){x[, row.names(annotSample)]})
QuantityListLiPRaw <- lapply(QuantityListLiPRaw, \(x){x[, row.names(annotSample)]})


## Filtering peptide and protein lists, removing peptides with NAs
QuantityList <- preprocessQuantityMatrix(quantityList = QuantityListRaw,
                                         annotPP = annotPepProt)

QuantityListFTHT <- preprocessQuantityMatrix(quantityList = QuantityListRaw,
                                             annotPP = annotPepProt,
                                             mode = "FTHTjoin")


QuantityListHTonly <- preprocessQuantityMatrix(quantityList = QuantityListRaw,
                                               annotPP = annotPepProt,
                                               mode = "HTonly")

QuantityListHTonly <- lapply(QuantityListHTonly, \(x){
    x <- x[annotPepProt[row.names(x), "isTryptic"] == "Specific-C"|
               annotPepProt[row.names(x), "isTryptic"] == "Specific-N",]
    return(x)
})  # keeping only HT peptides, not those matching different digestion sites

QuantityListLiPonly <- preprocessQuantityMatrix(quantityList = QuantityListLiPRaw,
                                                annotPP = annotPepProt,
                                                mode = "LiPonly", 
                                                filterTryptic = F)

## Batch correct data
QuantityListBC <- lapply(QuantityList, \(x){
    removeBatchEffect(x, annotSample$Batch)
})

QuantityListFTHTBC <- lapply(QuantityListFTHT, \(x){
    removeBatchEffect(x, annotSample$Batch)
})

QuantityListLiPonlyBC <- lapply(QuantityListLiPonly, \(x){
    removeBatchEffect(x, annotSample$Batch)
})

QuantityListHtonlyBC <- lapply(QuantityListHTonly, \(x){
    removeBatchEffect(x, annotSample$Batch)
})




# Overview figures of data
## Correlations
CorCols <- function(Peps, df1, df2,use = "pairwise.complete.obs"){
    j <- intersect(colnames(df1), colnames(df2))
    sapply(Peps, function(i){
        cor(as.numeric(df1[i,j]), as.numeric(df2[i,j]), use = use)
    })
}

## Plotting correlations of LiP peptides
plotDensity <- function(plotData){
    ggplot(plotData, aes(x = value, color = variable, fill = variable, alpha = variable)) +
        geom_density(linewidth = 0.7) +
        ylab("Density") +
        xlab("Correlation coefficient") +
        xlim(-1, 1) +
        geom_vline(xintercept = median(plotData$TrPProt), color = "#0057b7", linewidth = 1) +
        geom_vline(xintercept = median(plotData$TrPPep), color = "#ffc900", linewidth = 1) +
        scale_color_manual(values = c("#0057b7", "#ffc900")) +
        scale_fill_manual(values = c("#0057b7", "#ffc900")) +
        scale_alpha_manual(values = c(0.8, 0.4)) +
        theme_bw(base_size = 11)
}
plotScatter <- function(plotData){
    ggplot(plotData, aes(x = TrPProt, y = TrPPep)) +
        geom_point(alpha = 0.1, size = 2) +
        geom_abline(intercept = 0, color = "lightblue", linewidth = 1) +
        xlab("Correlations LiPPep & TrPProt") +
        ylab("Correlations LiPPep & TrPPep") +
        xlim(-1, 1) +
        ylim(-1, 1) +
        theme_bw(base_size = 11)
}
   
### with FT peps
plotData <- data.frame(TrPProt = CorCols(row.names(QuantityListBC$LiPPep), 
                                         QuantityListBC$LiPPep, 
                                         QuantityListBC$TrPProt),
                       TrPPep = CorCols(row.names(QuantityListBC$LiPPep), 
                                        QuantityListBC$LiPPep, 
                                        QuantityListBC$TrPPep))
#### Extended Data Figure 1b
plotScatter(plotData)

plotData <- melt(plotData)
#### Extended Data Figure 1a
plotDensity(plotData)

### HT peps with matched FT peps
plotData <- data.frame(TrPProt = CorCols(row.names(QuantityListHtonlyBC$LiPPep),
                                         QuantityListHtonlyBC$LiPPep,
                                         QuantityListHtonlyBC$TrPProt),
                       TrPPep = CorCols(row.names(QuantityListHtonlyBC$LiPPep),
                                        QuantityListHtonlyBC$LiPPep,
                                        QuantityListHtonlyBC$TrPPep))
plotData <- melt(plotData)

#### Extended Data Figure 5a
plotDensity(plotData)




## Plotting single example peptides
### Example of protein abundance effect
#### Figure 3a
i <- "AAPSPIAYVAEIR"

### Example of PK-independent peptide effect
#### Figure 3b & Figure 8d
i <- "VSALSGFEGDATPFTDVTVEAVSK"

### Example of PK-dependent peptide effect
#### Figure 3c
i <- "DTFIILR"

### Example of protein abundance effect & PK-independent peptide effect resulting in implausible coefficients if RUV is run without constraints
#### Figure 6c
i <- "RALIDSPCSEFPR"

### Example of protein abundance effect resulting in implausible coefficients if RUV and contrast model are run in a combined model
#### Figure 6d
i <- "GLPLEAVTTIAK"

plotData <- QuantityListBC

### TrPPep vs TrPProt
ggplot(mapping = aes(x = as.numeric(plotData$TrPProt[i, ]),
                     y = as.numeric(plotData$TrPPep[i, ]))) +
    geom_point(size = 3, mapping = aes(color = annotSample[colnames(plotData$LiPPep), "Strain"],
                                       alpha = annotSample[colnames(plotData$LiPPep), "Strain"])) +
    scale_color_manual(name = "Strain", values = c("#440154", "#25848e", "#bbdf27", "#440154")) +
    scale_alpha_manual(name = "Strain", values = c(1, 1, 1, 0.3)) +
    geom_line(stat="smooth",method = "lm", formula = y ~ x,
              linewidth = 1, linetype ="dashed", alpha = 0.5) +
    xlab("TrPProt quantity") +
    ylab("TrPPep quantity") ++
    theme_bw(base_size = 11) +
  theme(aspect.ratio = 1, legend.position = "none")

### LiPPPep vs TrPProt
ggplot(mapping = aes(x = as.numeric(plotData$TrPProt[i, ]),
                     y = as.numeric(plotData$LiPPep[i, ]))) +
    geom_point(size = 3, mapping = aes(color = annotSample[colnames(plotData$LiPPep), "Strain"],
                                       alpha = annotSample[colnames(plotData$LiPPep), "Strain"])) +
    scale_color_manual(name = "Strain", values = c("#440154", "#25848e", "#bbdf27", "#440154")) +
    scale_alpha_manual(name = "Strain", values = c(1, 1, 1, 0.3)) +
    geom_line(stat="smooth",method = "lm", formula = y ~ x,
              linewidth = 1, linetype ="dashed", alpha = 0.5) +
    xlab("TrPProt quantity") +
    ylab("LiPPep quantity") +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1, legend.position = "none")

### LiPPPep vs TrPPep
ggplot(mapping = aes(x = as.numeric(plotData$TrPPep[i, ]),
                     y = as.numeric(plotData$LiPPep[i, ]))) +
    geom_point(size = 3, mapping = aes(color = annotSample[colnames(plotData$LiPPep), "Strain"],
                                       alpha = annotSample[colnames(plotData$LiPPep), "Strain"])) +
    scale_color_manual(name = "Strain", values = c("#440154", "#25848e", "#bbdf27", "#440154")) +
    scale_alpha_manual(name = "Strain", values = c(1, 1, 1, 0.3)) +
    geom_line(stat="smooth",method = "lm", formula = y ~ x,
              linewidth = 1, linetype ="dashed", alpha = 0.5) +
    xlab("TrPPep quantity") +
    ylab("LiPPep quantity") +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1, legend.position = "none")

### LiPPPep vs LiPProt
ggplot(mapping = aes(x = as.numeric(plotData$LiPProt[i, ]),
                     y = as.numeric(plotData$LiPPep[i, ]))) +
    geom_point(size = 3, mapping = aes(color = annotSample[colnames(plotData$LiPPep), "Strain"],
                                       alpha = annotSample[colnames(plotData$LiPPep), "Strain"])) +
    scale_color_manual(name = "Strain", values = c("#440154", "#25848e", "#bbdf27", "#440154")) +
    scale_alpha_manual(name = "Strain", values = c(1, 1, 1, 0.3)) +
    geom_line(stat="smooth",method = "lm", formula = y ~ x,
              linewidth = 1, linetype ="dashed", alpha = 0.5) +
    xlab("LiPProt quantity") +
    ylab("LiPPep quantity") +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1, legend.position = "none")




## Variation of HT and FT peptides
### Creating HT and FT specific LiPPep df
LiPHT <- QuantityListLiPonlyBC$LiPPep[annotPepProt[row.names(QuantityListLiPonlyBC$LiPPep), "isTryptic"] == "Specific-C"|
                                          annotPepProt[row.names(QuantityListLiPonlyBC$LiPPep), "isTryptic"] == "Specific-N", ]
LiPFT <- QuantityListLiPonlyBC$LiPPep[annotPepProt[row.names(QuantityListLiPonlyBC$LiPPep), "isTryptic"] == "Specific", ]

## Getting strain specific HT and FT sd
sd_LiPHT <- data.frame(JB50 = apply(LiPHT[,grepl("JB50", colnames(LiPHT))], 1, sd),
                       JB759 = apply(LiPHT[,grepl("JB759", colnames(LiPHT))], 1, sd),
                       JB760 = apply(LiPHT[,grepl("JB760", colnames(LiPHT))], 1, sd),
                       PYK1_A = apply(LiPHT[,grepl("PYK1_A", colnames(LiPHT))], 1, sd))
sd_LiPHT$mean <- rowMeans(sd_LiPHT)

sd_LiPFT <- data.frame(JB50 = apply(LiPFT[,grepl("JB50", colnames(LiPFT))], 1, sd),
                       JB759 = apply(LiPFT[,grepl("JB759", colnames(LiPFT))], 1, sd),
                       JB760 = apply(LiPFT[,grepl("JB760", colnames(LiPFT))], 1, sd),
                       PYK1_A = apply(LiPFT[,grepl("PYK1_A", colnames(LiPFT))], 1, sd))
sd_LiPFT$mean <- rowMeans(sd_LiPFT)

#### Extended Data Figure 5g
ggplot(mapping = aes(y = c(sd_LiPFT$mean, sd_LiPHT$mean),
                     x = c(rep("SD LiPPep FT", nrow(sd_LiPFT)),
                           rep("SD LiPPep HT", nrow(sd_LiPHT))))) +
    geom_violin(alpha = 0.5, fill = "#a7a7a7") +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.3,
                 colour = "black") +
    xlab("") +
    ylab("Mean standard deviation per LiP peptide") +
    theme_bw(base_size = 15)




# LiP/TrP ratio correction approach
## Estimating ratios
rLiPPepTrPPep <- QuantityListBC$LiPPep- QuantityListBC$TrPPep
rLiPPepTrPProt <- QuantityListBC$LiPPep- QuantityListBC$TrPProt

## Running contrast model of LiPAnalyzeR on ratios 
CON_LiPPepTrPPepRatio <- runModel(list(Y = rLiPPepTrPPep),
                                  annotS = annotSample,
                                  formulaRUV = NULL, 
                                  formulaContrast = "Y~Strain")

CON_LiPPepTrPProtRatio <- runModel(list(Y = rLiPPepTrPProt),
                                   annotS = annotSample,
                                   formulaRUV = NULL,
                                   formulaContrast = "Y~Strain")

ResJB759_TrPPepRatio <- summarizeModelResults(resModel = CON_LiPPepTrPPepRatio, 
                                              evalCovariable = "StrainJB759_Contrast",
                                              correctPval  = "protein-wise",
                                              annotPP = annotPepProt)

ResJB759_TrPProtRatio <- summarizeModelResults(resModel = CON_LiPPepTrPProtRatio, 
                                               evalCovariable = "StrainJB759_Contrast",
                                               correctPval  = "protein-wise",
                                               annotPP = annotPepProt)

## Correlations ratio LipPep/TrPPep with TrPPep
cRLiPPepTrPPep_TrPPep <- CorCols(row.names(rLiPPepTrPPep),
                                 rLiPPepTrPPep,
                                 QuantityListBC$TrPPep)

## Correlations ratio LipPep/TrPProt with TrPProt
cRLiPPepTrPProt_TrPProt <- CorCols(row.names(rLiPPepTrPProt),
                                   rLiPPepTrPProt,
                                   QuantityListBC$TrPProt)

### Plotting correlations
plotData <- melt(data.frame(TrPProt = cRLiPPepTrPProt_TrPProt,
                            TrPPep = cRLiPPepTrPPep_TrPPep))

#### Figure 5a
ggplot(plotData, aes(x = value, color = variable, fill = variable, alpha = variable)) +
    geom_density(linewidth = 0.7) +
    ylab("Density") +
    xlab("Correlation coefficient") +
    geom_vline(xintercept = median(cRLiPPepTrPProt_TrPProt), color = "#0057b7", linewidth = 1) +
    geom_vline(xintercept = median(cRLiPPepTrPPep_TrPPep), color = "#ffc900", linewidth = 1) +
    xlim(-1, 1) +
    scale_color_manual(values = c("#0057b7", "#ffc900")) +
    scale_fill_manual(values = c("#0057b7", "#ffc900")) +
    scale_alpha_manual(values = c(0.8, 0.4)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")




# LiPAnalyzeR
## Running only RUV step
### RUV with constraints removing only TrPPep
RUV_TrPPepB <- runModel(quantityList = QuantityList,
                        annotS = annotSample,
                        formulaRUV = "Y~XPep+Batch", 
                        lowRUV = c(-Inf, 0),
                        upRUV = c(Inf, Inf),
                        formulaContrast = NULL,
                        addRUVbounds = T)

### RUV with constraints removing only TrPProt
RUV_TrPProtB <-runModel(quantityList = QuantityList,
                        annotS = annotSample,
                        formulaRUV = "Y~XProt+Batch", 
                        lowRUV = c(-Inf, 0),
                        upRUV = c(Inf, Inf),
                        formulaContrast = NULL,
                        addRUVbounds = T)

### RUV with constraints removing TrPPep & TrPProt
RUV_TrPPepTrPProtB <- runModel(quantityList = QuantityList,
                               annotS = annotSample,
                               formulaRUV = "Y~XPep+XProt+Batch",
                               formulaContrast = NULL,
                               addRUVbounds = T)

### RUV with constraints removing TrPPep & TrPProt in HTonly mode 
RUV_TrPPepTrPProtB_HTonly <- runModel(quantityList = QuantityListHTonly,
                                      annotS = annotSample,
                                      formulaRUV = "Y~XPep+XProt+Batch",
                                      formulaContrast = NULL,
                                      addRUVbounds = T)

### RUV with constraints removing LiPProt in LiPonly mode 
RUV_LiPProtB <- runModel(quantityList = QuantityListLiPonly,
                         annotS = annotSample,
                         formulaRUV = "Y~XProt+Batch", 
                         lowRUV = c(-Inf, 0),
                         upRUV = c(Inf, Inf),
                         formulaContrast = NULL,
                         addRUVbounds = T)

### RUV step without constraints removing only TrPPep
RUVnC_TrPPepB <- runModel(quantityList = QuantityList,
                          annotS = annotSample,
                          formulaRUV = NULL,
                          formulaContrast = "Y~XPep+Batch")

### RUV step without constraints removing only TrPProt
RUVnC_TrPProtB <- runModel(quantityList = QuantityList,
                           annotS = annotSample,
                           formulaRUV = NULL,
                           formulaContrast = "Y~XProt+Batch")

### RUV step without constraints removing TrPPep & TrPProt
RUVnC_TrPPepTrPProtB <- runModel(quantityList = QuantityList,
                                 annotS = annotSample,
                                 formulaRUV = NULL,
                                 formulaContrast = "Y~XPep+XProt+Batch")

## Running complete LiPAnalyzeR models
### default run, RUV with constraints removing TrPPep & TrPProt and CON with Strain
RUV_TrPPepTrPProtB_CON_Strain <- runModel(quantityList = QuantityList,
                                          annotS = annotSample,
                                          formulaRUV = "Y~XPep+XProt+Batch",
                                          formulaContrast = "Y~Strain",
                                          addRUVbounds = T)

### default run, alternating which strain is set as reference level in dummy coding
annotSample_refJB50 <- annotSample
annotSample_refJB50$Strain <- factor(annotSample_refJB50$Strain, levels = c("JB50", "JB759", "JB760", "PYK1_A"))
annotSample_refJB759 <- annotSample
annotSample_refJB759$Strain <- factor(annotSample_refJB759$Strain, levels = c("JB759", "JB50", "JB760", "PYK1_A"))
annotSample_refJB760 <- annotSample
annotSample_refJB760$Strain <- factor(annotSample_refJB760$Strain, levels = c("JB760", "JB50", "JB759", "PYK1_A"))
annotSample_refPYK1_A <- annotSample
annotSample_refPYK1_A$Strain <- factor(annotSample_refPYK1_A$Strain, levels = c("PYK1_A", "JB50", "JB759", "JB760"))

RUV_TrPPepTrPProtB_CON_Strain_refJB50 <- runModel(quantityList = QuantityList,
                                                  annotS = annotSample_refJB50,
                                                  formulaRUV = "Y~XPep+XProt+Batch",
                                                  formulaContrast = "Y~Strain",
                                                  addRUVbounds = T)

RUV_TrPPepTrPProtB_CON_Strain_refJB759 <- runModel(quantityList = QuantityList,
                                                   annotS = annotSample_refJB759,
                                                   formulaRUV = "Y~XPep+XProt+Batch",
                                                   formulaContrast = "Y~Strain",
                                                   addRUVbounds = T)

RUV_TrPPepTrPProtB_CON_Strain_refJB760 <- runModel(quantityList = QuantityList,
                                                   annotS = annotSample_refJB760,
                                                   formulaRUV = "Y~XPep+XProt+Batch",
                                                   formulaContrast = "Y~Strain",
                                                   addRUVbounds = T)

RUV_TrPPepTrPProtB_CON_Strain_refPYK1A <- runModel(quantityList = QuantityList,
                                                   annotS = annotSample_refPYK1_A,
                                                   formulaRUV = "Y~XPep+XProt+Batch",
                                                   formulaContrast = "Y~Strain",
                                                   addRUVbounds = T)

### RUV without constraints removing TrPPep and TrPProt, CON modeling strain
RUVnC_TrPPepTrPProtB_CON_Strain <- runModel(quantityList = QuantityList,
                                            annotS = annotSample,
                                            formulaRUV = "Y~XPep+XProt+Batch",
                                            formulaContrast = "Y~Strain",
                                            lowRUV = c(-Inf, -Inf, -Inf),
                                            upRUV = c(Inf, Inf, Inf),
                                            addRUVbounds = T)

### RUV without constraints, modeling TrPPep, TrPProt and strain in one model
RUVnC_TrPPepTrPProtBStrain <- runModel(quantityList = QuantityList,
                                       annotS = annotSample,
                                       formulaRUV = NULL,
                                       formulaContrast = "Y~XPep+XProt+Batch+Strain")

### RUV LiPonly join, removing LiPProt and modeling strain in CON model
RUV_LiPProtB_CON_Strain_LiPonly <- runModel(quantityList = QuantityListLiPonly,
                                            annotS = annotSample,
                                            formulaRUV = "Y~XProt+Batch",
                                            formulaContrast = "Y~Strain",
                                            lowRUV = c(-Inf, 0),
                                            upRUV = c(Inf, Inf),
                                            addRUVbounds = T)

### RUV FTHT join, removing TrPProt and modeling strain in CON model
RUV_TrPProtB_CON_Strain_FTHTjoin <- runModel(quantityList = QuantityListFTHT,
                                             annotS = annotSample,
                                             formulaRUV = "Y~XProt+Batch",
                                             formulaContrast = "Y~Strain",
                                             lowRUV = c(-Inf, 0),
                                             upRUV = c(Inf, Inf),
                                             addRUVbounds = T)

### RUV HTonly join, removing TrPPep & TrPProt and modeling strain in CON model
RUV_TrPPepTrPProtB_CON_Strain_HTonly <- runModel(quantityList = QuantityListHTonly,
                                                 annotS = annotSample,
                                                 formulaRUV = "Y~XProt+XPep+Batch",
                                                 formulaContrast = "Y~Strain",
                                                 addRUVbounds = T)




## Summarizing model results
### Summarizing results and hits for every strain in one
sumModelsAllStrains <- function(resModel, S1 = "StrainJB759_Contrast", 
                                S2 = "StrainJB760_Contrast", S3 = "StrainPYK1_A_Contrast", 
                                listNames = c("JB759", "JB760","PYK1A"),
                                corPval = "protein-wise", annotPP){
    resList <- list(N1 = summarizeModelResults(resModel = resModel,
                                               evalCovariable = S1,
                                               correctPval  = corPval,
                                               annotPP = annotPP),
                    N2 = summarizeModelResults(resModel = resModel, 
                                               evalCovariable = S2,
                                               correctPval  = corPval,
                                               annotPP = annotPP),
                    N3 = summarizeModelResults(resModel = resModel,
                                               evalCovariable = S3,
                                               correctPval  = corPval,
                                               annotPP = annotPP))
    names(resList) <- listNames
    return(resList)
}

### LiPAnalyzeR default
Res_LiPAnalyzeR <- sumModelsAllStrains(RUV_TrPPepTrPProtB_CON_Strain, 
                                       annotPP = annotPepProt)
### LiPAnalyzeR default, using different strains as reference strains
Res_LiPAnalyzeR_refJB50 <- sumModelsAllStrains(RUV_TrPPepTrPProtB_CON_Strain_refJB50, 
                                               annotPP = annotPepProt)

Res_LiPAnalyzeR_refJB759 <- sumModelsAllStrains(RUV_TrPPepTrPProtB_CON_Strain_refJB759, 
                                                S1 = "StrainJB50_Contrast", 
                                                listNames = c("JB50", "JB760","PYK1A"),
                                                annotPP = annotPepProt)

Res_LiPAnalyzeR_refJB760 <- sumModelsAllStrains(RUV_TrPPepTrPProtB_CON_Strain_refJB760, 
                                                S1 = "StrainJB50_Contrast", S2 = "StrainJB759_Contrast",
                                                listNames = c("JB50", "JB759","PYK1A"),
                                                annotPP = annotPepProt)

Res_LiPAnalyzeR_refPYK1A <- sumModelsAllStrains(RUV_TrPPepTrPProtB_CON_Strain_refPYK1A, 
                                                S1 = "StrainJB50_Contrast", S2 = "StrainJB759_Contrast",
                                                S3 = "StrainJB760_Contrast",
                                                listNames = c("JB50", "JB759","JB760"),
                                                annotPP = annotPepProt)

### LiPAnalyzeR LiPonly 
Res_LiPAnalyzeR_LiPonly <- sumModelsAllStrains(RUV_LiPProtB_CON_Strain_LiPonly, 
                                               annotPP = annotPepProt)

### LiPAnalyzeR HT default: BVLS_TrPProtB_OLSStrain
Res_LiPAnalyzeR_FTHT <- sumModelsAllStrains(RUV_TrPProtB_CON_Strain_FTHTjoin,
                                            annotPP = annotPepProt)

### LiPAnalyzeR HTonly default: BVLS_TrPPepTrPProtB_OLSStrain
Res_LiPAnalyzeR_HTonly <- sumModelsAllStrains(RUV_TrPPepTrPProtB_CON_Strain_HTonly,
                                              annotPP = annotPepProt)

### LiPAnalyzeR, no constraints in RUV
Res_LiPAnalyzeR_RUVnC_CON <- sumModelsAllStrains(RUVnC_TrPPepTrPProtB_CON_Strain,
                                                 annotPP = annotPepProt)

### LiPAnalyzeR, no constraints in RUV, RUV and contrast in one model
Res_LiPAnalyzeR_RUVnC_oneModel <- sumModelsAllStrains(RUVnC_TrPPepTrPProtBStrain,
                                                      annotPP = annotPepProt)




## Plotting correlations of residuals with TrPPep/TrPProt
### Correlations to residuals of RUV corrected for only TrPPep or TrPProt
peps <- row.names(RUV_TrPPepTrPProtB_CON_Strain$modelCoeff)[RUVnC_TrPPepTrPProtB_CON_Strain$modelCoeff$XProt_RUV > (-1e+06) &
                                                                RUVnC_TrPPepTrPProtB_CON_Strain$modelCoeff$XPep_RUV > (-1e+06)]
plotData <- melt(data.frame(RUV_TrPProt = CorCols(peps,
			    			  RUV_TrPPepB$modelResid, 
                               			  QuantityListBC$TrPProt) ,
                            RUV_TrPPep = CorCols(peps, 
                                                 RUV_TrPProtB$modelResid,
                                                 QuantityListBC$TrPPep)))
#### Figure 5c
ggplot(plotData, aes(x = value, color = variable, fill = variable, alpha = variable)) +
    geom_density(linewidth = 0.7) +
    ylab("Density") +
    xlab("Correlation coefficient") +
    xlim(-1, 1) +
    geom_vline(xintercept = median(plotData[plotData$variable == "RUV_TrPProt", "value"]), 
               color = "#002144", linewidth = 1) +
    geom_vline(xintercept = median(plotData[plotData$variable == "RUV_TrPPep", "value"]), 
               color = "#a07d00", linewidth = 1) +
    scale_color_manual(values = c("#002144", "#a07d00")) +
    scale_fill_manual(values = c("#002144", "#a07d00")) +
    scale_alpha_manual(values = c(0.8, 0.4)) +
    theme_bw(base_size = 11)

### Correlations to residuals of RUV corrected for both TrPPep or TrPProt
plotData <- melt(data.frame(RUV_TrPProt = CorCols(row.names(RUV_TrPPepTrPProtB$modelResid),
                                                  RUV_TrPPepTrPProtB$modelResid, 
                                                  QuantityList$TrPProt),
                            RUV_TrPPep = CorCols(row.names(RUV_TrPPepTrPProtB$modelResid), 
                                                 RUV_TrPPepTrPProtB$modelResid,
                                                 QuantityList$TrPPep)))

#### Figure 5e
ggplot(plotData, aes(x = value, color = variable, fill = variable, alpha = variable)) +
    geom_histogram(linewidth = 0.7, position = "identity", binwidth = 0.05) +
    ylab("Density") +
    xlab("Correlation coefficient") +
    xlim(-1, 1) +
    geom_vline(xintercept = median(plotData[plotData$variable == "RUV_TrPProt", "value"]), 
               color = "#002144", linewidth = 1) +
    geom_vline(xintercept = median(plotData[plotData$variable == "RUV_TrPPep", "value"]), 
               color = "#a07d00", linewidth = 1, linetype = "dashed") +
    scale_color_manual(values = c("#002144", "#a07d00")) +
    scale_fill_manual(values = c("#002144", "#a07d00")) +
    scale_alpha_manual(values = c(0.8, 0.4)) +
    theme_bw(base_size = 11) 




## Plotting TrPPep and TrPProt coefficients from different models
### Coefficients of TrPPep & TrPProt from RUV without constraints
plotData <- RUVnC_TrPPepTrPProtB$modelCoeff
P1 <- ggplot(mapping = aes(x = plotData$XProt_Contrast,
                           y = plotData$XPep_Contrast, 
                           color = plotData$XProt_Contrast<0|
                               plotData$XPep_Contrast<0)) +
    geom_point(size = 1, alpha = 0.1) +
    scale_color_manual(values = c("#000000", "#c94202")) +
    xlim(-18, 22) +
    ylim(-18, 22) +
    ylab("Coef RUV - TrPPep") +
    xlab("Coef RUV - TrPProt") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
P2 <- ggplot(mapping = aes(x = plotData$XPep_Contrast)) +
    geom_histogram(breaks = seq(-18, 22, 0.5)) +
    xlim(-18, 22) +
    coord_flip() +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = -90))
P3 <- ggplot(mapping = aes(x = plotData$XProt_Contrast)) +
    geom_histogram(breaks = seq(-18, 22, 0.5)) +
    xlim(-18, 22) +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
#### Figure 6a
P3 + plot_spacer() + P1 + P2 +
    plot_layout(nrow = 2, ncol = 2, widths = c(2, 1), heights = c(1, 2))

### Coefficients TrPPep & TrPProt from RUV with constraints
plotData <- RUV_TrPPepTrPProtB$modelCoeff

P1 <- ggplot(plotData, aes(x = XProt_RUV,
                           y = XPep_RUV)) +
    geom_point(size = 1, alpha = 0.1) +
    xlim(0,5) +
    ylim(0,5) +
    ylab("Coef RUV - TrPPep") +
    xlab("Coef RUV - TrPProt") +
    theme_bw(base_size = 11)
P2 <- ggplot(plotData, aes(x = XPep_RUV)) +
    geom_histogram(breaks = seq(0, 5, 0.1)) +
    xlim(0, 5) +
    coord_flip() +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = -90))
P3 <- ggplot(plotData, aes(x = XProt_RUV)) +
    geom_histogram(breaks = seq(0, 5, 0.1)) +
    xlim(0, 5) +
    xlab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
#### Figure 6b
P3 + plot_spacer() + P1 + P2 +
    plot_layout(nrow = 2, ncol = 2, widths = c(2, 1), heights = c(1, 2))


### Coefficients TrPPep & TrPProt from RUV with constraints run in HTonly mode
plotData <- RUV_TrPPepTrPProtB_HTonly$modelCoeff

P1 <- ggplot(plotData, aes(x = XProt_RUV,
                           y = XPep_RUV)) +
    geom_point(size = 1, alpha = 0.1) +
    xlim(0, 5) +
    ylim(0, 5) +
    ylab("Coef RUV - TrPPep") +
    xlab("Coef RUV - TrPProt") +
    theme_bw(base_size = 11)
P2 <- ggplot(plotData, aes(x = XPep_RUV)) +
    geom_histogram(breaks = seq(0, 5, 0.1)) +
    xlim(0, 5) +
    coord_flip() +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = -90))
P3 <- ggplot(plotData, aes(x = XProt_RUV)) +
    geom_histogram(breaks = seq(0, 5, 0.1)) +
    xlim(0, 5) +
    xlab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
#### Figure 7b
P3 + plot_spacer() + P1 + P2 +
    plot_layout(nrow = 2, ncol = 2, widths = c(2, 1), heights = c(1, 2))


### TrPPep coefficient RUV with vs RUV without constraints
#### Extended Data Figure 2c
ggplot(mapping = aes(y = RUV_TrPPepTrPProtB$modelCoeff[peps, "XPep_RUV"],
                     x = RUVnC_TrPPepTrPProtB$modelCoeff[peps, "XPep_Contrast"])) +
    geom_point(size = 3, alpha = 0.1) +
    xlab("Coef RUV (no contraints) - TrPPep") +
    ylab("Coef RUV - TrPPep") +
    theme_bw(base_size = 15) +
    theme(legend.position = "none")

### TrPProt coefficient RUV with vs RUV without constraints
#### Extended Data Figure 2d
ggplot(mapping = aes(y = RUV_TrPPepTrPProtB$modelCoeff[peps, "XProt_RUV"],
                     x = RUVnC_TrPPepTrPProtB$modelCoeff[peps, "XProt_Contrast"])) +
    geom_point(size = 3, alpha = 0.1) +
    xlab("Coef RUV (no contraints) - TrPProt") +
    ylab("Coef RUV - TrPProt") +
    theme_bw(base_size = 15) +
    theme(legend.position = "none")


### LiPAnalyzeR run in (RUV 1): combined model without constraints; (RUV2): two models without constraints; (RUV3): two models with constraints
#### Extended Data Figure 4a, 4c, 4e
nameS <- "StrainJB759_Contrast" # Extended Data Figure 4a
#nameS <- "StrainJB760_Contrast"  # Extended Data Figure 4c
#nameS <- "StrainPYK1_A_Contrast"  # Extended Data Figure 4e
P1 <- ggplot(mapping = aes(y = RUV_TrPPepTrPProtB_CON_Strain$modelCoeff[peps, nameS],
                           x = RUVnC_TrPPepTrPProtBStrain$modelCoeff[peps, nameS])) +
    geom_point(size = 3, alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, color = "lightblue", linewidth = 1) +
    #geom_smooth(formula = "y~x", method = "lm", color = "#c94202", linewidth = 1) +
    xlab("RUV (1)") +
    ylab("RUV (3)") +
    theme_bw(base_size = 15) +
    theme(aspect.ratio = 1)

P2 <- ggplot(mapping = aes(y = RUV_TrPPepTrPProtB_CON_Strain$modelCoeff[peps, nameS],
                           x = RUVnC_TrPPepTrPProtB_CON_Strain$modelCoeff[peps, nameS])) +
    geom_point(size = 3, alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, color = "lightblue", linewidth = 1) +
    #geom_smooth(formula = "y~x", method = "lm", color = "#c94202", linewidth = 1) +
    xlab("RUV (2)") +
    ylab("RUV (3)") +
    theme_bw(base_size = 15) +
    theme(aspect.ratio = 1)

P1 + P2 + plot_layout(ncol = 2)





## Plotting strain coefficients
### LiPAnalyzeR run on half-tryptic peptides in FTHTjoin and in HTonly mode
HTpeps <- intersect(row.names(Res_LiPAnalyzeR_FTHT$JB759),
                    row.names(Res_LiPAnalyzeR_HTonly$JB759))
doScatter <- function(xAxis, yAxis, xName, yName, xLim, yLim){
    ggplot(mapping = aes(x = xAxis, y = yAxis)) +
        geom_point(alpha = 0.3, size = 1.5) +
        geom_abline(intercept = 0, slope = 1, color = "lightblue", linewidth = 0.7) +
        geom_smooth(formula = "y~x", method = "lm", color = "#c94202", linewidth = 0.7) +
        xlab(xName) +
        ylab(yName) +
        ylim(yLim) +
        xlim(xLim) +
        theme_bw(base_size = 11) +
        theme(aspect.ratio = 1)
}

#### Extended Data Figure 5d, 5e, 5f
Strain <- "JB759" # Extended Data Figure 5d
Strain <- "JB760" # Extended Data Figure 5e
Strain <- "PYK1A" # Extended Data Figure 5f

P1 <- doScatter(Res_LiPAnalyzeR_HTonly[[Strain]][HTpeps, "Coefficient"],
                Res_LiPAnalyzeR_FTHT[[Strain]][HTpeps, "Coefficient"],
                "HT only - Coefficient", "HT - Coeffficient",
                c(-5, 5), c(-5,5))
P2 <- doScatter(-log10(Res_LiPAnalyzeR_HTonly[[Strain]][HTpeps, "Pvalue"]),
                -log10(Res_LiPAnalyzeR_FTHT[[Strain]][HTpeps, "Pvalue"]),
                "HT only - -log10(p-value)","HT - -log10(p-value)",
                c(0, 10), c(0,10))

P1 + P2 + plot_layout(ncol = 2)

### Comparison strain coefficients estimated in contrast model with (1) RUV LiPAnalyzeR or (2) LiP/TrP ratio for correcting for unwanted variation
Res_CON_Strain_RatioTrPPep <- sumModelsAllStrains(CON_LiPPepTrPPepRatio,
                                                  annotPP = annotPepProt)
Res_CON_Strain_RatioTrPProt <- sumModelsAllStrains(CON_LiPPepTrPProtRatio,
                                                   annotPP = annotPepProt)

compareHits <- function(Res1, Res2, pCol = "Padj", tPval = 0.05){
    isect <- intersect(row.names(Res1), row.names(Res2))
    Res1 <- Res1[isect, pCol] < tPval
    Res2 <- Res2[isect, pCol] < tPval
    comVec <- sapply(seq(1, length(isect)), \(i){
        if(!Res1[i] & !Res2[i]){"NN"}
        else if(Res1[i] & !Res2[i]){"PN"}
        else if(!Res1[i] & Res2[i]){"NP"}
        else{"PP"}
    })
    names(comVec) <- isect
    return(comVec)
}

compare_LARdefault_TPepR <- list(JB759 = compareHits(Res_LiPAnalyzeR$JB759, 
                                                     Res_CON_Strain_RatioTrPPep$JB759),
                                 JB760 = compareHits(Res_LiPAnalyzeR$JB760, 
                                                     Res_CON_Strain_RatioTrPPep$JB760),
                                 PYK1A = compareHits(Res_LiPAnalyzeR$PYK1A, 
                                                     Res_CON_Strain_RatioTrPPep$PYK1A))

compare_LARdefault_TProtR <- list(JB759 = compareHits(Res_LiPAnalyzeR$JB759, 
                                                      Res_CON_Strain_RatioTrPProt$JB759),
                                  JB760 = compareHits(Res_LiPAnalyzeR$JB760, 
                                                      Res_CON_Strain_RatioTrPProt$JB760),
                                  PYK1A = compareHits(Res_LiPAnalyzeR$PYK1A, 
                                                      Res_CON_Strain_RatioTrPProt$PYK1A))

createPlot_comparingLARRatio <- function(x, y, plotTitle, xlabels){
    ggplot(mapping = aes(x = unname(x), y = y[names(x)])) +
        geom_violin() +
        geom_boxplot(width = 0.5, alpha = 0.5, fill = "lightgrey") +
        scale_x_discrete(name = "", labels = xlabels) +
        ylab("Correlation coefficient") +
        ggtitle(plotTitle) + 
        theme_bw(base_size = 15)
}

labels_TPepR <- c("LAR: NH\nTPepR: NH", "LAR: NH\nTPepR: H",
                  "LAR: H\nTPepR: NH", "LAR: H\nTPepR: H")
labels_TProtR <- c("LAR: NH\nTProtR: NH", "LAR: NH\nTProtR: H",
                   "LAR: H\nTProtR: NH", "LAR: H\nTProtR: H")

#### Extended Data Figure 1c
createPlot_comparingLARRatio(compare_LARdefault_TPepR$JB759, 
                             cRLiPPepTrPPep_TrPPep,
                             "JB759", labels_TPepR) +
    createPlot_comparingLARRatio(compare_LARdefault_TPepR$JB760, 
                                 cRLiPPepTrPPep_TrPPep,
                                 "JB760", labels_TPepR) +
    createPlot_comparingLARRatio(compare_LARdefault_TPepR$PYK1A, 
                                 cRLiPPepTrPPep_TrPPep,
                                 "PYK1A", labels_TPepR) + plot_layout(nrow = 3)

#### Extended Data Figure 1d
createPlot_comparingLARRatio(compare_LARdefault_TProtR$JB759, 
                             cRLiPPepTrPProt_TrPProt,
                             "JB759", labels_TProtR) +
    createPlot_comparingLARRatio(compare_LARdefault_TProtR$JB760, 
                                 cRLiPPepTrPProt_TrPProt,
                                 "JB760", labels_TProtR) +
    createPlot_comparingLARRatio(compare_LARdefault_TProtR$PYK1A, 
                                 cRLiPPepTrPProt_TrPProt,
                                 "PYK1A", labels_TProtR) + plot_layout(nrow = 3)


# Woods plots
## Creating annotPepProt for Plotting
annotPepProtPlot <- annotPepProt
annotPepProtPlot$Protein <- sapply(annotPepProtPlot$Protein, \(x){
    x <- unlist(strsplit(x, ";"))
    x <- gsub(">>", "", gsub("_.*", "", x))
    x <- unique(x)
    x <- paste0(x, collapse = ";")
    return(x)
})

annotPepProtPlot$startPosition <- sapply(annotPepProtPlot$startPosition, \(x){
    x <- unlist(strsplit(x, ";"))
    x <- unique(x)
    x <- paste0(x, collapse = ";")
    return(x)
})

annotPepProtPlot$endPosition <- sapply(annotPepProtPlot$endPosition, \(x){
    x <- unlist(strsplit(x, ";"))
    x <- unique(x)
    x <- paste0(x, collapse = ";")
    return(x)
})

## Overview FT and HT over protein residue
### Example protein
i <- "SPAC26A3.07c"
exampleProt <- annotPepProtPlot[sapply(annotPepProtPlot$Protein, \(x) grepl(i, x)),]
exampleProt$startPosition <- as.numeric(exampleProt$startPosition)
exampleProt$endPosition <- as.numeric(exampleProt$endPosition)
exampleProt$isTryptic <- ifelse(exampleProt$isTryptic == "Specific", "FT", "HT")
exampleProt <- exampleProt[order(exampleProt$endPosition, decreasing = F),]
exampleProt <- exampleProt[order(exampleProt$startPosition, decreasing = F),]
exampleProt$yPosition <- c(1:nrow(exampleProt))

#### Figure 7a
ggplot(exampleProt, aes(x = startPosition, xend = endPosition, y = yPosition, 
                        yend = yPosition, fill = isTryptic, color = isTryptic)) +
    geom_segment(linewidth = 3) +
    scale_color_manual(values = c("#000000", "#a6a6a6")) +
    scale_fill_manual(values = c("#000000", "#a6a6a6")) +
    xlab("Protein residue") +
    xlim(0, 174) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())


## Woods plot with strain coefficients estimated by LiPAnalyzeR on Y-axis
### Example for half-tryptic peptides analyzed with HTonly  
comRes_FTdefault_HTonly <- list(JB759 = rbind(Res_LiPAnalyzeR$JB759, Res_LiPAnalyzeR_HTonly$JB759),
                                JB760 = rbind(Res_LiPAnalyzeR$JB760, Res_LiPAnalyzeR_HTonly$JB760),
                                PYK1A = rbind(Res_LiPAnalyzeR$PYK1A, Res_LiPAnalyzeR_HTonly$PYK1A))


#### Figure 7c
i <- "SPAC14C4.14"
makeWoodsPlotSingleProtein(sumDf = comRes_FTdefault_HTonly$JB759,
                           annotPP = annotPepProtPlot, 
                           deltaColorIsTryptic = T,
                           protName = i)

#### Figure 7d
i <- "SPAP8A3.07c"
makeWoodsPlotSingleProtein(sumDf = comRes_FTdefault_HTonly$JB759,
                           annotPP = annotPepProtPlot, 
                           deltaColorIsTryptic = T,
                           protName = i, ylim = c(-2.1, 2.1))

makeWoodsPlotSingleProtein(sumDf = comRes_FTdefault_HTonly$JB760,
                           annotPP = annotPepProtPlot, 
                           deltaColorIsTryptic = T,
                           protName = i, ylim = c(-2.1, 2.1))


### Running LiPAnalyzeR in default and in LiPonly mode
i <- "SPAC23G3.01"
#### Figure 8b
makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR$JB759,
                           annotPP = annotPepProtPlot, 
                           deltaColorIsTryptic = T,
                           protName = i, ylim = c(-2.1, 2.1))

#### Figure 8c
makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR_LiPonly$JB759,
                           annotPP = annotPepProtPlot, 
                           deltaColorIsTryptic = T,
                           protName = i, ylim = c(-2.1, 2.1))


### PYK1 with different reference strains
i <- "SPAC4H3.10c"
#### Extended Data Figure 6b (reference strain JB50)
plotList_refJB50 <- list(makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR$JB759,
                                                    annotPP = annotPepProtPlot, protName = i,
                                                    ylim = c(-2.25, 2.25)),
                         makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR$JB760,
                                                    annotPP = annotPepProtPlot, protName = i,
                                                    ylim = c(-2.25, 2.25)),
                         makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR$PYK1A,
                                                    annotPP = annotPepProtPlot, protName = i,
                                                    ylim = c(-2.25, 2.25)))
plotList_refJB50[[1]]$labels$title <- "Reference strain JB50, analyzed strain JB759"
plotList_refJB50[[2]]$labels$title <- "Reference strain JB50, analyzed strain JB760"
plotList_refJB50[[3]]$labels$title <- "Reference strain JB50, analyzed strain PYK1A"

plotList_refJB50[[1]] + theme(legend.position = "none") + 
    plotList_refJB50[[2]] + theme(legend.position = "none") +
    plotList_refJB50[[3]] + theme(legend.position = "none") +
    plot_layout(nrow = 3) 

#### Extended Data Figure 6c (reference strain PYK1A)
plotList_refPYK1A <- list(makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR_refPYK1A$JB50,
                                                     annotPP = annotPepProtPlot, protName = i,
                                                     ylim = c(-1.25, 3.5)),
                          makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR_refPYK1A$JB759,
                                                     annotPP = annotPepProtPlot, protName = i,
                                                     ylim = c(-1.25, 3.5)),
                          makeWoodsPlotSingleProtein(sumDf = Res_LiPAnalyzeR_refPYK1A$JB760,
                                                     annotPP = annotPepProtPlot, protName = i,
                                                     ylim = c(-1.25, 3.5)))
plotList_refPYK1A[[1]]$labels$title <- "Reference strain JB50 (pyk1a), analyzed strain JB50"
plotList_refPYK1A[[2]]$labels$title <- "Reference strain JB50 (pyk1a), analyzed strain JB759"
plotList_refPYK1A[[3]]$labels$title <- "Reference strain JB50 (pyk1a), analyzed strain JB760"

plotList_refPYK1A[[1]] + theme(legend.position = "none") + 
    plotList_refPYK1A[[2]] + theme(legend.position = "none") +
    plotList_refPYK1A[[3]] + theme(legend.position = "none") +
    plot_layout(nrow = 3) 
