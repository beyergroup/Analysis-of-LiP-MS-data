library(corrplot)
library(LiPAnalyzeR)
library(ggplot2)
library(ggpubr)
library(limma)
library(patchwork)
library(reshape2)

# Loading data & some general sorting
LiPSpectro <- read.csv("../Documents/LiPAnalyzer/Data/RMBY/ParentalStrains/SpectroOut/270722_BYRM_fastaCombined_LiP_Spectro_LiPAnalyzer_LN.csv")
LiPSpectro <- LiPSpectro[LiPSpectro$PG.ProteinGroups != "iRT_Peptides_Fusion",]
TrPSpectro <- read.csv("../Documents/LiPAnalyzer/Data/RMBY/ParentalStrains/SpectroOut/270722_BYRM_fastaCombined_Trp_Spectro_LiPAnalyzer_LN.csv")
TrPSpectro <- TrPSpectro[TrPSpectro$PG.ProteinGroups != "iRT_Peptides_Fusion",]

annotSample <- readRDS("../Documents/LiPAnalyzer/Data/RMBY/ParentalStrains/SampleInfo/info_SumSamples.rds")
row.names(annotSample) <- gsub("_", "_B", row.names(annotSample))

annotLiP <- readRDS("../Documents/LiPAnalyzer/Data/RMBY/ParentalStrains/SampleInfo/info_LiP_wOrgName.rds")
annotTrP <- readRDS("../Documents/LiPAnalyzer/Data/RMBY/ParentalStrains/SampleInfo/info_Trp_wOrgName.rds")
annotSampleReps <- annotLiP[,c(2:5)]
annotSampleReps$RunOrderLiP <- gsub("Q_", "", annotLiP$ID)
annotSampleReps$RunOrderTrP <- gsub("Q_", "", annotTrp$ID)


# Extracting peptide and protein quantities 
## Extracting LiPPep, TrPPep and TrPProt into a list
QuantityList_sReps <- extractSpectroData(spectroLiP = LiPSpectro,
                                         spectroTrP = TrPSpectro)
QuantityListLiP_sReps <- extractSpectroData(spectroLiP = LiPSpectro,
                                            LiPonly = T)

### Adjusting colnames of data
colnames(QuantityList_sReps$LiPPep) <- row.names(annotLiP)[match(colnames(QuantityList_sReps$LiPPep),
                                                                 annotLiP$OrgSampleNames)]
colnames(QuantityList_sReps$TrPPep) <- row.names(annotTrP)[match(colnames(QuantityList_sReps$TrPPep),
                                                                 annotTrP$OrgSampleNames)]
colnames(QuantityList_sReps$TrPProt) <- row.names(annotTrP)[match(colnames(QuantityList_sReps$TrPProt),
                                                                  annotTrP$OrgSampleNames)]

colnames(QuantityListLiP_sReps$LiPPep) <- row.names(annotLiP)[match(colnames(QuantityListLiP_sReps$LiPPep),
                                                                    annotLiP$OrgSampleNames)]
colnames(QuantityListLiP_sReps$LiPProt) <- row.names(annotLiP)[match(colnames(QuantityListLiP_sReps$LiPProt),
                                                                     annotLiP$OrgSampleNames)]

# Further sorting of the data
## Creating data.frame with information to all peptides in data
annotPepProt <- getPepProtAnnotSpectro(spectroOut = TrPSpectro, 
                                       spectroOut2 = LiPSpectro)

## Adding together technical replicates, mean if two measurements, if one is NA other LiPProt
SumTechReps <- function(SpectroList){
    lapply(SpectroList, \(x){
        colnames(x) <- sapply(colnames(x), \(y) {
            paste(strsplit(y, "_")[[1]][c(1,3)], collapse = "_")
        })
        x <- melt(as.matrix(x))
        x <- dcast(x, Var1~Var2, fun.aggregate = mean, na.rm = T)
        row.names(x) <- x$Var1
        x <- x[,-1]
        x[is.nan(as.matrix(x))] <- NA
        return(x)
    })
}
QuantityList <- SumTechReps(QuantityList_sReps)
QuantityListLiP <- SumTechReps(QuantityListLiP_sReps)

## Filtering quantities
QuantityList_sReps <- preprocessQuantityMatrix(quantityList = QuantityList_sReps, 
                                               annotS = annotSampleReps,
                                               annotPP = annotPepProt)
QuantityList_sReps <- lapply(QuantityList_sReps, \(x){
    x[, row.names(annotSampleReps)]
})

QuantityListLiP_sReps <- preprocessQuantityMatrix(quantityList = QuantityListLiP_sReps,
                                                  annotS = annotSampleReps,
                                                  annotPP = annotPepProt,
                                                  filterTryptic = F)
QuantityListLiP_sReps <- lapply(QuantityListLiP_sReps, \(x){
    x[, row.names(annotSampleReps)]
})

QuantityList <- preprocessQuantityMatrix(quantityList = QuantityList, 
                                         annotS = annotSample,
                                         annotPP = annotPepProt)
QuantityList <- lapply(QuantityList, \(x){
    x[, row.names(annotSample)]
})

QuantityListLiP <- preprocessQuantityMatrix(quantityList = QuantityListLiP,
                                            annotS = annotSample,
                                            annotPP = annotPepProt,
                                            filterTryptic = F)
QuantityListLiP <- lapply(QuantityListLiP, \(x){
    x[, row.names(annotSample)]
})

## Batch correct data
QuantityList_sRepsBC <- lapply(QuantityList_sReps, \(x){
    removeBatchEffect(x, annotSampleReps$Batch)
})

QuantityListLiP_sRepsBC <- lapply(QuantityListLiP_sReps, \(x){
    removeBatchEffect(x, annotSampleReps$Batch)
})

QuantityListBC <- lapply(QuantityList, \(x){
    removeBatchEffect(x, annotSample$Batch)
})

QuantityListLiPBC <- lapply(QuantityListLiP, \(x){
    removeBatchEffect(x, annotSample$Batch)
})





# LiPAnalyzeR
## Running models on technical replicates or 'biological' groups
### RUV without constraints, modeling TrPPep, TrPProt and strain in one model
RUVnC_TrPPepTrPProtStrain_R1 <- runModel(quantityList = QuantityList_sReps,
                                         annotS = annotSampleReps[annotSampleReps$Replicate == "R1", ],
                                         formulaRUV = NULL,
                                         formulaContrast = "Y~XPep+XProt+Strain")
RUVnC_TrPPepTrPProtStrain_R2 <- runModel(quantityList = QuantityList_sReps,
                                         annotS = annotSampleReps[annotSampleReps$Replicate == "R2", ],
                                         formulaRUV = NULL,
                                         formulaContrast = "Y~XPep+XProt+Strain")

RUVnC_TrPPepTrPProtStrain_G1 <- runModel(quantityList = QuantityList,
                                         annotS = annotSample[annotSample$Batch %in% c("01","02", "03", "04", "05"), ],
                                         formulaRUV = NULL,
                                         formulaContrast = "Y~XPep+XProt+Strain")
RUVnC_TrPPepTrPProtStrain_G2 <- runModel(quantityList = QuantityList,
                                         annotS = annotSample[annotSample$Batch %in% c("06","07", "08", "09", "10"), ],
                                         formulaRUV = NULL,
                                         formulaContrast = "Y~XPep+XProt+Strain")

## OLS TrpPep TrpProt, OLS Strain
### RUV without constraints, modeling TrPPep, TrPProt and contrast model with strain, two models
RUVnC_TrPPepTrPProt_CON_Strain_R1 <- runModel(quantityList = QuantityList_sReps,
                                              annotS = annotSampleReps[annotSampleReps$Replicate == "R1", ],
                                              formulaRUV = "Y~XPep+XProt",
                                              formulaContrast = "Y~Strain",
                                              lowRUV = c(-Inf, -Inf, -Inf),
                                              upRUV = c(Inf, Inf, Inf))
RUVnC_TrPPepTrPProt_CON_Strain_R2 <- runModel(quantityList = QuantityList_sReps,
                                              annotS = annotSampleReps[annotSampleReps$Replicate == "R2", ],
                                              formulaRUV = "Y~XPep+XProt",
                                              formulaContrast = "Y~Strain",
                                              lowRUV = c(-Inf, -Inf, -Inf),
                                              upRUV = c(Inf, Inf, Inf))

RUVnC_TrPPepTrPProt_CON_Strain_G1 <- runModel(quantityList = QuantityList,
                                              annotS = annotSample[annotSample$Batch %in% c("01","02", "03", "04", "05"), ],
                                              formulaRUV = "Y~XPep+XProt",
                                              formulaContrast = "Y~Strain",
                                              lowRUV = c(-Inf, -Inf, -Inf),
                                              upRUV = c(Inf, Inf, Inf))
RUVnC_TrPPepTrPProt_CON_Strain_G2 <- runModel(quantityList = QuantityList,
                                              annotS = annotSample[annotSample$Batch %in% c("06","07", "08", "09", "10"), ],
                                              formulaRUV = "Y~XPep+XProt",
                                              formulaContrast = "Y~Strain",
                                              lowRUV = c(-Inf, -Inf, -Inf),
                                              upRUV = c(Inf, Inf, Inf))

### RUV with constraints, modeling TrPPep, TrPProt and contrast model with strain, two models
RUV_TrPPepTrPProt_CON_Strain_R1 <- runModel(quantityList = QuantityList_sReps,
                                            annotS = annotSampleReps[annotSampleReps$Replicate == "R1", ],
                                            formulaRUV = "Y~XPep+XProt",
                                            formulaContrast = "Y~Strain")
RUV_TrPPepTrPProt_CON_Strain_R2 <- runModel(quantityList = QuantityList_sReps,
                                            annotS = annotSampleReps[annotSampleReps$Replicate == "R2", ],
                                            formulaRUV = "Y~XPep+XProt",
                                            formulaContrast = "Y~Strain")

RUV_TrPPepTrPProt_CON_Strain_G1 <- runModel(quantityList = QuantityList,
                                            annotS = annotSample[annotSample$Batch %in% c("01","02", "03", "04", "05"), ],
                                            formulaRUV = "Y~XPep+XProt",
                                            formulaContrast = "Y~Strain")
RUV_TrPPepTrPProt_CON_Strain_G2 <- runModel(quantityList = QuantityList,
                                            annotS = annotSample[annotSample$Batch %in% c("06","07", "08", "09", "10"), ],
                                            formulaRUV = "Y~XPep+XProt",
                                            formulaContrast = "Y~Strain")


## Comparing different models
### Scatterplots
DoScatter <- function(x, y, nameX, nameY, xlim = NULL){
    if(is.null(xlim)){
        min <- min(min(na.omit(x)), min(na.omit(y)))
        max <- max(max(na.omit(x)), max(na.omit(y)))
    }
    ggplot(mapping =  aes(x = x, y = y)) +
        geom_point(size = 2, alpha = 0.3) +
        geom_abline(slope = 1, color = "lightblue", linewidth = 1) +
        geom_smooth(formula = "y~x", method = "lm", color = "#c94202", linewidth = 1) +
        ylab(nameY) +
        xlab(nameX) +
        xlim(min, max) +
        theme_bw(base_size = 11) +
        theme(aspect.ratio = 1)
}

P1 <- DoScatter(RUVnC_TrPPepTrPProtStrain_R1$modelCoeff$`StrainRM11-1a_Contrast`,
                RUVnC_TrPPepTrPProtStrain_R2$modelCoeff$`StrainRM11-1a_Contrast`,
                "RUV (1) R1", "RUV (1) R2")
P2 <- DoScatter(RUVnC_TrPPepTrPProt_CON_Strain_R1$modelCoeff$Strain_Contrast,
                RUVnC_TrPPepTrPProt_CON_Strain_R2$modelCoeff$Strain_Contrast,
                "RUV (2) R1", "RUV (2) R2")
P3 <- DoScatter(RUV_TrPPepTrPProt_CON_Strain_R1$modelCoeff$Strain_Contrast,
                RUV_TrPPepTrPProt_CON_Strain_R2$modelCoeff$Strain_Contrast,
                "RUV (3) R1", "RUV (3) R2")
#### Figure 6e
P1 + P2 + P3

P1 <- DoScatter(RUVnC_TrPPepTrPProtStrain_G1$modelCoeff$`StrainRM11-1a_Contrast`,
                RUVnC_TrPPepTrPProtStrain_G2$modelCoeff$`StrainRM11-1a_Contrast`,
                "RUV (1) G1", "RUV (1) G2")
P2 <- DoScatter(RUVnC_TrPPepTrPProt_CON_Strain_G1$modelCoeff$Strain_Contrast,
                RUVnC_TrPPepTrPProt_CON_Strain_G2$modelCoeff$Strain_Contrast,
                "RUV (2) G1", "RUV (2) G2")
P3 <- DoScatter(RUV_TrPPepTrPProt_CON_Strain_G1$modelCoeff$Strain_Contrast,
                RUV_TrPPepTrPProt_CON_Strain_G2$modelCoeff$Strain_Contrast,
                "RUV (3) G1", "RUV (3) G2")
#### Extended Data Figure 3b
P1 + P2 + P3

### Creating correlation plots for technical replicates
replicatesCorPlot <- function(cName1, cName2, col = colorRampPalette(c("#ffffff", "#c94202"))){
    peps <- union(row.names(RUVnC_TrPPepTrPProtStrain_R1$modelCoeff)[RUVnC_TrPPepTrPProtStrain_R1$modelCoeff$XProt_Contrast > (-1e+06) &
                                                                         RUVnC_TrPPepTrPProtStrain_R1$modelCoeff$XPep_Contrast> (-1e+06)],
                  row.names(RUVnC_TrPPepTrPProtStrain_R2$modelCoeff)[RUVnC_TrPPepTrPProtStrain_R2$modelCoeff$XProt_Contrast > (-1e+06) &
                                                                         RUVnC_TrPPepTrPProtStrain_R2$modelCoeff$XPep_Contrast > (-1e+06)])
    mat <- data.frame("RUV(1)R1" = RUVnC_TrPPepTrPProtStrain_R1$modelCoeff[peps, cName1],
                      "RUV(2)R1" = RUVnC_TrPPepTrPProt_CON_Strain_R1$modelCoeff[peps, cName2],
                      "RUV(3)R1" = RUV_TrPPepTrPProt_CON_Strain_R1$modelCoeff[peps, cName2],
                      "RUV(1)R2" = RUVnC_TrPPepTrPProtStrain_R2$modelCoeff[peps, cName1],
                      "RUV(2)R2" = RUVnC_TrPPepTrPProt_CON_Strain_R2$modelCoeff[peps, cName2],
                      "RUV(3)R2" = RUV_TrPPepTrPProt_CON_Strain_R2$modelCoeff[peps, cName2])
    cor <- cor(mat, use = "complete.obs")
    corrplot(cor, type = "upper", method = "color", addCoef.col = 'black', cl.pos = "n",
             tl.col = "black", tl.srt = 45, is.corr = F, col = col(200))
}

#### Extended Data Figure 3a
replicatesCorPlot("StrainRM11-1a_Contrast", "Strain_Contrast")
replicatesCorPlot("XPep_Contrast", "XPep_RUV")
replicatesCorPlot("XProt_Contrast", "XProt_RUV")


### Creating correlation plots for 'biological' groups
biologicalCorPlot <- function(cName1, cName2, col = colorRampPalette(c("#ffffff", "#c94202"))){
    peps <- union(row.names(RUVnC_TrPPepTrPProtStrain_G1$modelCoeff)[RUVnC_TrPPepTrPProtStrain_G1$modelCoeff$XProt_Contrast > (-1e+06) &
                                                                         RUVnC_TrPPepTrPProtStrain_G1$modelCoeff$XPep_Contrast> (-1e+06)],
                  row.names(RUVnC_TrPPepTrPProtStrain_G2$modelCoeff)[RUVnC_TrPPepTrPProtStrain_G2$modelCoeff$XProt_Contrast > (-1e+06) &
                                                                         RUVnC_TrPPepTrPProtStrain_G2$modelCoeff$XPep_Contrast > (-1e+06)])
    mat <- data.frame("RUV(1)G1" = RUVnC_TrPPepTrPProtStrain_G1$modelCoeff[peps, cName1],
                      "RUV(2)G1" = RUVnC_TrPPepTrPProt_CON_Strain_G1$modelCoeff[peps, cName2],
                      "RUV(3)G1" = RUV_TrPPepTrPProt_CON_Strain_G1$modelCoeff[peps, cName2],
                      "RUV(1)G2" = RUVnC_TrPPepTrPProtStrain_G2$modelCoeff[peps, cName1],
                      "RUV(2)G2" = RUVnC_TrPPepTrPProt_CON_Strain_G2$modelCoeff[peps, cName2],
                      "RUV(3)G2" = RUV_TrPPepTrPProt_CON_Strain_G2$modelCoeff[peps, cName2])
    cor <- cor(mat, use = "complete.obs")
    corrplot(cor, type = "upper", method = "color", addCoef.col = 'black', cl.pos = "n",
             tl.col = "black", tl.srt = 45, is.corr = F, col = col(200))
}

#### Extended Data Figure 3c
biologicalCorPlot("StrainRM11-1a_Contrast", "Strain_Contrast")
biologicalCorPlot("XPep_Contrast", "XPep_RUV")
biologicalCorPlot("XProt_Contrast", "XProt_RUV", colorRampPalette(c("#0306ad", "#ffffff", "#c94202")))

                  
