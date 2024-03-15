library(limma)
library(LiPAnalyzeR)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(patchwork)

# Loading data
QuantityList <- readRDS("../Documents/LiPAnalyzer/Data/CSF_PD/SpectroOut/290224/SpectroMat_HG_FT_max20NAsPerCohort_BatchCorrection.rds")
QuantityListLiPonly <- readRDS("../Documents/LiPAnalyzer/Data/CSF_PD/SpectroOut/290224/SpectroMatLiP_HG_FTHT_max20NAsPerCohort_BatchCorrection.rds")
QuantityListHTonly <- readRDS("~/Documents/LiPAnalyzer/Data/CSF_PD/SpectroOut/290224/SpectroMat_HG_HTonly_matchedTrPViaLiPAnalyzeR_max20NAsPerCohort_BatchCorrection.rds")
annotSample <- readRDS("../Documents/LiPAnalyzer/Data/CSF_PD/SpectroOut/SampleInfo.rds")
annotPepProt <- readRDS("../Documents/LiPAnalyzer/Data/CSF_PD/SpectroOut/290224/PepProtInfo.rds")

## Reducing sample info file to healthy donors
annotSample <- annotSample[colnames(QuantityList$LiPPep), ]

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
plotData <- data.frame(TrPProt = CorCols(row.names(QuantityList$LiPPep), 
                                         QuantityList$LiPPep, 
                                         QuantityList$TrPProt),
                       TrPPep = CorCols(row.names(QuantityList$LiPPep), 
                                        QuantityList$LiPPep, 
                                        QuantityList$TrPPep))
#### Figure 4b
plotScatter(plotData)

plotData <- melt(plotData)
#### Figure 4a
plotDensity(plotData)

### HT peps with matched FT peps
plotData <- data.frame(TrPProt = CorCols(row.names(QuantityListHTonly$LiPPep),
                                         QuantityListHTonly$LiPPep,
                                         QuantityListHTonly$TrPProt),
                       TrPPep = CorCols(row.names(QuantityListHTonly$LiPPep),
                                        QuantityListHTonly$LiPPep,
                                        QuantityListHTonly$TrPPep))
plotData <- melt(plotData)
#### Extended Data Figure 5b
plotDensity(plotData)

### LiP Prot with  TrP Prot
TrPProt_PerProt <- do.call(rbind, lapply(split(as.data.frame(QuantityList$TrPProt), 
                                               annotPepProt[row.names(QuantityList$TrPProt), "Protein"], drop = F), \(x) x[1,]))
LiPProt_PerProt <- do.call(rbind, lapply(split(as.data.frame(QuantityListLiPonly$LiPProt), 
                                               annotPepProt[row.names(QuantityListLiPonly$LiPProt), "Protein"], drop = F), \(x) x[1,]))
isect <- intersect(row.names(TrPProt_PerProt), row.names(LiPProt_PerProt))

CorPlot <- CorCols(intersect(row.names(TrPProt_PerProt), row.names(LiPProt_PerProt)),
                             TrPProt_PerProt, LiPProt_PerProt)
#### Extended Data Figure 6a
ggplot(mapping =  aes(x = CorPlot)) +
    geom_density(linewidth = 0.7, color = "#0057b7", fill = "#0057b7", alpha = 0.8) +
    ylab("Density") +
    xlab("Correlation coefficient") +
    geom_vline(xintercept = median(CorPlot), color = "#0057b7", linewidth = 1) +
    xlim(-1, 1) +
    scale_alpha_manual(values = c(0.8, 0.4)) +
    theme_bw(base_size = 11))




## Variation of HT and FT peptides
### Creating HT and FT specific LiPPep df
LiP_noNA <- QuantityListLiPonly$LiPPep[apply(QuantityListLiPonly$LiPPep, 1, \(x){!is.na(sum(x))}),]
LiPHT <- QuantityListLiPonly$LiPPep[annotPepProt[row.names(QuantityListLiPonly$LiPPep), "isTryptic"] == "Specific-C"|
                                    annotPepProt[row.names(QuantityListLiPonly$LiPPep), "isTryptic"] == "Specific-N", ]
LiPFT <- QuantityListLiPonly$LiPPep[annotPepProt[row.names(QuantityListLiPonly$LiPPep), "isTryptic"] == "Specific", ]

LiPHT <- LiP_noNA[annotPepProt[row.names(LiP_noNA), "isTryptic"] == "Specific-C"|
                      annotPepProt[row.names(LiP_noNA), "isTryptic"] == "Specific-N", ]
LiPFT <- LiP_noNA[annotPepProt[row.names(LiP_noNA), "isTryptic"] == "Specific", ]

sd_LiPFTHT <- data.frame(SD = c(apply(LiPFT, 1, sd), apply(LiPHT, 1, sd)),
                         isTryptic = c(rep("FT", nrow(LiPFT)),
                                       rep("HT", nrow(LiPHT))))

#### Extended Data Figure 5h
ggplot(sd_LiPFTHT, aes(y = SD, x = isTryptic)) +
    geom_violin(alpha = 0.5, fill = "#a7a7a7") +
    stat_summary(fun = "mean",
                 geom = "crossbar", 
                 width = 0.3,
                 colour = "black") +
    xlab("") +
    ylab("Standard deviation per LiP peptide") +
    theme_bw(base_size = 15)



# Example of alternative splicing
## Preprocessing data
### Filtering for peptides from protein of interest
myPeps <- intersect(row.names(QuantityList$LiPPep),
                    row.names(annotPepProt[grepl("P00739", annotPepProt$AllProtein), ]))

### Getting annotPepProt data.frame only for peptides from the respective protein
annotmyPeps <- annotPepProt[myPeps, ]
annotmyPeps$start <- as.numeric(annotmyPeps$start) ## removing peptides with different/multiple start positions
annotmyPeps$end <- as.numeric(annotmyPeps$end)
annotmyPeps <- annotmyPeps[!is.na(annotmyPeps$end),]
myPeps <- row.names(annotmyPeps)

## Estimating correlations
correlateRows <- function(rows, df1, df2){
    sapply(rows, function(i){
        cor(as.numeric(df1[i, ]), as.numeric(df2[i, ]), use = "pairwise.complete.obs")
    })
}

CorProt <- data.frame(start = annotmyPeps$start, 
                      end = annotmyPeps$end,
                      LiPPep_TrPPep = correlateRows(myPeps, QuantityList$LiPPep, QuantityList$TrPPep),
                      LiPPep_TrPProt = correlateRows(myPeps, QuantityList$LiPPep, QuantityList$TrPProt),
                      TrPPep_TrPProt = correlateRows(myPeps, QuantityList$TrPPep, QuantityList$TrPProt))

### Woods plot with Y axis being isoforms or peptide specific correlation
AllIso <- unique(unlist(sapply(annotmyPeps$AllProtein, function(x){unlist(strsplit(x, ";"))})))
DataIsoPlot <- do.call(rbind, lapply(as.list(AllIso), function(x){
    inX <- sapply(annotmyPeps$AllProtein, function(y) {x %in% unlist(strsplit(y, ";"))})
    data.frame(start = annotmyPeps$start[inX],
               end = annotmyPeps$end[inX],
               Isoform = x)
}))
p1 <- ggplot(DataIsoPlot, aes(x = start, xend = end, y = Isoform, yend = Isoform, alpha = 0.1, color = Isoform)) +
    geom_segment(linewidth = 4) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(CorProt, aes(x = start, xend = end, y = LiPPep_TrPProt, yend =  LiPPep_TrPProt, alpha = 0.1)) +
    geom_segment(size = 4) +
    geom_hline(yintercept = 0, size = 1) +
    ylab("LiPPep & TrPProt") +
    scale_y_continuous(breaks = seq(-0.25, 1.0, 0.25), limits = c(-0.25, 1)) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

p3 <- ggplot(CorProt, aes(x = start, xend = end, y = TrPPep_TrPProt, yend = TrPPep_TrPProt, alpha = 0.1)) +
    geom_segment(size = 4) +
    geom_hline(yintercept = 0, size = 1) +
    ylab("TrPPep & TrPProt") +
    scale_y_continuous(breaks = seq(-0.25, 1.0, 0.25), limits = c(-0.25, 1)) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

p4 <- ggplot(CorProt, aes(x = start, xend = end, y = LiPPep_TrPPep, yend = LiPPep_TrPPep, alpha = 0.1)) +
    geom_segment(size = 4) +
    geom_hline(yintercept = 0, size = 1) +
    ylab("TrPPep & TrPProt") +
    scale_y_continuous(breaks = seq(-0.25, 1.0, 0.25), limits = c(-0.25, 1)) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none")

#### Figure 4c
p1/p2/p3/p4 + plot_layout(heights = c(2.5, 5, 5, 5))

### Plotting single example peptides
myPepsSort <- row.names(annotmyPeps)[order(annotmyPeps$start)]
col <- setNames(ifelse(QuantityList$TrPPep["LRTEGDGVYTLNDKK",] >16, "norm", "alt"), 
                colnames(QuantityList$LiPPep))
col <- col[order(unname(col))]
QuantityListPlot <- lapply(QuantityList, \(x){x[myPepsSort, names(col)]})

plotPoints <- function(xAxis, yAxis, xName, yName, title, col){
    ggplot(mapping = aes(x = xAxis,
                         y = yAxis,
                         color = col, alpha = col)) +
        geom_point(size = 2) +
        scale_color_manual(values = c("#c94202", "#000000")) +
        scale_alpha_manual(values = c(1, 0.3)) +
        ggtitle(title) +
        xlab(xName) +
        ylab(yName) +
        theme_bw(base_size = 11) +
        theme(legend.position = "none", aspect.ratio = 1)
}

plotList_TrPPepTrPProt <- lapply(as.list(myPepsSort), \(x){
    plotPoints(unname(QuantityListPlot$TrPProt[x, ]),
               unname(QuantityListPlot$TrPPep[x, ]),
               "TrPProt quantity", "TrPPep quantity", x, col)
})

plotList_LiPPepTrPProt <- lapply(as.list(myPepsSort), \(x){
    plotPoints(unname(QuantityListPlot$TrPProt[x, ]),
               unname(QuantityListPlot$LiPPep[x, ]),
               "TrPProt quantity", "LiPPep quantity", x, col)
})

plotList_LiPPepTrPPep <- lapply(as.list(myPepsSort), \(x){
    plotPoints(unname(QuantityListPlot$TrPPep[x, ]),
               unname(QuantityListPlot$LiPPep[x, ]),
               "TrPPep quantity", "LiPPep quantity", x, col)
})

#### Figure 4d
i <- 2
plotList_TrPPepTrPProt[[i]] + plotList_LiPPepTrPProt[[i]] + plotList_LiPPepTrPPep[[i]] + plot_layout(ncol = 3)

#### Figure 4e
i <- 9
plotList_TrPPepTrPProt[[i]] + plotList_LiPPepTrPProt[[i]] + plotList_LiPPepTrPPep[[i]] + plot_layout(ncol = 3)




# Estimating Ratios
rLiPPepTrPPep <- QuantityList$LiPPep- QuantityList$TrPPep
rLiPPepTrPProt <- QuantityList$LiPPep- QuantityList$TrPProt

## Correlations with RatLipPepTrPPep
cRLiPPepTrPPep_TrPPep <- CorCols(row.names(rLiPPepTrPPep),
                                 rLiPPepTrPPep,
                                 QuantityList$TrPPep)

cRLiPPepTrPPep_TrPProt <- CorCols(row.names(rLiPPepTrPPep),
                                  rLiPPepTrPPep,
                                  QuantityList$TrPProt)

## Correlations with RatLipPepTrPProt
cRLiPPepTrPProt_TrPPep <- CorCols(row.names(rLiPPepTrPProt),
                                  rLiPPepTrPProt,
                                  QuantityList$TrPPep)

cRLiPPepTrPProt_TrPProt <- CorCols(row.names(rLiPPepTrPProt),
                                   rLiPPepTrPProt,
                                   QuantityList$TrPProt)

### Plotting correlations
plotData <- melt(data.frame(TrPProt = cRLiPPepTrPProt_TrPProt,
                           TrPPep = cRLiPPepTrPPep_TrPPep))

#### Figure 5b
ggplot(plotData, aes(x = value, color = variable, fill = variable, alpha = variable)) +
    geom_density(linewidth = 0.7) +
    ylab("Density") +
    xlab("Correlation coefficient") +
    xlim(-1, 1) +
    geom_vline(xintercept = median(cRLiPPepTrPProt_TrPProt), color = "#0057b7", linewidth = 1) +
    geom_vline(xintercept = median(cRLiPPepTrPPep_TrPPep), color = "#ffc900", linewidth = 1) +
    scale_color_manual(values = c("#0057b7", "#ffc900")) +
    scale_fill_manual(values = c("#0057b7", "#ffc900")) +
    scale_alpha_manual(values = c(0.8, 0.4)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")


# LiPAnalyzeR
### setting sex for deviation coding
annotSample$Sex <- as.factor(annotSample$Sex)
contrasts(annotSample$Sex) <- contr.sum(2)

## Running only RUV step
### RUV with constraints removing only TrPPep
RUV_TrPPep <- runModel(quantityList = QuantityList,
                       annotS = annotSample,
                       formulaRUV = "Y~XPep", 
                       lowRUV = c(-Inf, 0),
                       upRUV = c(Inf, Inf),
                       formulaContrast = NULL)

### RUV with constraints removing only TrPProt
RUV_TrPProt <- runModel(quantityList = QuantityList,
                        annotS = annotSample,
                        formulaRUV = "Y~XProt", 
                        lowRUV = c(-Inf, 0),
                        upRUV = c(Inf, Inf),
                        formulaContrast = NULL)

### RUV with constraints removing TrPPep & TrPProt
RUV_TrPPepTrPProt <- runModel(quantityList = QuantityList,
                              annotS = annotSample,
                              formulaRUV = "Y~XPep+XProt",
                              formulaContrast = NULL)

### RUV with constraints removing TrPPep & TrPProt in HTonly mode 
RUV_TrPPepTrPProt_HTonly <- runModel(quantityList = QuantityListHTonly,
                                     annotS = annotSample,
                                     formulaRUV = "Y~XPep+XProt",
                                     formulaContrast = NULL)

### RUV with constraints removing LiPProt in LiPonly mode 
RUV_LiPProt <- runModel(quantityList = QuantityListLiPonly,
                        annotS = annotSample,
                        formulaRUV = "Y~XProt", 
                        lowRUV = c(-Inf, 0),
                        upRUV = c(Inf, Inf),
                        formulaContrast = NULL)

### RUV step without constraints removing only TrPPep
RUVnC_TrPPep <- runModel(quantityList = QuantityList,
                         annotS = annotSample,
                         formulaRUV = NULL,
                         formulaContrast = "Y~XPep")

### RUV step without constraints removing only TrPProt
RUVnC_TrPProt <- runModel(quantityList = QuantityList,
                          annotS = annotSample,
                          formulaRUV = NULL,
                          formulaContrast = "Y~XProt")

### RUV step without constraints removing TrPPep & TrPProt
RUVnC_TrPPepTrPProt <- runModel(quantityList = QuantityList,
                                annotS = annotSample,
                                formulaRUV = NULL,
                                formulaContrast = "Y~XPep+XProt")




## Plotting correlations of residuals with TrPPep/TrPProt
### Correlations to residuals of RUV corrected for only TrPPep or TrPProt
peps <- row.names(RUV_TrPPepTrPProt$modelCoeff)[RUV_TrPPepTrPProt$modelCoeff$XProt_RUV > (-1e+06) &
                                                RUV_TrPPepTrPProt$modelCoeff$XPep_RUV > (-1e+06)]
plotData <- melt(data.frame(RUV_TrPProt = CorCols(peps,
                                                  RUV_TrPPep$modelResid, 
                                                  QuantityList$TrPProt) ,
                            RUV_TrPPep = CorCols(peps, 
                                                 RUV_TrPProt$modelResid,
                                                 QuantityList$TrPPep)))
#### Figure 5d
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
plotData <-  melt(data.frame(RUV_TrPProt = CorCols(row.names(RUV_TrPPepTrPProt$modelResid),
                                                          RUV_TrPPepTrPProt$modelResid, 
                                                          QuantityList$TrPProt),
                                    RUV_TrPPep = CorCols(row.names(RUV_TrPPepTrPProt$modelResid), 
                                                         RUV_TrPPepTrPProt$modelResid,
                                                         QuantityList$TrPPep)))
#### Figure 5f
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
plotData <- RUVnC_TrPPepTrPProt$modelCoeff
P1 <- ggplot(mapping = aes(x = plotData$XProt_Contrast,
                           y = plotData$XPep_Contrast, 
                           color = plotData$XProt_Contrast<0|
                               plotData$XPep_Contrast<0)) +
    geom_point(size = 1, alpha = 0.1) +
    scale_color_manual(values = c("#000000", "#c94202")) +
    xlim(-5, 6) +
    ylim(-5, 6) +
    ylab("Coef RUV - TrPPep") +
    xlab("Coef RUV - TrPProt") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
P2 <- ggplot(mapping = aes(x = plotData$XPep_Contrast)) +
    geom_histogram(breaks = seq(-5, 6, 0.2)) +
    xlim(-5, 6) +
    coord_flip() +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = -90))
P3 <- ggplot(mapping = aes(x = plotData$XProt_Contrast)) +
    geom_histogram(breaks = seq(-5, 6, 0.2)) +
    xlim(-5, 6) +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
          #### Extended Data Figure 2a
P3 + plot_spacer() + P1 + P2 +
    plot_layout(nrow = 2, ncol = 2, widths = c(2, 1), heights = c(1, 2))

### Coefficients TrPPep & TrPProt from RUV with constraints
plotData <- RUV_TrPPepTrPProt$modelCoeff

P1 <- ggplot(plotData, aes(x = XProt_RUV,
                           y = XPep_RUV)) +
    geom_point(size = 1, alpha = 0.1) +
    xlim(0, 4) +
    ylim(0, 4) +
    ylab("Coef RUV - TrPPep") +
    xlab("Coef RUV - TrPProt") +
    theme_bw(base_size = 11)
P2 <- ggplot(plotData, aes(x = XPep_RUV)) +
    geom_histogram(breaks = seq(0, 4, 0.1)) +
    xlim(0, 4) +
    coord_flip() +
    ylab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = -90))
P3 <- ggplot(plotData, aes(x = XProt_RUV)) +
    geom_histogram(breaks = seq(0, 4, 0.1)) +
    xlim(0, 5) +
    xlab("Count") +
    theme_bw(base_size = 11) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
#### Extended Data Figure 2b
P3 + plot_spacer() + P1 + P2 +
    plot_layout(nrow = 2, ncol = 2, widths = c(2, 1), heights = c(1, 2))

