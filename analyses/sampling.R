# load packages
library(iNEXT)
library(ggplot2)

# load dataset
dataset <- read.delim("../table/dataset.tab", header = TRUE, stringsAsFactors = FALSE)

# Assign a one in individualCount for Farlowella mariaelenae, represented by just one lot
dataset[which(dataset$scientificName == "Farlowella mariaelenae"), "individualCount"] <- 1

abundances <- aggregate(individualCount ~ scientificName, data = dataset, FUN = sum, na.rm = TRUE)

### Assessing sampling efficiency

# iNEXT accumulation curve
specAccum <- iNEXT(x = abundances$individualCount, q = 0, datatype = "abundance")

# iNEXT plotting
specAccumPlot <- ggiNEXT(specAccum, grey = TRUE, type = 1) +
    ggtitle("Sampling in the Cusiana basin") +
    xlab("Number of individuals") +
    ylab("Richness") +
    guides(linetype = guide_legend(title = "Method"),
           colour = guide_legend(title = "none"),
           fill = guide_legend(title = "none"),
           shape = guide_legend(title = "none")) +
    geom_text(x = 35000, y = 200,
              label = paste("S.obs = ", specAccum$AsyEst[1, "Observed"], sep = ""),
              cex = 8, fontface = 3, hjust = 0) + 
    geom_text(x = 35000, y = 175,
              label = paste("S.est = ", round(specAccum$AsyEst[1, "Estimator"], 1), sep = ""),
              cex = 8, fontface = 3, hjust = 0) +
    theme(legend.position = "none")
specAccumPlot

# save plot
ggsave("accumulation.tiff", specAccumPlot)

### Histogram of elevations

# create a vector of breaks for the histogram 250 masl wide
br <- seq(floor(min(dataset$verbatimElevation)/250)*250, ceiling(max(dataset$verbatimElevation)/250)*250, by = 250)

# create a histogram object to be modified further with logarithmic y-label
histo <- hist(dataset$verbatimElevation, breaks = br, plot = FALSE)
logHisto <- histo
logHisto$counts <- log10(logHisto$counts)

# plot both histograms
tiff("elevationHistogram.tiff", width = 1500, height = 1500/2, pointsize = 24)
# save old parameters
oldPar <- par()
# set new outer margin parameter
par(oma = c(0, 0, 3, 0))

# split screen in one row and two columns for side-by-side plots
split.screen(c(1, 2))
# title
title("Distribution of sampling records\nalong the elevational gradient", outer = TRUE)
# histogram of absolute frequency
screen(1)
plot(histo, xlab = "Elevation (m asl)", ylab = "Absolute frequency",
     main = NULL, 
     col = "grey50")
# plot an asterisk above the bars with lowest abundance
text(x = c((2250+2000)/2, (2750+3000)/2),
     y = c(50, 50), labels = "*", cex = 3)

# histogram of log(absolute frequency)
screen(2)
plot(logHisto, ylim = c(0, max(logHisto$counts)),
     xlab = "Elevation (m asl)", ylab = "log (Frequency)",
     main = NULL, 
     col = "grey50")
# plot an asterisk above the bars with lowest abundance
text(x = c((2250+2000)/2, (2750+3000)/2), y = c(0.1, 0.1),
     labels = "*", cex = 3)
# close device
dev.off()
