source("https://bioconductor.org/biocLite.R")
install.packages("ggplot2")
biocLite(c("limma", "edgeR", "org.Mm.eg.db"))

library(limma)
library(edgeR)
library(org.Mm.eg.db)
library(ggplot2)

# Set working directory to the location with the raw data file
# setwd("G:\\My Drive\Genomics\lab")
# Load in the raw data through the read.csv function
rawData <- read.csv("seqSample.csv")
# Take a look to see what the raw data object is like
dim(rawData)
head(rawData)
str(rawData)

# Set the gene IDs (1st column) as the rownames and remove the gene column
geneVec <- rawData[,1] # save the first column as a vector of gene ids
rownames(rawData) <- rawData$X # make the first column the row names
rawData <- rawData[,-1] # remove the first column (X) from rawData now that it's saved as the rownames

# Set up the group identifier
group <- factor(c(rep(1:4, each = 3)),
                labels = c("WT_n", "WT_6h", "KO_n", "KO_6h"))

# Create a DGE object
# This object is a list-type object where we'll store many pieces of information
x <- DGEList(counts = rawData, genes = geneVec)
x$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(x),
                         keytype = "ENTREZID", column = "SYMBOL")

# Examine the distribution of the CPMs for different transcripts
# CPM normalized to transcript abundance
mean(rowSums(cpm(x)))
hist(rowSums(cpm(x)),
     xlim = range(0, 1000),
     breaks = 10000)

# Only keep transcripts above a certain abundance
keep <- rowSums( cpm(x) >1) >= 2
table(keep)
x <- x[keep, , keep.lib.sizes=FALSE]

# Calculate normalization factors and validate graphically
x <- calcNormFactors(x)
x$samples
plotMD(cpm(x, log = T), column = 1)
abline(h = 0, col = 'red', lty = 2, lwd = 2)

# Examine global grouping of samples
pch <- c(0,1,2,3)
colors <- rep(c("red", "blue", "green", "pink"), 3)
scatPlot <- plotMDS(x, col = colors[group], pch = pch[group])
# Lets take a closer look at the scatplot data object

# One of the benefits of R is graphics generation
# We can use a package called ggplot to create a nicer version of the previous graph
scatMat <- data.frame(cbind(as.numeric(scatPlot$x), as.numeric(scatY <- scatPlot$y)))
colnames(scatMat) <- c("x", "y")
scatMat$group <- group
ggplot(data = scatMat, aes(x = x, y = y, group = group)) +
  geom_point(aes(color = group, size = 2)) +
  scale_size(guide = 'none')

# Set up design matrix and model
designMat <- model.matrix(~0+group)
colnames(designMat) <- levels(group)

# Estimate dispersion parameters for your data
x <- estimateCommonDisp(x, designMat)
x <- estimateGLMTrendedDisp(x, designMat)
x <- estimateGLMTagwiseDisp(x, designMat)

fit <- glmQLFit(x, designMat, robust = T)
plotQLDisp(fit)

# Make comparisons to test for DE
comparison <- makeContrasts(WT_n - WT_6h, levels = designMat)
qlf <- glmQLFTest(fit, contrast = comparison)

# Check the top DE genes
topTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)

# There are a lot of DE genes, better to filter based on fold change cutoff
tr <- glmTreat(fit, contrast = comparison, lfc = log2(1.2), null = "worst.case")
summary(decideTests(tr))
topTags(tr)
plotMD(tr)

# Testing using an ANOVA framework --> does not increase power, but makes analysis more efficient
anvaCon <- makeContrasts(
  WTn_WT6h = WT_n - WT_6h,
  KOn_KO6h = KO_n - KO_6h,
  WTn_KOn = WT_n - KO_n,
  WT6h_KO6h = WT_6h - KO_6h,
  levels = designMat
)

anva <- glmQLFTest(fit, contrast = anvaCon)
topTags(anva)

# Gene ontology (pathway) analysis --> what everyone in the field wants to do
GOcon <- makeContrasts(WT_n - WT_6h, levels = designMat)
qlf <- glmQLFTest(fit, contrast = GOcon)
go <- goana(qlf, species = "Mm")
topGO(go, number = 30)

# kegg version
kegg <- kegga(qlf, species = "Mm")
topKEGG(kegg, number = 30)

# a lot of these functions are incorporated into bioconductor packages

### volcano plots
library(plotly)
# check website for the updated script with this section included
