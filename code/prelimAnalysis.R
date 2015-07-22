# File: prelimAnalysis.R
# Authour: Shyam Gopalakrishnan
# Date: 15th July 2015
# This script reads in the hapmap style
# genotype file and generates some simple
# graphs to parse through the information
# to get some sense of that the data looks
# like.

library(RColorBrewer)
colors = rep("", 8)
colors[1:3] = c("burlywood", "darkorange", "firebrick")
colors[4] = "darkgray" # LOW
colors[5] = "forestgreen" # OUTGROUP
colors[6:8]  = c("cadetblue", "dodgerblue", "navyblue")

# Read genotypes
genotypes <- read.table("/home/shyam/projects/rainforest/data/GBS_Erythrophleum_48X2Samples/TASSEL2/hapMap/HapMap.hmp.txt", header=TRUE, as.is=TRUE)
aux.info  <- genotypes[,1:11]
genotypes <- genotypes[,-c(1:12)]
names(genotypes) <- sub("_.+", "", names(genotypes))

# Read estimated assignments and positions
assigns   <- read.table("/home/shyam/projects/rainforest/data/Erythrophleum_Pop_ID.csv", header=T, as.is=TRUE, sep=',')
clades    <- sapply(names(genotypes),
                    function(x) {
                      if (x %in% assigns$POOLED_ID_POP   ) {
                        temp = assigns[which(assigns$POOLED_ID_POP == x), "GBS_CLADE"]
                      } else {
                        temp = NA
                      }
                      return (temp)
                    })

                        
# Compute missingness per sample and group
total.sites   <- nrow(genotypes)
ind.miss      <- apply(genotypes, 2, function(x) { sum(x=="N") } )
clade.miss    <- tapply(ind.miss, clades, mean)

# Plot missingness
require(lattice)
pdf(file="../analysis/missing.pdf")
ord <- order(ind.miss)
barplot(ind.miss[ord]/total.sites, col=colors[as.numeric(as.factor(clades))[ord]], las=3, axes=F, names.arg=NA, space=0, ylim=c(0,1),
        main='Per sample missingness', ylab='Proportion of missing genotypes')
axis(2, at=seq(0,1,0.2))
axis(1, at=1:length(ind.miss), cex.axis=0.6, las=3, labels=names(ind.miss)[ord])
legend("topleft", cex=0.8, pch=15, col=colors, legend=levels(as.factor(clades)))
bwplot((ind.miss/total.sites)~clades, ylab="Genotype missingness proportion", main="Missingness by clade", ylim=c(0.5,1))
plot(as.numeric(as.factor(clades)), (ind.miss/total.sites), axes=F, frame.plot=T, xlab="clade", ylab="missing proportion", main="missingness by clade", col=colors[as.numeric(as.factor(clades))], pch=19, ylim=c(0.5,1))
points(1:8, clade.miss/total.sites, pch='-', col=colors, cex=3)
axis(2)
axis(1, at=1:8, labels=names(clade.miss), cex.axis=0.7, tick=F)
dev.off()

# Compute heterozygosity per sample
homcodes       <- c('A', 'C', 'G', 'T')
total.sites    <- apply(genotypes, 2, function(x) { sum(x!="N") } )
ind.het        <- apply(genotypes, 2, function(x) { sum(x%in%homcodes)})
ind.het        <- (total.sites - ind.het)/total.sites
clade.het      <- tapply(ind.het, clades, mean)

pdf(file="../analysis/het.pdf")
ord  <- order(ind.het)
barplot(ind.het[ord], col=colors[as.numeric(as.factor(clades))[ord]], las=3, axes=F, names.arg=NA, space=0, ylim=c(0,1),
        main='Per sample heterozygosity', ylab='Heterozygosity (at variant sites)')
axis(2, at=seq(0,1,0.2))
axis(1, at=1:length(ind.het), cex.axis=0.6, las=3, labels=names(ind.het)[ord])
legend("topleft", cex=0.8, pch=15, col=colors, legend=levels(as.factor(clades)))
bwplot(ind.het~clades, ylab="heterozygosity", main="heterozygosity by clade")
plot(as.numeric(as.factor(clades)), ind.het, axes=F, frame.plot=T, xlab="clade", ylab="heterozygosity", main="heterozygosity by clade", col=colors[as.numeric(as.factor(clades))], pch=19, ylim=c(0,1))
points(1:8, clade.het, pch='-', col=colors, cex=3)
axis(2)
axis(1, at=1:8, labels=names(clade.miss), cex.axis=0.7, tick=F)
dev.off()
