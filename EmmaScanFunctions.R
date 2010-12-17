# Author: Keith Sheppard <keith.sheppard@jax.org>
#
# This script contains some convenience functions for scanning genotype and
# phenotype data with the emma package <http://mouse.cs.ucla.edu/emma/>
# inputs are:
#   - MPD Phenotypes in the "Individual" data point format.
#     See: <http://www.jax.org/phenome>
#   - SNP CSV data files where the 1st row is a header (strain names) and remaining
#     rows contain SNPs (a,g,c,t) case insensitive. Note: any position or
#     annotation columns must be removed before using this script since only
#     SNP columns are allowed
# Make sure that the strain order used for the phenotype files and genotype files
# is the same (it is OK for the phenotype file to have multiple rows for a
# single strain, but they must occur in the correct order)
#
# Here is a simple example of how you might invoke this script:
# source("EmmaScanFunctions.R")
# emmaScan(
#     genoFiles=dir("geno-data", full.names=T),
#     mpdPhenoFile="mpd-formatted-phenotypes.txt",
#     resultsDir="emma-scan-results")

library(emma)

# reads in the given CSV file with a header row for strains and a,g,c,t (case-insensitive)
# for all values and return a matrix of 1,0 values where 1 means that the leftmost
# strain is matched and 0 means there is a difference
readSnpsAsBinaryMatrix <- function(csvSNPFile)
{
    snpMat <- read.csv(csvSNPFile, colClasses="character")
    snpMatrixToBinary(snpMat)
}

readSNPs <- function(
        csvSNPFile,
        chromosomeColumn,
        basePairPositionColumn,
        preserveColumns,
        aAlleleColumn,
        bAlleleColumn,
        firstGenotypeColumn,
        lastGenotypeColumn,
        strainNames)
{
    csvSNPFile <- file(csvSNPFile, open = "rt")
    
    # take a look at the header and make sure the strains match up perfectly
    fileHeader <- read.csv(csvSNPFile, colClasses="character", nrows = 1)[1, ]
    if(is.null(lastGenotypeColumn))
    {
        lastGenotypeColumn <- length(fileHeader)
    }
    
    snpData <- list()
    snpData$allStrainNames <- fileHeader[firstGenotypeColumn : lastGenotypeColumn]
    snpData$keptStrainNames <- sort(snpData$allStrainNames[snpData$allStrainNames %in% strainNames])
    
    if(snpData$keptStrainNames < 4)
    {((( here we are
        stop("There are not enough matching strain names available for testing.",
             "We only have: ", paste())
    }
    
    snpData
}

# convert the given case insensitive g,a,c,t character matrix into a 0/1 integer matrix
snpMatrixToBinary <- function(snpMat)
{
    # the '1' strains match the 1st strain and the others will be '0'
    binMat <- t(apply(snpMat, 1, function(snpRow) {as.integer(toupper(snpRow) == toupper(snpRow[1]))}))

    # restore column names and return
    colnames(binMat) <- colnames(snpMat)
    binMat
}

# reads in an "MPD formatted" phenotype file with individual phenotype values
readMpdIndividualPhenotypes <- function(mpdPhenoFile)
{
    phenoMat <- read.delim(mpdPhenoFile, colClasses="character")
    phenoMat
}

# create the matrix that maps individual values to strains. With the returned
# matrix the rows represent the individuals and the columns represent what
# strain those individuals belong to.
mpdPhenosToStrainIncidenceMatrix <- function(mpdPhenos)
{
    phenoRle <- rle(mpdPhenos[["strain"]])
    strainIndividualCounts <- phenoRle$lengths
    numStrains <- length(strainIndividualCounts)
    
    incidenceMatrix <- NULL
    for(currStrainIndex in 1:numStrains)
    {
        currStrainIndCount <- strainIndividualCounts[currStrainIndex]
        
        currRow <- rep.int(0, numStrains)
        currRow[currStrainIndex] <- 1
        
        for(i in 1:currStrainIndCount)
        {
            incidenceMatrix <- rbind(incidenceMatrix, currRow)
        }
    }
    
    incidenceMatrix
}

# pull out only the phenotype values leaving all of the other MPD data behind
mpdPhenosToPhenoOnly <- function(mpdPhenos)
{
    rbind(NULL, as.numeric(mpdPhenos[["value"]]))
}

emmaScan <- function(
        mpdPhenoFile,
        genoFiles,
        resultsDir,
        chromosomeColumn = 4,
        basePairPositionColumn = 5,
        preserveColumns = 1 : 5,
        aAlleleColumn = 2,
        bAlleleColumn = 3,
        firstGenotypeColumn,
        lastGenotypeColumn = NULL,
        verbose = TRUE)
{
    # read in SNPs and calculate the kinship matrix
    #
    # snpMatrixList - a list of SNP matrices (one per genoFile)
    # allSnpsMatrix - combines all SNPs into a single matrix for calculating kinship
    snpMatrixList <- NULL
    allSnpsMatrix <- NULL
    for(genoFile in genoFiles)
    {
        print(paste("reading in genotypes from:", genoFile))
        
        currSnpMatrix <- readSnpsAsBinaryMatrix(genoFile)
        snpMatrixList <- append(snpMatrixList, list(currSnpMatrix))
        allSnpsMatrix <- rbind(allSnpsMatrix, currSnpMatrix)
    }
    
    print("calculating kinship")
    kinshipMat <- emma.kinship(allSnpsMatrix)
    
    # read in phenos along with an incidence matrix
    print("reading phenotypes")
    mpdPhenos <- readMpdIndividualPhenotypes(mpdPhenoFile)
    incMat <- mpdPhenosToStrainIncidenceMatrix(mpdPhenos)
    phenosOnly <- mpdPhenosToPhenoOnly(mpdPhenos)
    
    chr <- 0
    for(snpMatrix in snpMatrixList)
    {
        chr <- chr + 1
        
        print(paste("performing t-test for:", genoFiles[chr]))
        currScan <- emma.REML.t(ys=phenosOnly, xs=snpMatrix, K=kinshipMat, Z=incMat)
        
        print("exporting scan results")
        scanResultsFile <- paste(resultsDir, "/", "scan", chr, ".csv", sep="")
        
        # write some header info
        fileHandle <- file(scanResultsFile, "a")
        writeLines(
            c(paste("genotype data source:", genoFiles[chr]), paste("phenotype data source:", mpdPhenoFile)),
            fileHandle)
        close(fileHandle)
        
        # write the table
        write.csv(currScan, file=scanResultsFile, row.names=F, append=T)
        
        # create the scan image
        scanImageFile <- paste(resultsDir, "/", "scan", chr, ".png", sep="")
        png(filename=scanImageFile, width=800, height=300)
        plot(-log10(currScan$ps), type="l")
        dev.off()
    }
}
