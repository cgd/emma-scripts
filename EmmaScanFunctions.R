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
# Here is a simple example of how you might invoke the functions in this script:
# source("EmmaScanFunctions.R")
# scanResults <- emmaScan(
#   mpdPhenoFile = "mpd-formatted-phenotypes.txt",
#   genoFiles = dir("geno-data", full.names = T))

library(emma)

emmaScan <- function(
        mpdPhenoFile,
        genoFiles,
        chromosomeColumn = 4,
        basePairPositionColumn = 5,
        preserveColumns = 1 : 5,
        aAlleleColumn = 2,
        bAlleleColumn = 3,
        firstGenotypeColumn = 6,
        lastGenotypeColumn = NULL,
        sex = c("both", "male", "female"),
        cluster = NULL,
        verbose = TRUE,
        phenoAggregateFun = mean)
{
    sex <- match.arg(sex)
    
    # read in pheno data
    if(verbose)
    {
        cat("reading phenotypes\n")
    }
    mpdPhenos <- .readMpdIndividualPhenotypes(mpdPhenoFile)
    if(verbose)
    {
        cat("read ", nrow(mpdPhenos), " phenotypes\n", sep = "")
    }
    unfilteredPhenoStrains <- sort(unique(mpdPhenos[["strain"]]))
    if(sex != "both")
    {
        if(sex != "male" && sex != "female")
        {
            stop("bad value given for sex parameter")
        }
        
        if("sex" %in% names(mpdPhenos))
        {
            phenoSex <- toupper(mpdPhenos[["sex"]])
            if(sex == "male")
            {
                sexFilter <- phenoSex == "MALE" | phenoSex == "M"
                if(verbose)
                {
                    cat("removing ", sum(!sexFilter), " female phenotypes\n")
                }
            }
            else
            {
                sexFilter <- phenoSex == "FEMALE" | phenoSex == "F"
                if(verbose)
                {
                    cat("removing ", sum(!sexFilter), " male phenotypes\n")
                }
            }
            
            mpdPhenos <- mpdPhenos[sexFilter, ]
        }
        else
        {
            stop("failed to find \"sex\" column in phenotype file")
        }
    }
    
    if(!is.null(phenoAggregateFun))
    {
        aggPhenos <- aggregate(mpdPhenos[["value"]], by = list(mpdPhenos[["strain"]]), FUN = phenoAggregateFun)
        mpdPhenos <- data.frame(strain = aggPhenos[[1]], value = aggPhenos[[2]], stringsAsFactors = FALSE)
    }
    mpdPhenos <- mpdPhenos[order(mpdPhenos[["strain"]]), ]
    
    phenoStrains <- sort(unique(mpdPhenos[["strain"]]))
    
    if(verbose && !setequal(unfilteredPhenoStrains, phenoStrains))
    {
        cat("The following strains were removed during the sex filter phase: ",
                paste(setdiff(unfilteredPhenoStrains, phenoStrains), collapse = ", "),
                "\n")
    }
    
    # read in SNPs and calculate the kinship matrix
    #
    # snpMatrixList - a list of SNP matrices (one per genoFile)
    # allSnpsMatrix - combines all SNPs into a single matrix for calculating kinship
    #snpMatrixList <- NULL
    #allSnpsMatrix <- NULL
    snpDataList <- list()
    allSnpData <- NULL
    for(genoFile in genoFiles)
    {
        if(verbose)
        {
            cat("reading in genotypes from: ", genoFile, "\n", sep = "")
        }
        currSNPData <- .readSNPData(
                genoFile,
                chromosomeColumn,
                basePairPositionColumn,
                preserveColumns,
                aAlleleColumn,
                bAlleleColumn,
                firstGenotypeColumn,
                lastGenotypeColumn,
                phenoStrains)
        snpDataList[[genoFile]] <- currSNPData
        allSnpData <- .bindSnpData(allSnpData, currSNPData)
    }
    
    # subset phenotypes based on strains available in genotypes
    if(verbose)
    {
        removedGenoStrains <- setdiff(allSnpData$allStrainNames, phenoStrains)
        if(length(removedGenoStrains) >= 1)
        {
            cat("Removed ", length(removedGenoStrains), " strains from the genotype data ",
                    "because they were missing from the phenotype data: ",
                    paste(removedGenoStrains, collapse = ", "), "\n", sep = "")
        }
        
        removedPhenoStrains <- setdiff(phenoStrains, allSnpData$allStrainNames)
        if(length(removedPhenoStrains) >= 1)
        {
            cat("Removed ", length(removedPhenoStrains), " strains from the phenotype data ",
                    "because they were missing from the genotype data: ",
                    paste(removedPhenoStrains, collapse = ", "), "\n", sep = "")
        }
        
        cat("The scan will be performed using the following ",
                length(allSnpData$keptStrainNames), " strains: ",
                paste(allSnpData$keptStrainNames, collapse = ", "), "\n", sep = "")
    }
    mpdPhenos <- mpdPhenos[mpdPhenos[["strain"]] %in% allSnpData$keptStrainNames, ]
    
    if(verbose)
    {
        cat("calculating kinship\n")
    }
    kinshipMat <- emma.kinship(allSnpData$snpCalls)
    
    incMat <- .mpdPhenosToStrainIncidenceMatrix(mpdPhenos)
    phenosOnly <- .mpdPhenosToPhenoOnly(mpdPhenos)
    
    if(length(cluster) >= 2)
    {
        if(verbose)
        {
            cat("running EMMA scan using ", length(cluster), " cluster nodes\n", sep = "")
        }
        chunkedCalls <- .chunkRows(allSnpData$snpCalls, length(cluster))
        scanForApply <- function(snpMat, phenos, kin, inc)
        {
            .listToDataFrame(emma.REML.t(ys = phenos, xs = snpMat, K = kin, Z = inc))
        }
        chunkResults <- parLapply(
                cluster,
                chunkedCalls,
                scanForApply,
                phenosOnly,
                kinshipMat,
                incMat)
        results <- NULL
        for(res in chunkResults)
        {
            results <- rbind(results, res)
        }
    }
    else
    {
        if(verbose)
        {
            cat("running EMMA scan\n", sep = "")
        }
        results <- .listToDataFrame(emma.REML.t(
                        ys = phenosOnly,
                        xs = allSnpData$snpCalls,
                        K = kinshipMat,
                        Z = incMat))
    }
    
    if(length(preserveColumns) >= 1)
    {
        results <- cbind(allSnpData$preservedColumns, results)
    }
    results
}

.readSNPData <- function(
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
    fileHeader <- read.csv(csvSNPFile, colClasses="character", header = FALSE, nrows = 1)[1, ]
    if(is.null(lastGenotypeColumn))
    {
        lastGenotypeColumn <- length(fileHeader)
    }
    
    snpData <- list()
    snpData$allStrainNames <- fileHeader[firstGenotypeColumn : lastGenotypeColumn]
    snpData$keptStrainNames <- sort(snpData$allStrainNames[snpData$allStrainNames %in% strainNames])
    
    snpColsToKeep <- match(snpData$keptStrainNames, snpData$allStrainNames)
    snpColsToKeep <- snpColsToKeep + (firstGenotypeColumn - 1)
    
    fileData <- read.csv(csvSNPFile, colClasses = "character", header = FALSE)
    close(csvSNPFile)
    
    snpData$snpCalls <- as.matrix(fileData[ , snpColsToKeep, drop = FALSE])
    snpData$snpCalls <- .snpCallsToNumeric(
            fileData[ , aAlleleColumn],
            fileData[ , bAlleleColumn],
            snpData$snpCalls)
    
    snpData$preservedColumnNames <- fileHeader[preserveColumns]
    snpData$preservedColumns <- fileData[ , preserveColumns, drop = FALSE]
    names(snpData$preservedColumns) <- snpData$preservedColumnNames
    
    snpData$position <- data.frame(
            chromosome = fileData[ , chromosomeColumn],
            bpPosition = as.numeric(fileData[ , basePairPositionColumn]),
            stringsAsFactors = FALSE)
    
    snpData
}

# concatenate SNP data
.bindSnpData <- function(snpData1, snpData2)
{
    if(is.null(snpData1))
    {
        snpData2
    }
    else if(is.null(snpData2))
    {
        snpData1
    }
    else
    {
        snpData1$allStrainNames <- union(snpData1$allStrainNames, snpData$allStrainNames)
        
        if((!setequal(snpData1$keptStrainNames, snpData2$keptStrainNames)) &&
           (!setequal(snpData1$preservedColumnNames, snpData2$preservedColumnNames)))
        {
            stop("column missmatch between geno files")
        }
        else
        {
            snpData1$snpCalls <- rbind(snpData1$snpCalls, snpData2$snpCalls)
            snpData1$position <- rbind(snpData1$position, snpData2$position)
            snpData1$preservedColumns <- rbind(snpData1$preservedColumns, snpData2$preservedColumns)
            
            snpData1
        }
    }
}

.snpCallsToNumeric <- function(aAlleles, bAlleles, snpCalls)
{
    aAlleles <- toupper(aAlleles)
    bAlleles <- toupper(bAlleles)
    snpCalls <- toupper(snpCalls)
    
    numMat <- matrix(NaN, nrow = nrow(snpCalls), ncol = ncol(snpCalls))
    numMat[snpCalls == aAlleles] <- 1.0
    numMat[snpCalls == bAlleles] <- 0.0
    numMat[snpCalls == "H"]      <- 0.5
    
    numMat
}

# reads in an "MPD formatted" phenotype file with individual phenotype values
.readMpdIndividualPhenotypes <- function(mpdPhenoFile)
{
    phenoData <- read.delim(mpdPhenoFile, colClasses = "character")
    phenoData[["value"]] <- as.numeric(phenoData[["value"]])
    phenoData
}

# create the matrix that maps individual values to strains. With the returned
# matrix the rows represent the individuals and the columns represent what
# strain those individuals belong to.
.mpdPhenosToStrainIncidenceMatrix <- function(mpdPhenos)
{
    phenoRle <- rle(mpdPhenos[["strain"]])
    strainIndividualCounts <- phenoRle$lengths
    
    if(all(strainIndividualCounts) == 1)
    {
        NULL
    }
    else
    {
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
}

# pull out only the phenotype values leaving all of the other MPD data behind
.mpdPhenosToPhenoOnly <- function(mpdPhenos)
{
    rbind(NULL, as.numeric(mpdPhenos[["value"]]))
}

# returns integer values which can be used to order the given chromosome
# strings (ie the returned vector is suitable as an input to the order(...)
# function)
.chrNameValues <- function(chrNames)
{
    origNames <- chrNames
    
    chrNames <- toupper(chrNames)
    chrNames <- sub("(CHR|CHROMOSOME)", "", chrNames)
    chrInt <- suppressWarnings(as.integer(chrNames)) 
    chrInt[chrNames == "X"] <- .Machine$integer.max - as.integer(2)
    chrInt[chrNames == "Y"] <- .Machine$integer.max - as.integer(1)
    chrInt[chrNames == "M"] <- .Machine$integer.max
    
    if(any(is.na(chrInt)))
    {
        stop("unrecognized chromosome name(s): ",
             paste(unique(origNames[is.na(chrInt)]), collapse = ", "))
    }
    
    chrInt
}

.listToDataFrame <- function(theList)
{
    d <- NULL
    for(x in theList)
    {
        if(is.null(d))
        {
            d <- data.frame(x, stringsAsFactors = F)
        }
        else
        {
            d <- cbind(d, data.frame(x, stringsAsFactors = F))
        }
    }
    colnames(d) <- names(theList)
    
    d
}

.chunkRows <- function(m, chunkCount)
{
    rowCount <- nrow(m)
    chunkSize <- rowCount %/% chunkCount
    if((rowCount %% chunkCount) != 0)
    {
        chunkSize <- chunkSize + 1
    }
    
    chunks <- list()
    startRow <- 1
    while(startRow <= rowCount)
    {
        endRow <- min(startRow + chunkSize - 1, rowCount)
        chunks[[length(chunks) + 1]] <- m[startRow : endRow, ]
        startRow <- endRow + 1
    }
    
    chunks
}
