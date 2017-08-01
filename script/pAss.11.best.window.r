#!/usr/bin/env Rscript

#' # Accessing the distribution of sequences

#' ## Introduction
#' Based on the analysis of the sliding windows we isolate the smallest window with maximum diversity ie. most number of contigs. And output contigs in that region with >100 bp length. 
#'
#' **NOTE**: This is a generic summary, the same Rscript was applied all KOs. we chose K00001, alcohol dehydrogenase

#+ libraryload, message=FALSE, echo=FALSE


#' Install if missing
if (!require("Biostrings")){source("http://bioconductor.org/biocLite.R"); biocLite("Biostrings")}
if (!require("ggplot2"))    install.packages("ggplot2", repos="http://cran.bic.nus.edu.sg/")
if (!require("dplyr"))      install.packages("dplyr", repos="http://cran.bic.nus.edu.sg/")

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_bw())

#' ## Data processing
#+ Reading-in-KO-path, echo=FALSE
params = commandArgs(T)
#Depending on the KO that was supplied
if(length(params)==0){
    windowFile = "/export2/home/uesu2/reDiamond/out/pAss.10/K00001/K00001"
    msaFile    = "/export2/home/uesu/reDiamond/out/assm.0500/K00001.msa"
    outputDIR  = "."
    contigs = "/export2/home/uesu/simulation_fr_the_beginning/reAssemble/scg/data/newbler/K00927/454AllContigs.fna"
    params = c(windowFile, msaFile, outputDIR)
}
df        = read.table(params[1], h = T)
df[is.na(df)] = 0   #remove na
sequences = readDNAStringSet(filepath = params[2])
outputDIR = params[3]
fullSeq = params[4]
dir.create(outputDIR)

#+ koBackground
ko = tail(strsplit(params[1], '/')[[1]], n=1)
#koInfo = koname(ko, minimal=FALSE)

#' For this report, the KO chosen is `r ko`, `r unique(koInfo$ko.definition)` which belongs to the above following pathways:

#+
#koInfo$ko..pathway.name.


#' ## Choosing MSA loc to sample Diversity

#' We begin by finding a window (in thise case the input was scanned based on window size of 200bp) with the most number of sequences spanning the window.
#'
#' Because there neighbouring regions (+-5bp) with similar characteristics we choose the middle location.

#+ find-location
selectdf = df %>%
    mutate(width=end-start) %>%
    filter(seqInSameWindow == max(seqInSameWindow)) %>%
    filter(width == min(width))
#when there are multiple windows, (usually clustered together) with the same length
if(nrow(selectdf)>1)
    selectdf = selectdf[ceiling(nrow(selectdf)/2),]
#str(selectdf)

#' ## Output contigs of interest 
#' `findLoc` function generates the interval information for extraction the section on the contig which is of interest.

#+ findLoc-function
findLoc = function(sequence, seqname, start, end)
{
    seq            = gsub("-","",sequence)
    capturedSeq    = gsub("-","",substr(sequence, start,end))
    capturedLength = nchar(capturedSeq)
    leftLength  = nchar(gsub("-","", substr(sequence, 1, start)))
    rightLength = nchar(gsub("-","", substr(sequence,end,nchar(sequence))))
    matrix(c(seqname, capturedLength, leftLength, rightLength), nrow=1)
}

#+ ExampleOnly
setNames(as.data.frame(findLoc(sequences[1], as.character(sequences@ranges@NAMES)[1], selectdf$start, selectdf$end),stringsAsFactors=F), c("contig", "capturedLength", "leftLength","rightLength"))
#' We see in this example, the which we've created


#+ output
headerInfo = setNames(as.data.frame(t(
##################################################
    mapply(function(seq,seqname){
           findLoc(seq, seqname, selectdf$start, selectdf$end)
    }, sequences, as.character(sequences@ranges@NAMES))
##################################################
    ), stringsAsFactors=F), c("contig", "capturedLength", "leftLength","rightLength"))

#' Selects only contigs in the window of interest with > 100 basepairs.
#' Note, Spanning here referes to the # of sequences ON THE MSA which are spanning the window
#' the selectedSeqs below are sequences within this window, not neccessarily spanning it which are > 100 bp long

#+
selectedSeq = filter(headerInfo,as.integer(capturedLength)> 100)
#' Selected only the capturedLength > 100

#+ ## Output the selected sequences

#+
sequences2Output <-
lapply(selectedSeq$contig, function(contigName)
       {
           contigInfo = subset(selectedSeq, contig == contigName)
           seqWOGaps = gsub("-","",sequences[which(sequences@ranges@NAMES == contigName)])

           rightLength = ifelse(contigInfo$rightLength == 0, nchar(seqWOGaps), as.integer(nchar(seqWOGaps)) - as.integer(contigInfo$rightLength))
           region = substr(seqWOGaps, as.integer(contigInfo$leftLength), rightLength)
           list(
                sequence    = region, 
                header      = with(selectdf,sprintf("%s ## spanning:%s msaStart:%s msaEND:%s max10BPwindow:%s", contigName, seqInSameWindow, start, end, maxSeq.10bp))
                )
})

# adds the contig's position into the header as well, now its just the location
contigsSequences = readDNAStringSet(file=fullSeq)
contigIndex = names(contigsSequences) %>% gsub("(contig\\d{5}).+", "\\1", .)

locationsDF = lapply(sequences2Output, function(contig){
    #browser()
    thisContig = gsub("(contig\\d{5}).+", "\\1", contig$header)
    fullSequence = contigsSequences[which(contigIndex %in% thisContig)] 
    isStraight = grepl("ntRev", contig$header)
    if(isStraight){
        mi = regexpr(contig$sequence, fullSequence %>% as.character)
        if(mi == -1){
            warning(sprintf("Can't find MDR substring in %s", thisContig))
            mdr = data.frame(
                contig = thisContig,
                start = NA,
                length = NA,
                end = NA
            )
        }else{
            mdr = data.frame(
                contig = thisContig,
                start = mi[1],
                length = attr(mi, "match.length")
            ) %>% mutate(end = start + length)
        }
    }else{
        mi= regexpr(contig$sequence, fullSequence %>% reverse %>% complement %>% as.character)
        if(mi == -1){
            warning(sprintf("Can't find MDR substring in %s", thisContig))
            mdr = data.frame(
                contig = thisContig,
                start = NA,
                length = NA,
                end = NA
            )
        }else{
            mdr = data.frame(
                contig = thisContig,
                start = mi[1],
                length = attr(mi, "match.length")
            ) %>% mutate(end = start + length)
        }
    }
}) %>% do.call(rbind,.)

#Generate DNAStringSet Obj
output = DNAStringSet(lapply(sequences2Output, function(x) DNAString(x$sequence)))
#Annotated with headers
output@ranges@NAMES = lapply(sequences2Output, function(x){ 
    contigName = x$header %>% gsub("(contig\\d{5}).+", "\\1", .)
    locDF = filter(locationsDF, contig == contigName)
    sprintf("%s contigStart: %s contigEnd: %s", x$header, locDF$start, locDF$end)
}) %>% do.call(c,.)
#output to file
writeXStringSet(output, sprintf("%s/%s.fna", outputDIR, ko))

#+ plot, message=FALSE
p0 = qplot(output@ranges@width, geom="histogram")+
ggtitle(sprintf("Selected Contig lengths Total:%s contigs",length(sequences2Output)))
ggsave(p0, file=sprintf("%s/%s.pdf", outputDIR,ko))
#p0

#+ sessionInfo
sessionInfo()
