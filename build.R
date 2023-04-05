library(dplyr)
library(Biostrings)

vector <- as.character(readDNAStringSet('vector.fasta'))

set.seed(1)
minFragLength <- 15
maxFragLength <- 50
numFragments  <- 1000
minITRLength  <- 100
minRecombPos  <- 15
fullITRseq    <- 'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCCAGCTG'

ITRseqs <- bind_rows(lapply(c(1:numFragments), function(x){
  n <- 0
  while(n < minRecombPos) n <- sample(1:nchar(fullITRseq), 1)
  
  read <- substr(fullITRseq, 1, n)
  map <- paste0('1..', n, '[', paste0(stringr::str_locate(vector, read), collapse = '+'), '];')
  
  unknownSeg <- FALSE
  
  while(nchar(read) < minITRLength){
    segStart <- sample(1:(nchar(vector) - maxFragLength), 1)
    segEnd <- segStart + sample(minFragLength:maxFragLength, 1)
    seg <- substr(vector, segStart, segEnd)
    
    strand <- ifelse(sample(1:2, 1) == 1, '+', '-')

    loc <-  paste0(stringr::str_locate(vector, seg)[1], strand, stringr::str_locate(vector, seg)[2])
    
    # 1 out of 20 chance that the segment is unknown.
    if(sample(1:20, 1) == 1 & unknownSeg == FALSE){
      unknownSeg <- TRUE
      loc <- 'x'
      seg <- paste0(rep('N', nchar(seg)), collapse = '')
    }
    
    map <- paste0(map, n+1, '..', n+nchar(seg), '[', loc, '];')
    
    n <- n + nchar(seg)
    
    if(strand == '-') seg <- as.character(reverseComplement(DNAStringSet(seg)))
    read <- paste0(read, seg)
  }
  
  tibble(map = sub(';$', '', map), seq = read)
}))

readr::write_tsv(ITRseqs, 'simulatedITRs.tsv')


