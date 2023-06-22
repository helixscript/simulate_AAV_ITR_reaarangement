library(dplyr)
library(Biostrings)
library(ShortRead)

# Read in 5K sites synthetic reads; trim simulated over-reading from R1 reads.
system('cutadapt -a TGCTAGAGATTTTC -o trimmed_R1.fastq.gz ~/projects/AAVengeR_synReads/data/combined_R1.fastq.gz')
R1 <- readFastq('trimmed_R1.fastq.gz')
I1 <- readFastq('~/projects/AAVengeR_synReads/data/combined_I1.fastq.gz')
R2 <- readFastq('~/projects/AAVengeR_synReads/data/combined_R2.fastq.gz')

set.seed(1)
minFragLength <- 15
maxFragLength <- 50
numFragments  <- 5000
minITRLength  <- 75
minRecombPos  <- 15
fullITRseq    <- 'TCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCCAGCTG'

vector <- as.character(readDNAStringSet('ITR_scramble_vector.fasta'))

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




i <- Biostrings::vcountPattern('TGCTA', R1@sread)
# Remove TGCTAGAGATTTTC

o <- stringr::str_extract(as.character(I1@id), '\\d+_chr[\\d+XY]+')

i <- 1
ITRseqs$id <- NA

R2.2 <- Reduce('append', lapply(split(R2, o), function(x){
          o <- x@sread
          o <- subseq(o, 15, width(o))
          o <- DNAStringSet(paste0(ITRseqs[i,]$seq, as.character(o)))
  
          x@sread <- o
          x@quality <- ShortRead::FastqQuality(sapply(width(x), function(a) paste0(rep('D', a), collapse = '')))
  
          a <- unlist(strsplit(unique(stringr::str_extract(as.character(x@id), '\\d+_chr[\\d+XY]+')), '_'))
          b <- unique(stringr::str_extract(as.character(x@id), 'subject_\\d+'))
          
          ITRseqs[i,]$id <<- paste0(b, ',', a[2], '+', a[1])
          i <<- i + 1
          x
        }))

writeFastq(I1,   'ITR_scramble_I1.fastq.gz')
writeFastq(R1,   'ITR_scramble_R1.fastq.gz')
writeFastq(R2.2, 'ITR_scramble_R2.fastq.gz')

readr::write_tsv(ITRseqs, 'ITR_scramble_truth.tsv')
