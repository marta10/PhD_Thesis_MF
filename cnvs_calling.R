library(ExomeDepth)

data(exons.hg19)




my.bam <- list.files('/Volumes/Expansion_Drive/BAMs/batch_2/', pattern = ".bam$", full.names = TRUE)
my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = my.bam)


###
my.reference.bam <- list.files('/Volumes/Seagate Expansion Drive/BAMs_Jan_2012/HYP/batch_5/', pattern = ".bam$", full.names = TRUE)

my.reference.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = my.reference.bam)

my.reference.counts.dafr <- as(my.reference.counts[, colnames(my.reference.counts)], 'data.frame')

my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

print(head(my.counts.dafr))



test <- new('ExomeDepth',
            test = my.counts.dafr$UK10K_HYP5315266.bam,
            reference = my.counts.dafr$UK10K_HYP5315268.bam,
            formula = 'cbind(test, reference) ~ 1',
            subset.for.speed = seq(1, nrow(my.counts.dafr), 100))

show(test)


##building a reference set

my.test <- my.counts.dafr$UK10K_HYP5269570.bam

my.ref.samples <- c('UK10K_HYP5269571.bam', 'UK10K_HYP5269572.bam', 'UK10K_HYP5269573.bam', 'UK10K_HYP5269574.bam', 'UK10K_HYP5269575.bam', 'UK10K_HYP5269576.bam', 'UK10K_HYP5269577.bam', 'UK10K_HYP5269578.bam', 'UK10K_HYP5269581.bam', 'UK10K_HYP5269585.bam', 'UK10K_HYP5269589.bam', 'UK10K_HYP5269595.bam', 'UK10K_HYP5269597.bam', 'UK10K_HYP5269598.bam', 'UK10K_HYP5269601.bam', 'UK10K_HYP5269602.bam', 'UK10K_HYP5269604.bam', 'UK10K_HYP5269605.bam', 'UK10K_HYP5269606.bam', 'UK10K_HYP5269607.bam', 'UK10K_HYP5269608.bam')
my.ref.samples <- c('UK10K_HYP5315266.bam', 'UK10K_HYP5315268.bam', 'UK10K_HYP5315271.bam', 'UK10K_HYP5315273.bam')
my.reference.set <- as.matrix(my.counts.dafr[, my.ref.samples])
my.choice <- select.reference.set (test.counts = my.test,
                                   reference.counts = my.reference.set,
                                   bin.length = (my.counts.dafr$end - my.counts.dafr$start)/1000,
                                   n.bins.reduced = 10000)

print(my.choice[[1]])

my.matrix <- as.matrix (my.counts.dafr[, my.choice$reference.choice, drop = FALSE])
my.reference.selected <- apply (X = my.matrix,
                                MAR = 1,
                                FUN = sum)

## CNVs calling

all.exons <- new('ExomeDepth',
                 test = my.test,
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = my.counts.dafr$space,
                      start = my.counts.dafr$start,
                      end = my.counts.dafr$end,
                      name = my.counts.dafr$names)

print(head(all.exons@CNV.calls))

output.file <- 'cnv_calls_HYP5315275.csv'
write.csv(file = output.file,
          x = all.exons@CNV.calls,
          row.names = FALSE)

## annotating CNVs using Conrad et al Nature 2010 data for common CNVs
data(Conrad.hg19)
head(Conrad.hg19.common.CNVs)

all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = Conrad.hg19.common.CNVs,
                           min.overlap = 0.0001,
                           column.name = 'Conrad.hg19')

print(head(all.exons@CNV.calls, n=30))


## plotting

plot (all.exons,
      sequence = '19',
      xlim = c(11200038 - 1000, 11244505 + 1000),
      count.threshold = 20,
      main = 'LDLR gene',
      with.gene = TRUE)

pdf(image = "/Users/martafutema/Documents/UK10K_FH_analysis/")

#getting error while plotting: Error in 0:size : NA/NaN argument



## looping over muliple samples

data(Conrad.hg19)

exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome,
                              IRanges(start=exons.hg19$start, end=exons.hg19$end),
                              names = exons.hg19$name)

my.matrix <- as.matrix(my.counts.dafr[, grep(names(my.counts.dafr), pattern = 'UK10K_HYP.*')])

nsamples <- ncol(my.matrix)

#### start looping over each sample

for (i in 1:nsamples) {
  
  my.choice <- select.reference.set (test.counts = my.matrix [,i],
                                     reference.counts = my.matrix [,-i],
                                     bin.length = (my.counts.dafr$end - my.counts.dafr$start)/1000,
                                     n.bins.reduced = 10000)
  
  my.reference.selected <- apply(X = my.matrix[, my.choice$reference.choice, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
  
  message('Now creating the ExomeDepth object')
  all.exons <- new('ExomeDepth',
                   test = my.matrix[,i],
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = my.counts.dafr$space,
                        start = my.counts.dafr$start,
                        end = my.counts.dafr$end,
                        name = my.counts.dafr$names)
  
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.hg19.GRanges,
                             min.overlap = 0.0001,
                             column.name = 'exons.hg19')
  
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = Conrad.hg19.common.CNVs,
                             min.overlap = 0.5,
                             column.name = 'Conrad.hg19')
  
  names <- names(my.counts.dafr[6:ncol(my.counts.dafr)])
  names <- gsub(names, pattern = ".bam", replacement="")
  output.file <- paste(names[i], 'CNVcalls.csv', sep='')
  dat <- cbind(names[i], all.exons@CNV.calls)
  
  write.table(dat, "CNV_calls_batch2.csv", row.names = FALSE, quote=F, sep= ",", append=T, col.names=F)
 # write.table(file = output.file, x = all.exons@CNV.calls, row.names = FALSE, quote=F, sep= ",", append=T, col.names=F)
  
 pdf(paste(names[i], 'LDLR_CNVs.pdf', sep=''))
  plot (all.exons,
        sequence = '19',
        xlim = c(11200038 - 1000, 11244505 + 1000),
        count.threshold = 20,
        main = 'LDLR gene',
        with.gene = TRUE)
  dev.off()
  
}




