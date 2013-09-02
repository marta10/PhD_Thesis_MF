support <- read.table('support/support_with_extra.tab', header = TRUE)
good.samples <- as.character(subset(support, unknown_FH == 1, 'ID', drop = TRUE))
source('/cluster/project4/vyp/vincent/Software/pipeline/GATK/process_multiVCF.R')

message('Read the all pass table for filtering')
all.pass <- read.csv('../../UK10K/all_cohorts/data/all_pass.exome_summary.csv', na.string = c('NA', ''), stringsAsFactors = FALSE, col.names =  c('Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','ESP6500_ALL','X1000g2012apr_ALL','dbSNP137','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred','LJB_GERP++', 'cg69', 'Omim', 'Chr','Start','End','Ref','Obs','Call', 'QUAL', 'Depth', 'junk1', 'junk2'))

all.pass$signature <- paste(all.pass$Chr, all.pass$Start, all.pass$Ref, all.pass$Obs, sep = '_')


message('Read the case control data')
cases <- read.table('data/flagged_UK10K_FH_cases.csv', sep = '\t', header = TRUE, na.string = c('NA', ''), stringsAsFactors = FALSE)
controls <- read.csv('../UK10K_controls_for_FH/annovar_UK10K_controls_for_FH.tab.exome_summary.csv', na.string = c('NA', ''), stringsAsFactors = FALSE)

controls <- annotate.standard.annovar.output ( controls )

names(controls)[30:31] <- c('AN', 'AC')

cases <- subset(cases, ID %in% good.samples)
n.cases <- length(unique(cases$ID))

SNP.only <- TRUE
if (SNP.only) {
message('Number of calls in controls/cases ', nrow(controls), ' ', nrow(cases))
controls <- subset(controls, nchar(Obs) == 1)
cases <- subset(cases, nchar(Obs) == 1)
message('Number of calls in controls/cases ', nrow(controls), ' ', nrow(cases))
}

cases$Gene <- gsub(cases$Gene, pattern = '\\(.*', replacement = '')
controls$Gene <- gsub(controls$Gene, pattern = '\\(.*', replacement = '')

cases$signature <- paste(cases$Chr, cases$Start, cases$Ref, cases$Obs, sep = '_')
controls$signature <- paste(controls$Chr, controls$Start, controls$Ref, controls$Obs, sep = '_')


my.sigs <- sort(table(subset(cases$signature, cases$rare)), decreasing = TRUE)
bad.sigs <- names(subset(my.sigs, my.sigs > 10))
extra.exclude <- c(bad.sigs, '16_67229794_CAG_-', 'X_54566663_AAA_-', 'X_142121809_AGA_-', '1_52499071_GCCCGCCACCTGAGGTCCCGCGATCGG_-', '17_77807917_-_GCCGCC', '17_77807926_-_GCCGCC', '17_77807917_T_TGCCGCC', '10_33136819_AA_-', '10_46999591_-_ATGAGGGAG')

exclude.common <- subset(controls, AC > 50, 'signature', drop = TRUE) ##variants that should really be annotated as common                    
cases.filtered <- subset(cases, signature %in% all.pass$signature & ! signature %in% exclude.common & ! signature %in% extra.exclude)
controls.filtered <- subset(controls, ! signature %in% exclude.common & ! signature %in% extra.exclude)

cases.filtered <- subset(cases.filtered, (splicing | non.syn | lof) & ((SegDup < 0.91) | is.na(SegDup)) )
controls.filtered <- subset(controls.filtered,  (splicing | non.syn | lof) & (( SegDup < 0.91) | is.na(SegDup)))

### LDLR check                                                                                                                                
my.prob <- c(119, 1300)/(1300+119)
n.LDLR.controls <- sum(subset(controls, HUGO.no.splice.info == 'LDLR' & somewhat.rare, 'AC', drop = TRUE))
n.LDLR.cases <- nrow(subset(cases.filtered, HUGO.no.splice.info == 'LDLR'))
LDLR.pvalue <- chisq.test(c(22, 54), p = c(119, 1300)/(1300+119))

all.genes <- unique(c(cases.filtered$HUGO.no.splice.info, controls$HUGO.no.splice.info))

print(sum(subset(controls, Gene == 'ABCA1' & rare)$AC))


####                                                                                                                                          
controls.rare <- subset(controls.filtered, rare)
cases.rare <- subset(cases.filtered, rare)

my.count.controls <- tapply(controls.rare$AC, IND = controls.rare$HUGO.no.splice.info, FUN = sum)
my.count.cases <- table(factor(cases.rare$HUGO.no.splice.info, levels = all.genes))



my.frame <- data.frame(gene = names(my.count.cases),
count.cases = as.numeric(my.count.cases))


my.frame.temp <- data.frame(gene = names(my.count.controls),
count.controls = as.numeric(my.count.controls))

combined <- merge(my.frame, my.frame.temp, by = 'gene', all = TRUE)
combined$count.cases <- ifelse (is.na(combined$count.cases), 0, combined$count.cases)
combined$count.controls <- ifelse (is.na(combined$count.controls), 0, combined$count.controls)

my.mat <- as.matrix(combined[, 2:3])
combined$my.p <- ifelse ( my.mat[,1] > 0, pbinom(q = my.mat[,1]-1, size = my.mat[,1] + my.mat[,2], prob = n.cases/(n.cases+1300), lower.tail = FALSE), 1)
combined <- combined[ order(combined$my.p), ]
##############                                                                                                                                

print(subset(combined, my.p < 10^(-4)))
print(subset(combined, gene == 'CH25H'))
write.csv(x = combined, file = 'data/pvalues_functional_multisamples_messy.csv', row.names = FALSE)

