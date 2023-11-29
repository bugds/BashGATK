library(ExomeDepth)



setwd(curr.recalibrated.folder)

exons.hg38 <- read.table("exome_table_file_path", sep = "\t", header = T)
genes.hg38 <- read.table("gene_table_file_path", sep = "\t", header = T)
genes.hg38.GRanges <- GenomicRanges::GRanges(
    seqnames = genes.hg38$chrom,
    IRanges::IRanges(
        start = genes.hg38$exonStarts,
        end = genes.hg38$exonEnds
    ),
    names = genes.hg38$name2
)

dgv.hg38 <- read.csv("DGV_merged", header = TRUE, sep = ',')
dgv.hg38.GRanges <- GenomicRanges::GRanges(
    seqnames = dgv.hg38$chromosome,
    IRanges::IRanges(
        start=dgv.hg38$start,
        end=dgv.hg38$end
    ),
    names = dgv.hg38$dgv_name
)

bams <- c(bam_files)

ref.fasta <- "ref_fasta_path"

counts <- getBamCounts(
    bed.frame = exons.hg38,
    bam.files = bams,
    include.chr = F,
    referenceFasta = ref.fasta
)

counts.df <- as(counts, 'data.frame')

test.set <- tests

for (t in test.set) {
    curr.control.set <- test.set[test.set != t]
    control.samples <- as.matrix(counts.df[, curr.control.set])

    choice <- select.reference.set(
        test.counts <- counts[[t]],
        reference.count <- control.samples,
        bin.length <- (counts.df$end - counts.df$start)/1000,
        n.bins.reduced <- 10000
    )

    ref.mtrx <- as.matrix(counts.df[, choice$reference.choice, drop = FALSE])
    reference.selected <- apply(X = ref.mtrx, MAR = 1, FUN = sum)

    all.exons <- new(
        'ExomeDepth',
        test = counts[[t]],
        reference = reference.selected,
        formula = 'cbind(test, reference) ~ 1'
    )

    all.exons <- CallCNVs(
        x = all.exons,
        transition.probability = 10^-4,
        chromosome = counts$chromosome,
        start = counts$start,
        end = counts$end,
        name = counts$exon
    )

    all.exons <- AnnotateExtra(
        x = all.exons,
        reference.annotation = dgv.hg38.GRanges,
        min.overlap = 0.5,
        column.name = 'dgv.hg38'
    )

    all.exons <- AnnotateExtra(
        x = all.exons,
        reference.annotation = genes.hg38.GRanges,
        min.overlap = 0.0001,
        column.name = 'genes.hg38'
    )

    output.file <- paste("exome_calls_path", paste(t, 'exome_calls.csv', sep = '_'), sep = '/')

    write.csv(
        file = output.file,
        x = all.exons@CNV.calls,
        row.names = FALSE
    )

    df <- read.csv(output.file)
    df <- df[order(df$BF, decreasing=TRUE),]
    df <- df[order(df$dgv.hg38, decreasing=TRUE),]
    top10 <- head(df, n = 10)
    write.csv(df, output.file)

    pdf(file= paste("exome_calls_path", paste(t, 'pdf', sep = '.'), sep = '/')) 
    for(i in 1:nrow(top10)) {
        row <- top10[i,]
        plot(
            all.exons,
            sequence = row$chromosome,
            xlim = c(row$start - 100000, row$end + 100000),
            count.threshold = 20,
            main = paste(row$type, 'in', row$chromosome, sep = ' '),
            cex.lab = 0.8,
            with.gene = TRUE
        )
    }
    dev.off()

    output.rds <- paste("exome_calls_path", paste(t, 'rds', sep = '.'), sep = '/')

    saveRDS(all.exons, output.rds)
}