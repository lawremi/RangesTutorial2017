library(RangesTutorial2017)
calls <- readVcf(system.file("extdata", "lumpy.vcf.gz",
                             package="RangesTutorial2017"))
truth <- readVcf(system.file("extdata", "truth.vcf.bgz",
                             package="RangesTutorial2017"))

truth <- subset(truth, SVTYPE=="DEL")
calls <- subset(calls, SVTYPE=="DEL")

seqlevelsStyle(calls) <- "NCBI"
truth <- keepStandardChromosomes(truth,
                                 pruning.mode="coarse")

calls <- as(calls, "VRanges")
truth <- as(truth, "VRanges")

ref(truth) <- "."

ref(truth) <- "."

calls <- calls[called(calls)]

hits <- findOverlaps(truth, calls)

hits <- findOverlaps(truth, calls)

hits <- as(hits, "List")
call_rl <- extractList(ranges(calls), hits)
dev <- abs(start(truth) - start(call_rl)) +
    abs(end(truth) - end(call_rl))

dev_ord <- order(dev)
keep <- phead(dev_ord, 1L)
truth$deviance <- drop(dev[keep])
truth$call <- drop(hits[keep])

library(ggplot2)
rdf <- as.data.frame(truth)
ggplot(aes(x=deviance),
       data=subset(rdf, deviance <= 500)) +
    stat_ecdf() + ylab("fraction <= deviance")

truth$called <- with(truth,
                     !is.na(deviance) & deviance <= 300)

mean(truth$called)

mean(truth$called)

calls$fp <- TRUE
calls$fp[subset(truth, called)$call] <- FALSE

mean(calls$fp)

bed <-
    system.file("extdata", "altRegions.GRCh38.bed.gz",
                package="RangesTutorial2017")
altRegions <- import(bed)
seqlevelsStyle(altRegions) <- "NCBI"
altRegions <-
    keepStandardChromosomes(altRegions,
                            pruning.mode="coarse")

calls$inAlt <- calls %over% altRegions
xtabs(~ inAlt + fp, calls)
