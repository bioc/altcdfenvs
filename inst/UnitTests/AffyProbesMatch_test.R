
library(RUnit)

library(altcdfenvs)

pm_1 <- "AATAATAATAATAATAATAATAAGC"
mm_1 <- "AATAATAATAATTATAATAATAAGC"

pm_2 <- "CCACCACCACCACCACCACCACCTG"
mm_2 <- "CCACCACCACCAGCACCACCACCTG"

pm_3 <- "GGAGGGAGGGAGGGAGGGAGGGACT"
mm_3 <- "GGAGGGAGGGAGCGAGGGAGGGACT"

pm_4 <- "TTATTATTATTATTATTATTATTGC"
mm_4 <- "TTATTATTATTAATATTATTATTGC"


probetable <-
  data.frame(sequence = I(c(pm_1, pm_2, pm_3, pm_4)),
             x = c(10, 11, 13, 13),
             y = c(10, 12, 11, 12),
             Probe.Set.Name = I(c("12_at", "12_at", "m4_at", "4_at")),
             Probe.interrogation.Position = c(100, 120, 130, 140),
             Target.Strandeness = factor(rep("Antisense", 4)))

class(probetable) <- c("data.frame", "probetable")

## test mmProbe
mmp <- mmProbes(probetable)
checkIdentical(c(mm_1, mm_2, mm_3, mm_4), mmp)



target_1 <- paste(pm_1, pm_2, sep="GCGCG")
target_2 <- paste(pm_1, pm_4, sep="GCGCG")
target_3 <- paste("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", mm_3, sep="")

targets <- list(t1 = DNAString(target_1),
                t2 = DNAString(target_2),
                t3 = DNAString(target_3))

## test match
apm <- matchAffyProbes(probetable, targets,
                       chip_type = "foo")


checkIdentical(apm@pm[[1]], as.integer(c(1,2)))
checkIdentical(apm@pm[[2]], as.integer(c(1,4)))
checkIdentical(apm@pm[[3]], integer(0))

checkIdentical(apm@mm[[2]], integer(0))
checkIdentical(apm@mm[[3]], as.integer(c(3)))


## test merge

checkIdentical(combine(matchAffyProbes(probetable, targets[1:2], "foo"),
                       matchAffyProbes(probetable, targets[3], "foo")),
               matchAffyProbes(probetable, targets, "foo"))



## toHypergraph

hg <- toHypergraph(apm)


## test build env

altCdf <- buildCdfEnv.biostrings(apm, nrow.chip = 15, ncol.chip = 15)

##

checkIdentical(hg, toHypergraph(altCdf))


altenv <- as(altCdf, "environment")


#

cdfenv <- new.env(hash = TRUE, parent=emptyenv())
