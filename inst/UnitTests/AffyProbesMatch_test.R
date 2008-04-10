library(RUnit)

library(altcdfenvs)

pm_1 <- "AATAATAATAATAATAATAATAAGC"
mm_1 <- "AATAATAATAATTATAATAATAAGC"

pm_2 <- "CCACCACCACCACCACCACCACCTG"
mm_2 <- "CCACCACCACCAGCACCACCACCTG"

pm_3 <- "TTATTATTATTATTATTATTATTGC"
mm_3 <- "TTATTATTATTAATATTATTATTGC"

pm_4 <- "GGAGGGAGGGAGGGAGGGAGGGACT"
mm_4 <- "GGAGGGAGGGAGCGAGGGAGGGACT"

probetable <-
  data.frame(sequence = I(c(pm_1, pm_2, pm_3, pm_4)),
             x = c(100, 200, 300, 400),
             y = c(100, 200, 300, 400),
             Probe.Set.Name = I(c("12_at", "12_at", "3_at", "m3_at")),
             Probe.interrogation.Position = c(100, 200, 300, 400),
             Target.Strandeness = factor(rep("Antisense", 4)))

class(probetable) <- c("data.frame", "probetable")

## test mmProbe
mmp <- mmProbes(probetable)
checkIdentical(c(mm_1, mm_2, mm_3, mm_4), mmp)



target_1 <- paste(pm_1, pm_2, sep="GCGCG")
target_2 <- paste(pm_1, pm_3, sep="GCGCG")
target_3 <- paste("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", mm_4, sep="")

targets <- list(t1 = DNAString(target_1),
                t2 = DNAString(target_2),
                t3 = DNAString(target_3))

## test match
apm <- matchAffyProbes(probetable, targets,
                       chip_type = "foo")


checkIdentical(apm@pm[[1]], as.integer(c(1,2)))
checkIdentical(apm@pm[[2]], as.integer(c(1,3)))
checkIdentical(apm@pm[[3]], integer(0))

checkIdentical(apm@mm[[2]], integer(0))
checkIdentical(apm@mm[[3]], as.integer(c(4)))


## test merge

checkIdentical(combine(matchAffyProbes(probetable, targets[1:2], "foo"),
                       matchAffyProbes(probetable, targets[3], "foo")),
               matchAffyProbes(probetable, targets, "foo"))



##

