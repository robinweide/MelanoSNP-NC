file$ncrank <- rank(-file$V3)
file$crank <- rank(-file$V4)
file$caddrank <- rank(-file$V5)
file$dannrank <- rank(-file$V6)
file$meanRanks=rowMeans(file[,c("ncrank","crank","caddrank","dannrank")], na.rm=TRUE)
file$rank <- rank(file$meanRanks)
file_sorted <- file[order(file[,14] ),]
colnames(file_sorted) <- c("coord", "mut", "F2nc", "F2c", "CADD", "DANN", "Case", "Control", "Rank-F2nc", "Rank-F2c", "Rank-CADD", "Rank-DANN", "Rank-mean", "Rank")
output <- file_sorted[,c(1,2,7,8,14)]
write.table(output, "rareNcCandidates.tsv", sep="\t")

