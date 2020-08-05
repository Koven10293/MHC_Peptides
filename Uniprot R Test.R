setwd("E:/Users/joshu/Assignments/BCH3990/R Working Directories")

install.packages("UniprotR")

library("UniprotR")

test_seq <- GetSequences("Q53EL6")

test_seq$Sequence
