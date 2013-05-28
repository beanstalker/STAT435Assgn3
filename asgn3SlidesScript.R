setwd("~/Documents/University/Otago/Fourth Year/STAT435/Assignment3/jack")
load("~/Documents/University/Otago/Fourth Year/STAT435/Assignment3/STAT435-VDVdata.RData")
load("~/Documents/University/Otago/Fourth Year/STAT435/Assignment3/rdata.RData")


install.packages("devtools")
require(devtools)
install_github("slidify", "ramnathv")
install_github("slidifyLibraries", "ramnathv")

library(knitr)
library(markdown)
library(slidify)

author()

source("http://bioconductor.org/biocLite.R")
biocLite("pamr")
biocLite("gplots")
biocLite("multtest")
biocLite("survival")

slidify("index.Rmd")

save.image("~/Documents/University/Otago/Fourth Year/STAT435/Assignment3/rdata.RData")