---
title       : STAT435
subtitle    : Assignment 3
author      : Jack Hervey
job         : STAT435
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : []            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
---

## Libraries and Data
```
source("http://bioconductor.org/biocLite.R")
biocLite("pamr")
biocLite("gplots")
biocLite("multtest")
biocLite("survival")
```

```r
library(survival)
library(pamr)
library(gplots)
library(multtest)
```

```
load("./STAT435-VDVdata.RData")
```

---

## Rescale data to zero mean and standard deviation 1

```r
yy <- t(apply(vdv.int, 1, function(x) (x - mean(x))/sd(x)))
yy[yy < -3] <- -3
yy[yy > 3] <- 3
```


---

## Heatmap of clustered data

```r
heatmap(yy, distfun = function(x) as.dist(1 - cor(t(x))), col = c(rgb(0, seq(1, 
    0, l = 25), 0), rgb(seq(0, 1, l = 25), 0, 0)), labRow = "", labCol = "")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 


---

## Heatmap of clustered data
- There appears to be approximately 3 major gene clusters (side dendrogram).

- There are four distinct groups of patients (top dendrogram).

---

## Split data into 4 groups
Number of women assigned to each group:

```r
ctint <- cutree(hclust(as.dist(1 - cor(yy))), 4)
table(ctint)
```

```
## ctint
##   1   2   3   4 
##  40  79 114  62
```

This is similar to the number of women assigned to the different intrinsic subtypes:

```r
table(vdv.clin$Intrinsic.Subtype)
```

```
## 
##       Basal        HER2        LumA        LumB Normal-like 
##          53          35         123          69          15
```


---

## Instrinsic subtypes of the four groups

```r
table(vdv.clin$Intrinsic.Subtype, ctint, dnn = c("Subtype", "Group"))
```

```
##              Group
## Subtype        1  2  3  4
##   Basal        1  0 50  2
##   HER2         9  0 24  2
##   LumA         3 65 19 36
##   LumB        27  7 13 22
##   Normal-like  0  7  8  0
```


---

## Relapse Rates

```r
table(vdv.clin$RFS.event, vdv.clin$Intrinsic.Subtype, dnn = c("Relapse rate", 
    "Subtype"))
```

```
##             Subtype
## Relapse rate Basal HER2 LumA LumB Normal-like
##            0    26   17   99   28           7
##            1    27   18   24   41           8
```

```r
table(vdv.clin$RFS.event, ctint, dnn = c("Relapse rate", "Group"))
```

```
##             Group
## Relapse rate  1  2  3  4
##            0 17 67 58 35
##            1 23 12 56 27
```


---

## Death Rates

```r
table(vdv.clin$OS.event, vdv.clin$Intrinsic.Subtype, dnn = c("Death rate", "Subtype"))
```

```
##           Subtype
## Death rate Basal HER2 LumA LumB Normal-like
##          0    28   20  109   46          13
##          1    25   15   14   23           2
```

```r
table(vdv.clin$OS.event, ctint, dnn = c("Death rate", "Group"))
```

```
##           Group
## Death rate  1  2  3  4
##          0 24 73 72 47
##          1 16  6 42 15
```


---

## Singular value decomposition metagene

```r
xx <- svd(yy)
attach(xx)
mg <- svd(yy)$v[, 1]
mg.rnk <- rank(mg)/length(mg)
mg.col <- greenred(length(mg))[rank(mg)]
ord <- order(mg.rnk)
```


---

## Metagene and heatmap

```r
heatmap.2(yy[, ord], trace = "none", scale = "none", Colv = T, Rowv = T, ColSideColors = mg.col[ord], 
    col = greenred(50), labCol = paste("Sample", 1:ncol(yy))[ord], labRow = paste("Gene", 
        1:ncol(yy)))
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


---

## Boxplot showing distribution of metagene values accross intrinsic subtype

```r
boxplot(mg.rnk ~ vdv.clin$Intrinsic.Subtype, xlab = "Intrinsic.Subtype", ylab = "Metagene Value")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 

---

## Boxplot showing distribution of metagene values accross ER status

```r
boxplot(mg.rnk ~ vdv.clin$ER, xlab = "ER", ylab = "Metagene Value")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 

---

## Boxplot showing distribution of metagene values accross grade

```r
boxplot(mg.rnk ~ vdv.clin$Grade, xlab = "Grade", ylab = "Metagene Value")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


---

## Information
- Subtype and estrogen receptor status appear to gave the greatest impact on proliferation, followed by grade. Tumour grade shows an obvious upward trend, as expected, but the ranges are very broad.

---

## Classification and survival analysis
Split van de Vijver data into a training set and test set.

```r
tr <- sample(seq(1, 295), 200)
ts <- setdiff(1:295, tr)
tr.dat <- list()
tr.dat$x <- vdv.genes[, tr]
tr.dat$y <- vdv.clin$RFS.event[tr]
ts.dat <- list()
ts.dat$x <- vdv.genes[, ts]
ts.dat$y <- vdv.clin$RFS.event[ts]
```


---

## Relapse and non-relapse class between training and test sets

```r
table(tr.dat$y)
```

```
## 
##   0   1 
## 122  78
```

```r
table(ts.dat$y)
```

```
## 
##  0  1 
## 55 40
```


---

## PAMR and CV plot

```r
train <- pamr.train(tr.dat)
cv <- pamr.cv(train, tr.dat)
```


```r
pamr.plotcv(cv)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 


---

## Geneplots

```r
pamr.geneplot(train, tr.dat, threshold = 2.6)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 


---

## Predictions

```r
pred.ts <- pamr.predict(train, ts.dat$x, threshold = 1.5)
table(ts.dat$y, pred.ts)
```

```
##    pred.ts
##      0  1
##   0 36 19
##   1 13 27
```

```r
pred.tr <- pamr.predict(train, tr.dat$x, threshold = 1.5)
table(tr.dat$y, pred.tr)
```

```
##    pred.tr
##      0  1
##   0 87 35
##   1 26 52
```


---
  
## Fisher test on test set

```r
fisher.test(table(ts.dat$y, pred.ts))
```

```
## 
## 	Fisher's Exact Test for Count Data
## 
## data:  table(ts.dat$y, pred.ts)
## p-value = 0.001885
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##   1.53 10.27
## sample estimates:
## odds ratio 
##      3.874
```


---
  
## Fisher test on training set

```r
fisher.test(table(tr.dat$y, pred.tr))
```

```
## 
## 	Fisher's Exact Test for Count Data
## 
## data:  table(tr.dat$y, pred.tr)
## p-value = 2.095e-07
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  2.580 9.622
## sample estimates:
## odds ratio 
##      4.927
```


---

## Sensitivity and specificty
Test set


- Sensitivity: 67.5 %
- Specificity: 65.4545 %
- PPV: 58.6957 %
- NPV: 73.4694 %

Training set


- Sensitivity: 66.6667 %
- Specificity: 71.3115 %
- PPV: 59.7701 %
- NPV: 76.9912 %

As expected the training set performs much better, as it was used to create the predictor.

---

## Survival analysis

```r
tr.clin <- vdv.clin[tr, ]
ts.clin <- vdv.clin[ts, ]
```


---

## Survival plot of training set

```r
plot(survfit(Surv(tr.clin$RFS.months, tr.clin$RFS.event) ~ pred.tr))
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25.png) 


---

## Survival plot of training set

```r
survdiff(Surv(tr.clin$RFS.months, tr.clin$RFS.event) ~ pred.tr)
```

```
## Call:
## survdiff(formula = Surv(tr.clin$RFS.months, tr.clin$RFS.event) ~ 
##     pred.tr)
## 
##             N Observed Expected (O-E)^2/E (O-E)^2/V
## pred.tr=0 113       26       51      12.3      35.9
## pred.tr=1  87       52       27      23.2      35.9
## 
##  Chisq= 35.9  on 1 degrees of freedom, p= 2.09e-09
```


---

## Survival plot of test set

```r
plot(survfit(Surv(ts.clin$RFS.months, ts.clin$RFS.event) ~ pred.ts))
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27.png) 


---

## Survival plot of test set

```r
survdiff(Surv(ts.clin$RFS.months, ts.clin$RFS.event) ~ pred.ts)
```

```
## Call:
## survdiff(formula = Surv(ts.clin$RFS.months, ts.clin$RFS.event) ~ 
##     pred.ts)
## 
##            N Observed Expected (O-E)^2/E (O-E)^2/V
## pred.ts=0 49       13     23.1      4.45      10.6
## pred.ts=1 46       27     16.9      6.11      10.6
## 
##  Chisq= 10.6  on 1 degrees of freedom, p= 0.0011
```


---

## Survival plot using Luminal A subtype predictor

```r
int.pred <- ifelse(vdv.clin$Intrinsic.Subtype == "LumA", 0, 1)
table(int.pred, vdv.clin$RFS.event)
```

```
##         
## int.pred  0  1
##        0 99 24
##        1 78 94
```


---

## Survival plot using Luminal A subtype predictor

```r
fisher.test(table(int.pred, vdv.clin$RFS.event))
```

```
## 
## 	Fisher's Exact Test for Count Data
## 
## data:  table(int.pred, vdv.clin$RFS.event)
## p-value = 1.012e-09
## alternative hypothesis: true odds ratio is not equal to 1
## 95 percent confidence interval:
##  2.823 8.894
## sample estimates:
## odds ratio 
##      4.943
```


---

## Survival plot using Luminal A subtype predictor

```r
plot(survfit(Surv(vdv.clin$RFS.months, vdv.clin$RFS.event) ~ int.pred))
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31.png) 


---

## Survival plot using Luminal A subtype predictor

```r
survdiff(Surv(vdv.clin$RFS.months, vdv.clin$RFS.event) ~ int.pred)
```

```
## Call:
## survdiff(formula = Surv(vdv.clin$RFS.months, vdv.clin$RFS.event) ~ 
##     int.pred)
## 
##              N Observed Expected (O-E)^2/E (O-E)^2/V
## int.pred=0 123       24     56.3      18.5      35.5
## int.pred=1 172       94     61.7      16.9      35.5
## 
##  Chisq= 35.5  on 1 degrees of freedom, p= 2.5e-09
```


---

## Sensitivity and specificity of Luminal A subtype predictor
## Sensitivity and specificty


- Sensitivity: 54.6512 %
- Specificity: 80.4878 %
- PPV: 79.661 %
- NPV: 55.9322 %

The Luminal subtype predictor perorms much better than the predictor produced by PAMR. This is to be expected since the subtype has a real biological meaning and therefore related far better to the data collected than the PAMR predictor, which is based off a sample of the total data collected.

---
