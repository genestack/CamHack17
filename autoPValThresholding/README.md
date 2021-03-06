Automatic P-value thresholding
================
27/09/2017

Imagine a scenario of finding differentially expressed genes between two groups 100 individuals. For each gene, we do a t-test to see if the expression distribution in one group differs from the second:
- Some tests correspond to no differences (noises)
- Some tests correspond to true differences (signals)

The following function, `TestFDRThresholding`, simulates this scenario, where the number of signals and noises are specified by the parameters `number.of.signals` and `number.of.noises`. It then adjusts the P-values into False Discovery Rate (FDR) using the Benjamini-Hochberg method. The parameter `fdr.threshold` will be used to define the significant results and the function will report how many of those are true positives and how many are false positives.

The manual way
--------------

``` r
library(ggplot2)
set.seed(1)
TestFDRThresholding <- function(number.of.signals, number.of.noises, fdr.threshold) {
    # simulate a single t-test between two groups of case subjects and control subjects
    # returns test p-value
    TestCaseControl <- function(case.mean, number.of.samples) {
        cases = rnorm(n=number.of.samples, mean=case.mean, sd=1)
        controls = rnorm(n=number.of.samples, mean=0, sd=1)

        return(t.test(cases, controls, alternative="two.sided")$p.value)
    }

    # simulate t-test p-values for signals
    p.true = replicate(number.of.signals, TestCaseControl(1, 100))
    # simulate t-test p-values for noises
    p.noise = replicate(number.of.noises, TestCaseControl(0, 100))
    p.all = sort(c(p.true, p.noise))

    p.all.adjusted = p.adjust(p.all, method="BH")
    p.threshold = p.all[which.min(p.all.adjusted < fdr.threshold)]
    p.threshold # p-value equivalent of the maximum FDR threshold

    number.of.significant.results = sum(p.all < p.threshold)
    number.of.true.positives = sum(p.true < p.threshold)

    print(paste("# signals:", number.of.signals))
    print(paste("# noises:", number.of.noises))
    print(paste("# FDR threshold:", fdr.threshold))
    print(paste("# significant results:", number.of.significant.results))
    print(paste("# true positives:", number.of.true.positives, "out of", number.of.signals))
    print(paste("# false positives:", number.of.significant.results-number.of.true.positives))
    
    # show histogram
    dat = data.frame(p = c(p.true, p.noise),
        source = rep(c('signal','noise'), c(length(p.true), length(p.noise))))
    ggplot(dat, aes(x=p, fill=source)) +
        geom_histogram(alpha = 0.2, bins=100) +
        scale_fill_manual(name="source",values=c("red","blue"),labels=c("noise","signal")) +
        geom_vline(xintercept=p.threshold)
}
```

If we have 100 signals and 10000 noises and set `fdr.threshold = 0.1`, we cover all the true positives but get 13 false positives.

``` r
TestFDRThresholding(number.of.signals=100, number.of.noises=10000, fdr.threshold=0.1)
```

    ## [1] "# signals: 100"
    ## [1] "# noises: 10000"
    ## [1] "# FDR threshold: 0.1"
    ## [1] "# significant results: 113"
    ## [1] "# true positives: 100 out of 100"
    ## [1] "# false positives: 13"

![](autoPValThresholding_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-2-1.png)

If we lower the threshold down to `fdr.threshold = 0.01`, we still cover all the true positives **but no false positives**.

``` r
TestFDRThresholding(number.of.signals=100, number.of.noises=10000, fdr.threshold=0.01)
```

    ## [1] "# signals: 100"
    ## [1] "# noises: 10000"
    ## [1] "# FDR threshold: 0.01"
    ## [1] "# significant results: 100"
    ## [1] "# true positives: 100 out of 100"
    ## [1] "# false positives: 0"

![](autoPValThresholding_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png)

But what if we go stricter by setting `fdr.threshold = 0.01`? This time, we miss 1 true positive.

``` r
TestFDRThresholding(number.of.signals=100, number.of.noises=10000, fdr.threshold=0.001)
```

    ## [1] "# signals: 100"
    ## [1] "# noises: 10000"
    ## [1] "# FDR threshold: 0.001"
    ## [1] "# significant results: 99"
    ## [1] "# true positives: 99 out of 100"
    ## [1] "# false positives: 0"

![](autoPValThresholding_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-4-1.png)

The problem gets worse as you are dealing with more tests. If we have 10000 signals instead of noise, setting `fdr.threshold = 0.1` **would bring in 521, instead of 10 false positives**.

``` r
TestFDRThresholding(number.of.signals=10000, number.of.noises=10000, fdr.threshold=0.1)
```

    ## [1] "# signals: 10000"
    ## [1] "# noises: 10000"
    ## [1] "# FDR threshold: 0.1"
    ## [1] "# significant results: 10493"
    ## [1] "# true positives: 10000 out of 10000"
    ## [1] "# false positives: 493"

![](autoPValThresholding_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-5-1.png)

The automatic way
-----------------

TODO: optimal threshold ==
- Maximimises the number of true positives whilst minimising the number of false positives
- the point where the uniform distribution (corresponding to the P-values of the noises) becomes a spike (corresponding to a mixture of the P-values of the signals and the noises)
