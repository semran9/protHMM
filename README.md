
<!-- README.md is generated from README.Rmd. Please edit that file -->

# protHMM

<!-- badges: start -->
![R CMD check](https://github.com/semran9/protHMM/actions/workflows/R-CMD-check.yaml/badge.svg)

<!-- badges: end --> 

## Summary

The goal of protHMM is to help integrate profile hidden markov model
(HMM) representations of proteins into the machine learning and
bioinformatics workflow. protHMM ports a number of features from use in
Position Specific Scoring Matrices (PSSMs) to HMMs, along with
implementing features used with HMMs specifically, which to our
knowledge has not been done before. The adoption of HMM representations
of proteins derived from HHblits and HMMer also presents an opportunity
for innovation; it has been shown that HMMs can benefit from better
multiple sequence alignment than PSSMs and thus get better results than
corresponding HMMs using similar feature extraction techniques (Lyons et
al. 2015). protHMM implements 20 different feature extraction techniques
to provide a comprehensive list of feature sets for use in
bioinformatics tasks ranging from protein fold classification to
protein-protein interaction.

#### Installation

You can install the development version of protHMM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("semran9/protHMM")
```

## Functions List

hmm_ac()

hmm_bigrams()

hmm_cc()

chmm()

hmm_distance()

fp_hmm()

hmm_GA()

hmm_GSD()

IM_psehmm()

hmm_LBP()

hmm_LPC()

hmm_MA()

hmm_MB()

pse_hmm()

hmm_read()

hmm_SCSH()

hmm_SepDim()

hmm_Sigle_Average()

hmm_smooth()

hmm_svd()

hmm_trigrams()

## Example

``` r
## this shows the functionality of hmm_distance, which calculates a similarity score between two proteins
## other functions are documented fully in the protHMM vignette
library(protHMM)
## these proteins are from the same fold and similar; h should be low
h <- hmm_distance(system.file("extdata", "1DLHA2-7", package="protHMM"), system.file("extdata", "1TEN-7", package="protHMM"))
## these proteins are from different folds and not similar; h_2 should be high
h_2<- hmm_distance(system.file("extdata", "1DLHA2-7", package="protHMM"), system.file("extdata", "1TAHA-23", package="protHMM"))
h < h_2
#> [1] TRUE
```

## References

Lyons, J., Dehzangi, A., Heffernan, R., Yang, Y., Zhou, Y., Sharma, A.,
& Paliwal, K. K. (2015). Advancing the Accuracy of Protein Fold
Recognition by Utilizing Profiles From Hidden Markov Models. IEEE
Transactions on Nanobioscience, 14(7), 761–772.
