# Code accompanying the paper 'Embryonal precursors of Wilms tumor' (Coorens et al., 2019)

## General prep of SNV calls

For this paper, SNVs were called using CaVEMan against the reference genome (GRCh37). 
After standard CaVEMan filtering and flagging (ASMD>=140, CLPM=0), all mutation calls per patient were collated and recounted across all samples from the same patient using cgpVAF/AlleleCounter. 
This will give you a matrix of total depth (NR, muts x samples) and variants reads (NV), which are used as a basis for further analysis.

## Filtering out germline variants

We fitted a binomial distribution to the combined read counts of all normal samples from one patient per SNV site, with the total depth as the number of trials, and the total number of reads supporting the variant as number of successes.  
Germline and somatic variants were differentiated based on a one-sided exact binomial test. 
For this test, the null hypothesis is that the number of reads supporting the variants across copy number normal samples is drawn from a binomial with p=0.5 (p=0.95 for copy number equal to one), and the alternative hypothesis drawn from a distribution with p<0.5 (or p<0.95). 
Resulting p-values were corrected for multiple testing with the Benjamini-Hochberg method and a cut-off was set at q < 10-5 to minimize false positives as on average, roughly 40,000 variants were subjected to this statistical test. 
Variants for which the null hypothesis could be rejected were classified as somatic, otherwise as germline. 
When parental genomes were sequenced, de novo germline mutations were taken to be those variants classified as germline in the child, but absent in both parents. 
In addition, it's very advisable to look at the average depth per variant across the normal samples. 
The resulting distribution should roughly look normal; any outliers should be filtered out. 
These might come from repeat regions in the germline and wouldn't be picked up by the exact binomial. 
Be careful to split autosomal and XY-chromosomal variants in male samples.

The code for this exact binomial test can be found in germline_exact_binom.R. 

## Filtering out noise/testing for true presence and absence by shearwater-like filter

The remaining somatic variants were recounted using AlleleCounter across all the samples in this study. 
For each patient, the non-tumor samples in this study not belonging to that patient were used as a reference panel to obtain the locus-specific error rate. 
Presence of the variant in the sample was accepted if the multiple-testing corrected p-value was less than 0.001, again, to minimize the false positive rate. 
In effect, this will look at the mutation rate of the background population (representing sequencing noise/artefact etc) and test whether the observed NV/NR is likely to come from the background distribution. 
Presence of the variant in the sample was accepted if the multiple-testing corrected p-value was less than 0.001, again, to minimize the false positive rate. 

The code for this filter can be found in shearwater_like_filter.R.

## Filtering out possible tumour contamination by binomial mixture modelling

The possibility of low-level tumor contamination in the non-tumor samples was excluded with a binomial mixture model on the read counts of all variants in the tumor trunk in non-tumor samples. The underlying model is that the reads supporting tumor variants in normal tissues either come from contamination or from shared development. This mixture can be deconvoluted because the underlying distributions would differ in proportion and VAF (probability of success in the binomial distribution). For this reason, the binomially distributed read counts are best represented by a one-dimensional binomial mixture model.

The optimal number, proportion, and locations of the mixtures were determined using an expectation-maximization algorithm. A range of cluster numbers (1:3) was used in this algorithm, and the optimal was chosen using the Bayesian Information Criterion (BIC). Variants belonging to the lowest peak, which corresponds to variants coming from noise or contamination, were excluded from supporting clonal nephrogenesis.

VAF peaks for the LCM WGS samples were computed via a binomial mixture model as well, but in this case the binomial distributions used for fitting were truncated. This reflects that no variant will be called with fewer than four supporting reads and renormalises the binomial distribution accordingly. 

The code for this can be found in binom_mix_model.R

## Code for plotting figures

The code to plot the figures in the paper can be found in plot_figures.R This script requires you to download the supplementary material from the paper itself, as well as the LCM SNV calls and binomial mixtures in the .zip file under Data/LCM.

