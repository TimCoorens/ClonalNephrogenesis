# Code accompanying the paper 'Embryonal precursors of childhood kidney cancer' (Coorens et al., 2019)

## General prep of SNV calls

For this paper, SNVs were called using CaVEMan against the reference genome (GRCh37). 
After standard CaVEMan filtering and flagging (ASMD>=140, CLPM=0), all mutation calls per patient were collated and recounted across all samples from the same patient using cgpVAF/AlleleCounter. 
This will give you a matrix of total depth (NR, muts x samples) and variants reads (NV), which are used as a basis for further analysis.

## Filtering out germline variants

Germline variants were then filtered on by an exact binomial test on the aggregate (rowSums) of NVs and NRs across the normal samples of a patient.
The code for this exact binomial test can be found in germline_exact_binom.R. 
In addition, it's very advisable to look at the average depth per variant across the normal samples. 
The resulting distribution should roughly look normal; any outliers should be filtered out. 
These might come from repeat regions in the germline and wouldn't be picked up by the exact binomial. 
Be careful to split autosomal and XY-chromosomal variants in male samples.

## Filtering out noise/testing for true presence and absence by shearwater-like filter

The remaining somatic variants were recounted using AlleleCounter across all the samples in this study. 
Per patient, all other patients then provide an unmatched normal panel to use in the statistical framework of Shearwater.
In effect, this will look at the mutation rate of the background population (representing sequencing noise/artefact etc) and test whether the observed NV/NR is likely to come from the background distribution.
If this is very unlikely (log(q)<(-3)), the variant classified as truly present in this sample.
The code for this filter can be found in shearwater_like_filter.R.

## Filtering out possible tumour contamination by binomial mixture modelling

The possibility of low-level tumor contamination in the non-tumor samples was excluded with a binomial mixture model on the read counts of all variants in the tumor trunk in non-tumor samples. The underlying model is that the reads supporting tumor variants in normal tissues either come from contamination or from shared development. This mixture can be deconvoluted because the underlying distributions would differ in proportion and VAF (probability of success in the binomial distribution). For this reason, the binomially distributed read counts are best represented by a one-dimensional binomial mixture model.

The optimal number, proportion, and locations of the mixtures were determined using an expectation-maximization algorithm. A range of cluster numbers (1:3) was used in this algorithm, and the optimal was chosen using the Bayesian Information Criterion (BIC). Variants belonging to the lowest peak, which corresponds to variants coming from noise or contamination, were excluded from supporting clonal nephrogenesis 

The code for this can be found in binom_mix_model.R
