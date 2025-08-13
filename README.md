# msqrob2universe_paper

Repository with code to reproduce the manuscript for the
msqrob2universe paper. 

#### Data sets

*Spike-in*

- CPTAC: LFQ
- E Coli: LFQ (cf Adriaan's paper) **Stijn check email**
- MSstatsTMT: TMT (cf msqrob2TMT paper)
- DIA? (ask Ann Staes? or Leander's data?)

*Case study*

- Plubell mouse diet: TMT (cf msrob2TMT paper)
- Mouse macrophage: LFQ (cf Emmy's paper)
- DIA? (ask Ann Staes?)

#### To test

- **Stijn:** Check if we should replace the normalisation with the median
  centering with DESeq2-style normalization (= median of the log2
  ratios between the feature intensity in the sample and the geom mean
  in this feature over all samples). Test this in the msTrawler 
  spike-in (without their unrealistic normalisation approach) and the
  MsStatsTMT spike-in.
- Integration of iSEE to explore DA results

#### Demonstrate functionality

- **Stijn:** Adding omnibus test (cf Stijn and Lucas) for multiple contrasts,
  link with stageR. Koen showed that this improves statistical power
  when assessing interactions.
- Automatically create contrasts for all paiwrse comparisons for a
  given factor. Similarly, automatically generate contrast that compares one level against the
  average of all the other levels. 
- Implement plot PCA (also omicsGMF, NIPALS), detail plots, trended
  overdispersion (variance against covariate, eg .n)
- **Stijn:** Unlock the weights slot, provided either as a vector (in function of
sample size) or a matrix.
- Add missing values when low peptide counts for a protein in a
  sample.
- (omicsGMF for imputation)
