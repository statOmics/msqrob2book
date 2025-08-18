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

#### To do

- **Stijn:** Check if we should replace the normalisation with the median
  centering with DESeq2-style normalization (= median of the log2
  ratios between the feature intensity in the sample and the geom mean
  in this feature over all samples). Test this in the msTrawler 
  spike-in (without their unrealistic normalisation approach) and the
  MsStatsTMT spike-in.
- **Stijn:** Adding omnibus test (cf Stijn and Lucas) for multiple contrasts,
  link with stageR. Koen showed that this improves statistical power
  when assessing interactions.
- **Stijn:** Unlock the weights slot, provided either as a vector (in function of
sample size) or a matrix.
- Officially release msqrobRefit()
- (omicsGMF for imputation)

## New functionality

Implement and demonstrate the new functionality. 

#### Integration of iSEE to explore DA results

iSEE provides a user-friendly interface to create custom panels that
generate a static plot. I hence created a custom overdispersion plot.
I created a wrapper around iSEE to select a set in a QFeatures object
and a model column so and combined this with iSEEu to use their
volcano plot. Remaining issue:

- The wrapper function performs the set selection as iSEE only except
  objects that inherit from SE. Hence, it is not possible to make a
  detail plot since it involves multiple sets (hence multiple SE 
  objects). The work around is to reimplement an iSEE-like app that
  take a QFeatures object, but this needs carefull design and may put
  some maintenance burden on us.
- The overdispersion plot is static, meaning it cannot be coloured, 
  shaped, etc with respect to other variables. The solution is to
  implement dedicated panel classes (use iSEEu panels as template).
  
#### Improved class management

- I implemented a few more accessor methods
- I think msqrob2's interface should be simplified. In my opinion
  these methods should not be dispatched:
  + getCoef
  + getContrast
  + getDF
  + getDFPosterior
  + and co ...
  + topTable
  + varContrast
  If these slot accessors are important, maybe provide a generic function, eg
  getSlot(object, slotName)
- I added new accessors which may be more useful to users:
  + getCoefNames (rethink name): provides the available parameter 
    names when building contrasts
- Added function to automatically generate all pairwise and one-vs-all
  contrasts
  
#### plotting functions

- overdispersionPlot()
- detailPlot()
- Dimension reduction: is scater functionality not sufficient? 