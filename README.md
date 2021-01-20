# Autosomal Dominant Penetrance Calculator
###### README *Updated 20/01/2021*
###### The repository is maintained by Thomas P Spargo <thomas.spargo@kcl.ac.uk>


This README file details:
* The contents of this repository
* Operation instructions for the penetrance calculator function available within

The approach presented here is described in full within this publication: ___

The method is also integrated into this web-app: adpenetrance.rosalind.kcl.ac.uk


## Repository contents

###### Penetrance_function.R
* This R script presents the penetrance calculator function.
* The function can be applied in accordance with the guidance below.
  
###### Published_examples.R
* This R script presents the function and the input values for the case examples which were published alongside the original description of this method.


***
  
## Function documentation  
#### Description
This R function can be used to calculate the pentrance of a genetic variant, or set of variants, which are pathogenic for an autosomal dominant phenotype. The calculation is based on family structures observed in association with the variant, characteristics of the implicated disease, and average family size for the region from which other data are drawn. Penetrance can be estimated with or without confidence intervals. Please see ___ for a more comprehensive description of the method.


  
#### Usage examples

Estimate penetrance *without* confidence intervals:  
`Pen_calculator(N, MF, MS, PF)`

Estimate penetrance *with* confidence intervals:  
`Pen_calculator(N, MF, MU, PA, PF, MF_SE, MU_SE, Zout)`
  
  
#### Arguments
`N` - sibship size (Provide data on average sibship size for the sample from which variant frequencies are estimated or assign an estimate for the region from which these data are drawn). Must be provided.

`MF` - Variant frequency in *familial* disease state. Specify alongside `MS` and/or `MU`. Do not specify `MA` if used.

`MS` - Variant frequency in *sporadic* disease state. Specify alongside `MF` and/or `MU`. Do not specify `MA` if used.

`MU` - Variant frequency in *unaffected* disease state. Specify alongside `MF` and/or `MS` OR with `MA`.

`MA` - Variant frequency in *affected* disease state. Specify alongside `MU`. Do not specify `MF` or `MS` if used.

`PA` - Probability of being affected in the population. Specify if values are given for `MA` and/or `MU`.

`PF` - Probability of being familial if affected (i.e. the disease familiality rate). Specify if values are given for `MF` and/or `MS`.
  
  
`MF_SE` - Standard error in `MF`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MF` and the SE estimates for each state with variant frequency data provided.

`MS_SE` - Standard error in `MS`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MS` and the SE estimates for each state with variant frequency data provided.

`MU_SE` - Standard error in `MU`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MU` and the SE estimates for each state with variant frequency data provided.

`MA_SE` - Standard error in `MA`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MA` and the SE estimates for each state with variant frequency data provided.

`Zout` - Specify Z value for deriving confidence intervals of output from the standard error that is derived. Defaults to 1.96, estimating 95% confidence intervals.
  
  
#### Details  
The approach is also described in the publication and supplementary materials of _____.

__**Function Input**__
To calculate penetrance, the user must indicate:  
* **Sibship size** (`N`).
* **Variant frequency** for either two or three of the familial (`MF`), sporadic (`MS`) and unaffected (`MU`) states **OR** for the affected (`MA`) and unaffected (`MU`) states.
* **Weighting factors** (`PA` and/or `PF`) for the variant frequency estimates.

**Sibship size** estimates can be given as as the average number of siblings reported in the sample from which variant frequencies have been drawn or they can be approximated from reported data available in global databases. In the original publication, we drew upon the [World Bank, World Development Indicators](https://databank.worldbank.org/reports.aspx?source=World-Development-Indicators) database, estimating `N` as the Total Fertility Rate for the regions from which variant frequency had been sampled.
  
**Variant frequency** estimates can be given for the following family disease structures ('disease states'):
* familial - 2 or more family members are affected
* sporadic - 1 family member is affected
* unaffected - 0 family members are affected
* affected -  1 or more family members are affected

Note that the familial and sporadic states are subsumed within the affected state; variant frequency estimates for the 'affected' state cannot be provided alongside estimates for either or both of the the 'familial' or 'sporadic' states.
  
**Weighting factors** must be given, but `PA` and `PF` may not both be necessary in all calculations:
* `PA` must be specified if `MA` and/or `MU` are provided.
* `PF` must be specified if `MF` and/or `MS` are provided.

To calculate confidence intervals for the penetrance estimate, the user must indicate the **standard error** of the variant frequency estimates provided. These are the specified in the arguments `MF_SE`, `MS_SE`,`MU_SE`, `MA_SE` and should be given in all the states for which variant frequency estimates were provided.
  
  
__**Output**__
The results will be returned within two matrices

**Matrix 1**: `$Calculated` presents the rate at which one of the disease states represented within the input data occurred in people harbouring the variant modelled and the penetrance estimate to which this corresponds (at the defined sibship size).

The disease state rate result is presented in two forms:
1. the 'observed' rate of this state.
2. the 'expected' rate of this state.

The 'observed' rate is calculated as a weighted proportion of the variant frequencies defined in the user-input data.
The 'expected' rate is one of a series of disease rates stored in a lookup table generated within the function. Each of these has been derived as per the model equations (___) for a particular value of penetrance between 0 and 1 specified at intervals of .0001 and the sibship size `N`. The expected rate value in the output represents the value in the lookup table which most closely mached the observed rate. The penetrance estimate given in this output corresponds directly to the 'expected' rate value.

Note that the expected rate given should approximately equal the observed rate. An exception to this would be if penetrance is estimated to be 0 or 1. In this scenario, the observed and expected rates may deviate, as the observed rate could be lesser than or exceed the rate expected respectively at penetrance values of 0 and 1.

If error propagation is performed, the output of this matrix will also include estimates of disease state rate and penetrance at the confidence interval bounds. The standard error in the 'observed' rate estimate will also be given.
  
  
**Matrix 2**: `$Probabilities` shows the rates at which the familial, sporadic and unaffected states would be expected to occur in the population at the penetrance estimate yielded. If error propagation is performed, the output of this matrix will also include these estimates at the confidence interval bounds of penetrance.
