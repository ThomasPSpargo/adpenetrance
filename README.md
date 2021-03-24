# ADPenetrance: Penetrance calculation for autosomal dominant traits
###### *README Updated 18/03/2021*
###### The repository is maintained by Thomas P Spargo (<thomas.spargo@kcl.ac.uk>)


This README file details:
* The contents of this repository
* Operation instructions for the penetrance calculator function available within

The approach presented here is described in full within the following publication ([1](https://doi.org/10.1101/2021.03.16.21253691)) and available within a [web app](https://adpenetrance.rosalind.kcl.ac.uk).


## Repository contents

###### Penetrance_function.R
* Running this R script generates a function for this penetrance calculation approach.
* It can be applied in accordance with the guidance below.
  
###### Published_examples.R
* This R script presents the input values for the case studies described within the publication within which this method was first outlined ([1](https://doi.org/10.1101/2021.03.16.21253691)).
* This script calls the Penetrance_function.R script in order to generate the penetrance calculator function and will not operate without additionally downloading this file.


***
  
## Function documentation  
#### Description
This R function can be used to calculate the pentrance of a germline genetic variant, or set of variants, which are pathogenic for an autosomal dominant phenotype. Penetrance can be estimated with or without confidence intervals. The calculation is based on the rate at which one family disease structure ('disease state') occurs across a valid subset of disease states in people who harbour the assessed variant and on the average sibship size of those people sampled these data. 

The approach considers the following disease states:
* familial - two or more members of a family are affected
* sporadic - one member of a family is affected
* unaffected - no family members are affected (a control population)
* affected - one or more members of a family are affected 

To operate the approach, input data must be include any two or three of the familial, sporadic and unaffected states **OR** for the affected and unaffected states.

The approach is further described in the details section below, and a comprehensive outline of the method is given in [(1)](https://doi.org/10.1101/2021.03.16.21253691).


  
#### Usage examples

Estimate penetrance using input data from the familial and sporadic states *without* confidence intervals:  
`adpenetrance(N, MF, MS, PF)`

Estimate penetrance using input data from the familial and unaffected states *with* confidence intervals:  
`adpenetrance(N, MF, MU, PA, PF, MF_SE, MU_SE, Zout)`

Estimate penetrance according to the rate of familial disease across people harbouring the variant sampled from the familial and sporadic disease states (*without* confidence intervals):
`adpenetrance(N, RX, states="fs")`

Estimate penetrance according to the rate of familial disease across people harbouring the variant sampled from the familial and sporadic disease states (*with* confidence intervals):
`adpenetrance(N, RX, RX_SE, Zout, states="fs")`

  
  
#### Arguments
`N` - sibship size (Provide data on average sibship size for the sample from which variant characteristics are estimated or assign an estimate for the region from which these data are sampled). Must be provided.

`MF` - Variant frequency in *familial* disease state. Specify alongside `MS` and/or `MU`. Do not specify `MA` or `RX` if used.

`MS` - Variant frequency in *sporadic* disease state. Specify alongside `MF` and/or `MU`. Do not specify `MA` or `RX` if used.

`MU` - Variant frequency in *unaffected* disease state. Specify alongside `MF` and/or `MS` OR with `MA`. Do not specify `RX` if used.

`MA` - Variant frequency in *affected* disease state. Specify alongside `MU`. Do not specify `MF`,  `MS`, or `RX` if used.

`PA` - Probability of a person from the sampled population of being affected. Specify if values are given for `MA` and/or `MU`.

`PF` - Probability of being familial if affected (i.e. the disease familiality rate). Specify if values are given for `MF` and/or `MS`.
  
  
`MF_SE` - Standard error in `MF`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MF` and the SE estimates for each state with variant frequency data provided.

`MS_SE` - Standard error in `MS`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MS` and the SE estimates for each state with variant frequency data provided.

`MU_SE` - Standard error in `MU`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MU` and the SE estimates for each state with variant frequency data provided.

`MA_SE` - Standard error in `MA`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MA` and the SE estimates for each state with variant frequency data provided.

`Zout` - Specify Z value for deriving confidence intervals of output from the standard error that is derived. Defaults to 1.96, estimating 95% confidence intervals.


`RX` - Specify the rate of 'state X' among people harbouring the tested variant sampled from a valid set of disease states. State X can be either familial, sporadic or affected (see details below). Must also specify the `states` term where `RX` is given. Do not specify any of `MF`, `MS`, `MU`, or `MA` if used.

`RX_SE` - Standard error in `RX`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `RX`.

`states` - Indicates which states are represented within the `RX` calculation. Is a string variable and can be defined as: `"fsu"`,`"fs"`,`"fu"`,`"su"`,`"au"` (see details below). Must be provided where `RX` is defined.
  
  
#### Details  
The approach is also described in the manuscript of ([1](https://doi.org/10.1101/2021.03.16.21253691)) and it's supplementary materials.

#####Function Input

Input data for penetance calculation can be given as one of two main structures.
Within __data structure 1__, the user must indicate:  
* *Sibship size* (`N`).
* *Variant frequency* within populations screened for the variant, drawn from either two or three of the familial (`MF`), sporadic (`MS`) and unaffected (`MU`) states **OR** for the affected (`MA`) and unaffected (`MU`) states.
* *Weighting factors* (`PA` and/or `PF`) for the variant frequency estimates.

Within __data structure 2__, the user must indicate:  
* *Sibship size* (`N`).
* *Rate of state X* (`RX`), which is either familial, sporadic, or affected among people harbouring the assessed variant, drawn from either two or three of the familial, sporadic and unaffected states **OR** for the affected and unaffected states as indicated within the `states` term.
* *States* (`states`) that are included for calculating `RX`.

`N` is estimated equivalently across each data structure. This should represent the average size of sibships across the samples used to define variant characteristics. It can  be estimated for the sample either directly, based on the average sibship size among the described families, or indirectly, by designating an estimate representative of the sampled population (e.g. available within global databases). In the original publication, we drew upon the [World Bank, World Development Indicators](https://databank.worldbank.org/reports.aspx?source=World-Development-Indicators) database, approximating `N` as the Total Fertility Rate of the regions from which variant frequencies were estimated. 


#####Information for data structure 1
*Variant frequency* estimates can be given for the following family disease structures ('disease states'):
* familial (`MF`) - two or more family members are affected
* sporadic (`MS`) - one family member is affected
* unaffected (`MU`) - no family members are affected
* affected (`MA`) -  one or more family members are affected

Note that the familial and sporadic states are subsumed within the affected state; variant frequency estimates for the 'affected' state cannot be provided alongside estimates for either or both of the the 'familial' or 'sporadic' states.
  
*Weighting factors* must be given, but `PA` and `PF` are not both necessary in all disease-state combinations:
* `PA` must be specified if `MA` and/or `MU` are provided.
* `PF` must be specified if `MF` and/or `MS` are provided.

To additionally calculate confidence intervals for the penetrance estimate, the user must indicate the *standard error* of all the variant frequency estimates provided. These are the specified in the arguments `MF_SE`, `MS_SE`,`MU_SE`, `MA_SE` and should be given in all the states for which variant frequency estimates were provided. The `Zout` argument defines the level of confidence to be estimated, defaulting to a value of `1.96` which will give 95% confidence.



#####Information for data structure 2
The value of the `RX` argument indicates the rate at which 'state X' occurs across a valid set of disease states among people harbouring the assessed variant. The value of the `states` argument indicates from which states people have been sampled and which state is considered to be 'X' within the function:
* `"fsu"` - people are sampled from the familial, sporadic and unaffected states and state X is familial
* `"fs"` - people are sampled from the familial and sporadic states and state X is familial
* `"fu"` - people are sampled from the familial, unaffected states and state X is familial
* `"su"` - people are sampled from the sporadic and unaffected states and state X is sporadic
* `"au"` - people are sampled from the affected and unaffected states and state X is affected

The user should specify which states have been sampled in the `states` argument and the value of `RX` for this state combination (e.g. `states="fs", RX=0.5`, where `RX` is the rate of familial disease across the familial and sporadic states among people harbouring the assessed variant).

To calculate confidence intervals for the penetrance estimate, the user must indicate the standard error of `RX`, given in the `RX_SE` argument. The `Zout` argument defines the level of confidence to be estimated, defaulting to a value of `1.96` which will give 95% confidence.



 


  
#####Output
The results will be returned across two matrices.

**Matrix 1**: `$Calculated` presents the rate at which 'state X' occurred across all the included states included within the input data for people harbouring the assessed variant and the penetrance estimate to which this rate corresponds (at the defined `N`).

The disease state rate result is presented in two forms:
1. the 'observed' rate of this state.
2. the 'expected' rate of this state.

If data are specified using **structure 1**, then the 'observed' rate has been calculated as a weighted proportion of the variant frequencies defined in the input data. If data are specified using **structure 2** then the 'observed' rate is the value defined as `RX`.

The 'expected' rate is one of a series of 'state X' rates stored in a lookup table that is generated within the function. These are derived as per the disease model equations followed within this method ([1](https://doi.org/10.1101/2021.03.16.21253691), [2](https://doi.org/10.1159/000330167)) for a values of penetrance between 0 and 1 at intervals of .0001 and the sibship size `N`. The expected rate shown in the output represents the value in the lookup table which most closely maches the observed rate. The penetrance estimate provided corresponds directly to the 'expected' rate shown.

Note that the expected rate given should approximately equal the observed rate. An exception to this would be if penetrance is estimated to be 0 or 1. In this scenario, the observed and expected rates may deviate, as the observed rate could be less than or exceed the rate expected respectively at penetrance values of 0 and 1.

If data are given to allow estimation of error in the penetrance estimate, then this matrix will include estimates of disease state rate and penetrance at the confidence interval bounds. The standard error in the 'observed' rate will also be given.
  
  
**Matrix 2**: `$Probabilities` displays (with or without confidence intervals) the probability that a family of sibship size `N` and harbouring a variant with the obtained penetrance estimate will present as familial, sporadic, or unaffected.


#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691
2. Al-Chalabi, A. & Lewis, C. M. (2011). Modelling the Effects of Penetrance and Family Size on Rates of Sporadic and Familial Disease. *Human Heredity, 71*(4): 281-288. doi: 10.1159/000330167