# ADPenetrance: Penetrance calculation for autosomal dominant traits
###### *README Updated 01/10/2021*
###### The repository is maintained by Thomas Spargo (<thomas.spargo@kcl.ac.uk>) - please get in touch if you encounter any issues with the contents

This README details:
* The contents of this repository
* Operation instructions for the `adpenetrance` R function, for estimating penetrance in autosomal dominant traits, available within

The approach presented here is described in full within the following publication ([1](https://doi.org/10.1101/2021.03.16.21253691)) and is also available within a [web app](https://adpenetrance.rosalind.kcl.ac.uk).


## Repository contents

###### adpenetrance_function.R
* Running this R script generates the function for this penetrance calculation approach, `adpenetrance`.
* Use of `adpenetrance` also requires the subfunctions `adpenetrance.errorfit` and `adpenetrance.unadjusted` to be loaded - these are retrieved from GitHub when running the `adpenetrance_function.R` script.
* The `adpenetrance` approach can be applied in accordance with the guidance below.
  
#### subfunctions/


###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; adpenetrance_errorfit_function.R
* This script generates the function `adpenetrance.errorfit` which fits a polynomial regression model to predict error in unadjusted penetrance estimates made within `adpenetrance`.
* `adpenetrance.errorfit` will call the `adpenetrance.unadjusted` function in order to make unadjusted penetrance estimates for simulated variants of penetrance between 0 and 1 occuring within a simulated population. The polynomial regression model produced by `adpenetrance.errorfit` is fitted according to differences between the estimated and true penetrance values of each variant simulated.
* The unadjusted penetrance value obtained by `adpenetrance` is then adjusted by error predicted in this estimate under the fitted model to determine the final adjusted penetrance estimate.
* See [1](https://doi.org/10.1101/2021.03.16.21253691) for further details.


###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; adpenetrance_unadjusted_function.R
* This script generates a reduced version of the `adpenetrance` function which generates only unadjusted penetrance estimates.
* It is to be called by the `adpenetrance.errorfit` function to calculate error expected in unadjusted penetrance estimates across values of true penetrance as modelled in a simulated population defined within `adpenetrance.errorfit`.


  
#### case_studies/

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; case_data.csv
*  This contains the raw sample data for all case studies described in the publication associated with this repository ([1](https://doi.org/10.1101/2021.03.16.21253691)).


###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; published_case_studies.R
*  This script calls `adpenetrance_function.R` and `case_data.csv` to estimate penetrance for case studies represented within `case_data.csv`.


#### approach_validation/

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; adpenetrance_MLE_function.R
* This script generates the `adpenetrance.MLE` function, which adapts the primary `adpenetrance` function to produce additional unadjusted penetrance estimates via a maximum likelihood approach. These can be compared to the unadjusted penetrance estimates generated using the lookup table approach otherwise employed.

* The function is called within the `MLE_validation.R` script.
    

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; MLE_validation.R

* This script calls `adpenetrance_MLE_function.R` and `case_data.csv` to make penetrance estimates for those case studies presented in Table 2 of the associated manuscript ([1](https://doi.org/10.1101/2021.03.16.21253691)).

* This re-analysis compares unadjusted penetrance estimates obtained from the lookup table approach employed within `adpenetrance` to those obtained via a maximum likelihood method.

* Negligable differences are observed between the results of the two approaches.


##### approach_validation/simulation_studies/

* The following scripts each test the accuracy of adjusted penetrance estimates produced by `adpenetrance` in datasets of  simulated families harbouring hypothetical variants for which the true penetrance is known.

* Two distinct populations of simulated families are generated in these simulations. The first is based on the sibship distribution observed in a UK population 1974 birth cohort ([2](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable)). The second is based on the distribtion the 'Next Steps' cohort study of English families ([3](https://doi.org/10.1017/ehs.2020.54)).

* The accuracy of penetrance estimates made for hypothetical variants occuring in these two populations is compared according to the states modelled (see `states`) and whether the user supplies `adpenetrance` with information about the distirbution of sibships in their dataset (see `define_sibstructure`).

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; simulation_1_correct_parameter_specification.R

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 20 ground truth penetrance values when input parameters for the function (`N` and `RX`) are correctly specified.

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; simulation_2_N_parameter_incorrect.R

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 5 ground truth penetrance values when the mean sibship size parameter (`N`) is incorrectly specified, testing different degrees of misspecification in `N`.

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; simulation_3_obsRX_parameter_incorrect.R

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 5 ground truth penetrance values when the observed rate of state X parameter (`RX`) is incorrectly specified, testing different degrees of misspecification in `RX`.





***
  
## Function documentation  
#### Description
The `adpenetrance` R function can be used to calculate the pentrance of a germline genetic variant, or set of variants, which are pathogenic for an autosomal dominant phenotype. Penetrance can be estimated with or without confidence intervals. The calculation is based on the rate at which one family disease structure ('disease state') occurs across a valid subset of disease states in people who harbour the assessed variant and on the average sibship size of those people sampled these data. 

The approach considers the following disease states:
* familial - two or more first-degree family members are affected
* sporadic - one family member is affected (none of their first-degree relatives are affected)
* unaffected - no family members are affected (a control population)
* affected - one or more first-degree family members are affected 

To operate the approach, input data must be include any two or three of the familial, sporadic and unaffected states **OR** for the affected and unaffected states.

The approach is further described in the details section below, and a comprehensive outline of the method is given in ([1](https://doi.org/10.1101/2021.03.16.21253691)).


  
#### Usage examples

Estimate penetrance using input data from the familial and sporadic states *without* confidence intervals:  
`adpenetrance(N, MF, MS, PF)`

Estimate penetrance using input data from the familial and unaffected states *with* confidence intervals:  
`adpenetrance(N, MF, MU, PA, PF, MF_SE, MU_SE, Zout)`

Estimate penetrance using input data from the familial and sporadic states *without* confidence intervals, and supply the method with information about the distribution of sibship sizes across the sample data using a 2 column matrix:  
`adpenetrance(N, MF, MS, PF, define_sibstructure=matrix(c(0:4,c(0.18, 0.18, 0.37, 0.16, 0.11)),ncol=2))`

Estimate penetrance according to the rate of familial disease across people harbouring the variant sampled from the familial and sporadic disease states (*without* confidence intervals):
`adpenetrance(N, RX, states="fs")`

Estimate penetrance according to the rate of familial disease across people harbouring the variant sampled from the familial and sporadic disease states (*with* confidence intervals):
`adpenetrance(N, RX, RX_SE, Zout, states="fs")`



  
  
#### Arguments
`N` - sibship size (define the average sibship size for the sample from which variant characteristics are estimated or assign an estimate representative of this sample). Must be provided.

`MF` - Variant frequency in *familial* disease state. Specify alongside `MS` and/or `MU`. Do not specify `MA` or `RX` if used.

`MS` - Variant frequency in *sporadic* disease state. Specify alongside `MF` and/or `MU`. Do not specify `MA` or `RX` if used.

`MU` - Variant frequency in *unaffected* disease state. Specify alongside `MF` and/or `MS` OR with `MA`. Do not specify `RX` if used.

`MA` - Variant frequency in *affected* disease state. Specify alongside `MU`. Do not specify `MF`,  `MS`, or `RX` if used.

`PA` - Probability of a person from the sampled population of being affected. Specify if values are given for `MA` and/or `MU`.

`PF` - Probability of being familial if affected (i.e. the disease first-degree familiality rate). Specify if values are given for `MF` and/or `MS`.
  
  
`MF_SE` - Standard error in `MF`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MF` and the SE estimates for each state with variant frequency data provided.

`MS_SE` - Standard error in `MS`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MS` and the SE estimates for each state with variant frequency data provided.

`MU_SE` - Standard error in `MU`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MU` and the SE estimates for each state with variant frequency data provided.

`MA_SE` - Standard error in `MA`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `MA` and the SE estimates for each state with variant frequency data provided.

`Zout` - Specify Z value for deriving confidence intervals of output from the calculated standard error. Defaults to 1.96, estimating 95% confidence intervals.


`RX` - Specify the rate of 'state X' among people harbouring the tested variant sampled from a valid set of disease states. State X can be either familial, sporadic or affected (see details below). Must also specify the `states` term where `RX` is given. Do not specify any of `MF`, `MS`, `MU`, or `MA` if used.

`RX_SE` - Standard error in `RX`. Used to calculate confidence intervals of penetrance estimate. Specify alongside `RX`.

`states` - Indicates which states are represented within the `RX` calculation. Is a string variable and can be defined as: `"fsu"`,`"fs"`,`"fu"`,`"su"`,`"au"` (see details below). Must be provided where `RX` is defined.

`define_sibstructure` - Optionally supply either a vector detailing the sibship sizes of all sampled families or a summary of the sibship distribution of the sample (see details below). Passed to `adpenetrance.errorfit` subfunction.

`include_MLE` - Logical (defaults to `TRUE`). For `adpenetrance.MLE` only, indicate whether or not to make additional unadjusted penetrance estimates via a maximum likelihood approach. If `FALSE`, `adpenetrance.MLE` functions equivalently to `adpenetrance`.

  
  
#### Details  
The approach is comprehensively described in the manuscript of ([1](https://doi.org/10.1101/2021.03.16.21253691)) and its supplementary materials.

##### Function Input

Input data for penetance calculation can be given as one of two main structures.

In __both data structures__:

* *Sibship size* (`N`) must be indicated. This should represent the average size of sibships across the samples used to define variant characteristics. It can  be estimated for the sample either directly, based on the average sibship size among the described families, or indirectly, by designating an estimate representative of the sampled population (e.g. available within global databases). In the original publication, we drew upon the [World Bank, World Development Indicators](https://databank.worldbank.org/reports.aspx?source=World-Development-Indicators) database, approximating `N` as the Total Fertility Rate of the regions from which variant frequencies were estimated.

* The user can optionally indicate the *distribution of sibship sizes across sampled families* (`define_sibstructure`). This should be either a vector of integers containing the sibship sizes of each sampled family or a 2 column matrix or data frame where column 1 details the sibship sizes occuring in the sample and column 2 details the sample proportion to which each sib-size corresponds. This sib-structure information is passed to `adpenetrance.errorfit` and is used to tailor the sibship distribution of an internally-simulated population; where no information is passed to `define_sibstructure`, sibships in the simulated dataset follow a Poisson distribution, with lambda defined by `N`. This simulated population is used to fit a polynomial regression model to predict the difference between an unadjusted penetrance estimate and the true penetrance value. The unadjusted penetrance estimate made by `adpenetrance` is then adjusted by this predicted error to derive the adjusted penetrance estimate (See [1](https://doi.org/10.1101/2021.03.16.21253691) for further details). Supplying accurate information to `define_sibstructure` will allow for a more precise adjustment of the penetrance estimate.

Within __data structure 1__, the user must also indicate:  
* *Variant frequency* within populations screened for the variant, drawn from either two or three of the familial (`MF`), sporadic (`MS`) and unaffected (`MU`) states **OR** for the affected (`MA`) and unaffected (`MU`) states.
* *Weighting factors* (`PA` and/or `PF`) for the variant frequency estimates.

Within __data structure 2__, the user must also indicate:  
* *Rate of state X* (`RX`), which is either familial, sporadic, or affected among people harbouring the assessed variant, drawn from either two or three of the familial, sporadic and unaffected states **OR** for the affected and unaffected states as indicated within the `states` argument. State X is always the first state indicated within the states argument.
* *States* (`states`) that are included for calculating `RX`.


##### Information for data structure 1
*Variant frequency* estimates can be given for the following family disease structures ('disease states'):
* familial (`MF`) - two or more first-degree family members are affected
* sporadic (`MS`) - one family member is affected (none of their first-degree relatives are affected)
* unaffected (`MU`) - no family members are affected (a control population)
* affected (`MA`) - one or more first-degree family members are affected 

Note that the familial and sporadic states are subsumed within the affected state; variant frequency estimates for the 'affected' state cannot be provided alongside estimates for either or both of the the 'familial' or 'sporadic' states.
  
*Weighting factors* must be given, but `PA` and `PF` are not both necessary in all disease-state combinations:
* `PA` must be specified if `MA` and/or `MU` are provided.
* `PF` must be specified if `MF` and/or `MS` are provided.

To additionally calculate confidence intervals for the penetrance estimate, the user must indicate the *standard error* of all the variant frequency estimates provided. These are the specified in the arguments `MF_SE`, `MS_SE`, `MU_SE`, `MA_SE` and should be given in all the states for which variant frequency estimates were provided. The `Zout` argument defines the level of confidence to be estimated, defaulting to a value of `1.96` which will give 95% confidence.



##### Information for data structure 2
The value of the `RX` argument indicates the rate at which 'state X' occurs across a valid set of disease states among people harbouring the assessed variant. The value of the `states` argument indicates from which states people have been sampled and which state is considered to be 'X' within the function:
* `"fsu"` - people are sampled from the familial, sporadic and unaffected states and state X is familial
* `"fs"` - people are sampled from the familial and sporadic states and state X is familial
* `"fu"` - people are sampled from the familial, unaffected states and state X is familial
* `"su"` - people are sampled from the sporadic and unaffected states and state X is sporadic
* `"au"` - people are sampled from the affected and unaffected states and state X is affected

The user should specify which states have been sampled in the `states` argument and the value of `RX` for this state combination (e.g. if `states="fs"`, `RX` is the rate of familial disease across people harbouring the assessed variant sampled across the familial and sporadic states).

To calculate confidence intervals for the penetrance estimate, the user must indicate the standard error of `RX`, given in the `RX_SE` argument. The `Zout` argument defines the level of confidence to be estimated, defaulting to a value of `1.96` which will give 95% confidence.


  
##### Output
The results will be returned as a matrix: `$output`.

This presents the rate at which 'state X' occurred across all the states modelled within the input data for people harbouring the assessed variant, the unadjusted penetrance estimate to which this rate corresponds (at the defined `N`), and the adjusted penetrance estimate after correcting for systematic bias within the unadjusted estimate.

The disease state rate result is presented in two forms:
1. the 'observed' rate 
2. the 'expected' rate 

If data are specified using **structure 1**, then the 'observed' rate has been calculated as a weighted proportion of the variant frequencies defined in the input data. If data are specified using **structure 2** then the 'observed' rate is the value defined as `RX`.

The 'expected' rate is one of a series of 'state X' rates stored in a lookup table that is generated within the function. These are derived as per the disease model equations followed within this method ([1](https://doi.org/10.1101/2021.03.16.21253691), [4](https://doi.org/10.1159/000330167)) for values of penetrance between 0 and 1 at intervals of 0.0001 and the sibship size `N`. The expected rate shown in the output represents the value in the lookup table which most closely maches the observed rate.

The unadjusted penetrance estimate is obtained from the ookup table and corresponds directly to the 'expected' disease state rate. The adjusted penetrance estimate is derived from the unadjusted estimate and error predicted in this estimate under an nth degree polynomial regression model which is fitted by `adpenetrance.errorfit`: `adjusted penetrance = unadjusted penetrance + predicted error`. Refer to [1](https://doi.org/10.1101/2021.03.16.21253691) for further detail on fitting this model.

Note that the expected disease state rate should approximately equal the observed rate. An exception to this would be if the unadjusted penetrance is estimated to be 0 or 1. In this scenario, the observed and expected rates may deviate, as the observed rate could be less than or exceed the rate expected respectively at penetrance values of 0 and 1.

If data are given to allow estimation of error in the penetrance estimate, then the `output` matrix will include estimates of disease state rate and penetrance at the confidence interval bounds. The standard error in the 'observed' rate will also be given.

#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691

2. Office for National Statistics. Childbearing for women born in different years. 2020. Available from: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable

3. Sheppard P, Monden C. When does family size matter? Sibship size, socioeconomic status and education in England. Evol Hum Sci. 2020; 2, e51:1-21. doi: 10.1017/ehs.2020.54

4. Al-Chalabi, A. & Lewis, C. M. (2011). Modelling the Effects of Penetrance and Family Size on Rates of Sporadic and Familial Disease. *Human Heredity, 71*(4): 281-288. doi: 10.1159/000330167