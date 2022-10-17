# ADPenetrance: Penetrance calculation for autosomal dominant traits
###### *README updated 14/10/2022*
###### The repository is maintained by Thomas Spargo (<thomas.spargo@kcl.ac.uk>) - please reach out with any questions.

This README details:
* The contents of this repository
* Links to documentation for the R functions contained within

The approach presented here is described in full within the following publication ([1](https://doi.org/10.1101/2021.03.16.21253691)) and is also available within a [web app](https://adpenetrance.rosalind.kcl.ac.uk). The approach was developed on the basis of the disease model outlined here ([2](https://doi.org/10.1159/000330167)).

## Repository contents

###### adpenetrance_function.R
* Running this R script generates the function for this penetrance calculation approach, `adpenetrance`.
* Use of `adpenetrance` requires the subfunctions `adpenetrance.errorfit` and `adpenetrance.unadjusted` to be loaded - these are retrieved from GitHub when running the `adpenetrance_function.R` script.
* The `adpenetrance` approach can be applied in accordance with [this](https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance) documentation.

###### getResidualRisk.R
* Running this R script generates the `getResidualRisk` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/getResidualRisk).
* This is a simple function to facilitate calculation of the disease model parameter `g`, which can be passed to the `useG` argument of `adpenetrance` as an indication of disease risk for family members not inheriting the tested variant.

###### checkOnsetVariability.R
* Running this R script generates the `checkOnsetVariability` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/checkOnsetVariability).
* This function centres and overlays the age of onset distributions (as a density or cumulative density function) across two groups. It can be used to test whether there is similar variability in disease onset for people with and without a variant, which is an indication that penetrance estimates in age-dependent traits may be less affected by age of sampling.

###### simADPenetrance.R
* Running this R script generates the `simADPenetrance` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/simADPenetrance).
* simADPenetrance can be used to easily perform simulation studies to test possible effects of age of sampling upon penetrance estimation when there is unequal onset variability between groups (as indicated by the output of `checkOnsetVariability`).
* This function is dependent upon the main `adpenetrance` function and subfunctions, and subfunctions within the `subfunctions/` directory - each of these dependencies are retrieved from GitHub when running `simADPenetrance.R`.
* It is also dependent upon the R packages `ggplot2`, `plyr`, `reshape2`, which must be installed and loaded by the user.

#### subfunctions/


###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; adpenetrance_errorfit_function.R
* This script loads the `adpenetrance.errorfit` function which fits a polynomial regression model to predict error in unadjusted penetrance estimates made within `adpenetrance`.
* `adpenetrance.errorfit` will call the `adpenetrance.unadjusted` function in order to make unadjusted penetrance estimates for simulated variants of penetrance between 0 and 1 occuring within a simulated population. The polynomial regression model produced by `adpenetrance.errorfit` is fitted according to differences between the estimated and true penetrance values of each variant simulated.
* The unadjusted penetrance value obtained by `adpenetrance` is then adjusted by error predicted in this estimate under the fitted model to determine the final adjusted penetrance estimate.
* See [1](https://doi.org/10.1101/2021.03.16.21253691) and the ADPenetrance [documentation](https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance) for further details.

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; adpenetrance_unadjusted_function.R
* This script loads the `adpenetrance.unadjusted` function which performs steps 1-3 of the `adpenetrance` approach, providing unadjusted penetrance estimates.
* The function is to be called by both of the `adpenetrance` and `adpenetrance.errorfit` functions.
* Please see the the main ADPenetrance [documentation](https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance).

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; varChars.R
* This script loads the `varChars` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/varChars) which defines a matrix indicating disease onset characteristics of people with and without a variant.
* The function is to be called within the `simADPenetrance` function; and the output of the function is passed to `affAtAge` via `genFamily`.

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; affAtAge.R
* This script loads the `affAtAge` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/affAtAge) which pseudo-randomly assigns status of a person being affected (1) or not (0) by a disease based on age, variant status, and disease characteristics associated with variant status in a var.Chars matrix.
* The function is called recurrently by `genFamily` (within `simADPenetrance`) to generate assignments of disease affectedness for individual family members over time.
* It requires an output matrix from `varChars`.

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; genFamily.R
* This script loads the `genFamily` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/genFamily) which pseudo-randomly generates a vector detailing a family of disease states "Unaffected"/"Sporadic"/"Familial" across time (or at lifetime risk) according to various parameters and disease characteristics associated with the output of `varChars`.
* The function is called recurrently within `simADPenetrance` to generate populations of simulated families.
* It requires `affAtAge` and an output matrix from `varChars`.

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; bindRows.R
* This script loads the `bindRows` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/bindRows), which produces a matrix by combining a list of vectors or matrices with varying length/ncol, The final value/column of shorter list elements will be duplicated until max length/ncol is reached for all elements and rows can be bound.
* The function is called within `simADPenetrance` to combine output vectors from `genFamily`.

  
#### case_studies/

###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; case_data.csv
*  This contains the raw sample data for all case studies described in the publication associated with this repository ([1](https://doi.org/10.1101/2021.03.16.21253691)).


###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; published_case_studies.R
*  This script calls the `adpenetrance` and `getResidualRisk` functions and `case_data.csv` to estimate penetrance for case studies represented within `case_data.csv`.


#### approach_validation/

Please refer to the README documentation within the `approach_validation/` directory which describes the contents and steps taken for approach validation



___
  

#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691

2. Al-Chalabi, A. & Lewis, C. M. (2011). Modelling the Effects of Penetrance and Family Size on Rates of Sporadic and Familial Disease. *Human Heredity, 71*(4): 281-288. doi: 10.1159/000330167