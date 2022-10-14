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
*  This script calls the `adpenetrance` and `getResidualRisk` functions and `case_data.csv` to estimate penetrance for case studies represented within `case_data.csv`.


#### approach_validation/

Please refer to the README documentation within the `approach_validation/` directory which outlines the incluted contents and steps taken for approach validation



***
  

#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691

2. Al-Chalabi, A. & Lewis, C. M. (2011). Modelling the Effects of Penetrance and Family Size on Rates of Sporadic and Familial Disease. *Human Heredity, 71*(4): 281-288. doi: 10.1159/000330167