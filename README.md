# ADPenetrance: Penetrance calculation for autosomal dominant traits
___Updated 17/10/2022___

_The repository is maintained by Thomas Spargo (<thomas.spargo@kcl.ac.uk>) - please reach out with any questions._

This README outlines:
* The contents of this repository
* Links to documentation for the main R functions contained within

The approach followed here is described in full within the following publication ([1](https://doi.org/10.1101/2021.03.16.21253691)) and is also available within a [web app](https://adpenetrance.rosalind.kcl.ac.uk). ADPenetrance was developed on the basis of the disease model outlined here ([2](https://doi.org/10.1159/000330167)).

## Repository contents

__adpenetrance_function.R__
* Running this R script generates the function for this penetrance calculation approach, `adpenetrance`.
* Use of `adpenetrance` requires the subfunctions `adpenetrance.errorfit` and `adpenetrance.unadjusted` to be loaded - these are retrieved from GitHub when running the `adpenetrance_function.R` script.
* `adpenetrance` can be applied in accordance with [this](https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance) documentation.

__getResidualRisk.R__
* Running this R script generates the `getResidualRisk` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/getResidualRisk).
* This is a simple function to facilitate calculation of the disease model parameter `g`, which can be passed to the `useG` argument of `adpenetrance` as an indication of disease risk for family members not harbouring the tested variant.

__checkOnsetVariability.R__
* Running this R script generates the `checkOnsetVariability` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/checkOnsetVariability).
* This function centres and overlays the age of onset distributions (as a density or cumulative density function) across two groups. It can be used to test whether there is similar variability in disease onset for people with and without a variant, which is an indication that penetrance estimates in age-dependent traits may be less affected by age of sampling.

__simADPenetrance.R__
* Running this R script generates the `simADPenetrance` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/simADPenetrance).
* simADPenetrance can be used to easily perform simulation studies to test possible effects of age of sampling upon penetrance estimation when there is unequal onset variability between groups (as indicated by the output of `checkOnsetVariability`).
* This function is dependent upon the main `adpenetrance` function and subfunctions, and subfunctions within the `subfunctions/` directory - each of these dependencies are automatically retrieved from GitHub when running `simADPenetrance.R`.
* It is also dependent upon the R packages `ggplot2`, `plyr`, `reshape2`, which must be installed and loaded by the user.

#### approach_validation/

This directory contains scripts utilised as part of validation for this approach to penetrance calculation. 

Please refer to the README documentation within the `approach_validation/` directory which describes the contents and validation steps taken.

#### case_studies/

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __case_data.csv__
*  This contains the raw sample data for all case studies described in the publication associated with this repository ([1](https://doi.org/10.1101/2021.03.16.21253691)).


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __published_case_studies.R__
*  This script calls the `adpenetrance` and `getResidualRisk` functions and `case_data.csv` to estimate penetrance for case studies represented within `case_data.csv`.

#### subfunctions/
This directory contains several functions which are utilised internally when running `adpenetrance` or `simADPenetrance`.

Details regarding each function are providied within, and in the repository [wiki](https://github.com/ThomasPSpargo/adpenetrance/wiki).

___  

#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691

2. Al-Chalabi, A. & Lewis, C. M. (2011). Modelling the Effects of Penetrance and Family Size on Rates of Sporadic and Familial Disease. *Human Heredity, 71*(4): 281-288. doi: 10.1159/000330167