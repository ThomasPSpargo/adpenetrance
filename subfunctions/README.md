## ADPenetrance subfunctions
___Updated 17/10/2022___

_The repository is maintained by Thomas Spargo (<thomas.spargo@kcl.ac.uk>) - please reach out with any questions._
___

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __adpenetrance_errorfit_function.R__
* This script loads the `adpenetrance.errorfit` function which fits a polynomial regression model to predict error in unadjusted penetrance estimates made within `adpenetrance`.
* `adpenetrance.errorfit` will call the `adpenetrance.unadjusted` function in order to make unadjusted penetrance estimates for simulated variants of penetrance between 0 and 1 occuring within a simulated population. The polynomial regression model produced by `adpenetrance.errorfit` is fitted according to differences between the estimated and true penetrance values of each variant simulated.
* The unadjusted penetrance value obtained by `adpenetrance` is then adjusted by error predicted in this estimate under the fitted model to determine the final adjusted penetrance estimate.
* See [1](https://doi.org/10.1101/2021.03.16.21253691) and the ADPenetrance [documentation](https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance) for further details.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __adpenetrance_unadjusted_function.R__
* This script loads the `adpenetrance.unadjusted` function which performs steps 1-3 of the `adpenetrance` approach, providing unadjusted penetrance estimates.
* The function is to be called by both of the `adpenetrance` and `adpenetrance.errorfit` functions.
* Please see the the main ADPenetrance [documentation](https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __varChars.R__
* This script loads the `varChars` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/varChars) which defines a matrix indicating disease onset characteristics of people with and without a variant.
* The function is to be called within the `simADPenetrance` function; and the output of the function is passed to `affAtAge` via `genFamily`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __affAtAge.R__
* This script loads the `affAtAge` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/affAtAge) which pseudo-randomly assigns status of a person being affected (1) or not (0) by a disease based on age, variant status, and disease characteristics associated with variant status in a var.Chars matrix.
* The function is called recurrently by `genFamily` (within `simADPenetrance`) to generate assignments of disease affectedness for individual family members over time.
* It requires an output matrix from `varChars`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __genFamily.R__
* This script loads the `genFamily` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/genFamily) which pseudo-randomly generates a vector detailing a family of disease states "Unaffected"/"Sporadic"/"Familial" across time (or at lifetime risk) according to various parameters and disease characteristics associated with the output of `varChars`.
* The function is called recurrently within `simADPenetrance` to generate populations of simulated families.
* It requires `affAtAge` and an output matrix from `varChars`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __bindRows.R__
* This script loads the `bindRows` [function](https://github.com/ThomasPSpargo/adpenetrance/wiki/subfunctions/bindRows), which produces a matrix by combining a list of vectors or matrices with varying length/ncol, The final value/column of shorter list elements will be duplicated until max length/ncol is reached for all elements and rows can be bound.
* The function is called within `simADPenetrance` to combine output vectors from `genFamily`.

___
  
#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691