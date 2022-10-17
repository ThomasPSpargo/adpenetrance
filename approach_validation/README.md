## ADPenetrance approach validation steps
___Updated 17/10/2022___

_The repository is maintained by Thomas Spargo (<thomas.spargo@kcl.ac.uk>) - please reach out with any questions._

___

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __adpenetrance_MLE_function.R__
* This script generates the `adpenetrance.MLE` function, which adapts the primary `adpenetrance` function to produce additional unadjusted penetrance estimates via a maximum likelihood approach. These can be compared to the unadjusted penetrance estimates generated using the lookup table approach otherwise employed.

* The `adpenetrance.MLE` function is to be called within the `MLE_validation.R` script.
    

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __MLE_validation.R__

* This script calls `adpenetrance_MLE_function.R` and `case_data.csv` to make penetrance estimates for those case studies presented in Table 2 of the associated manuscript ([1](https://doi.org/10.1101/2021.03.16.21253691)).

* This re-analysis compares unadjusted penetrance estimates obtained from the lookup table approach employed within the main `adpenetrance` function to those obtained via a maximum likelihood method.

* Negligable differences are observed between the results of the two approaches.

* Running this script automatically retrieves the `adpenetrance.MLE`,`adpenetrance`, and `getResidualRisk` functions from GitHub.


##### simulation_studies/

* The following scripts each test the accuracy of adjusted penetrance estimates produced by `adpenetrance` in datasets of  simulated families harbouring hypothetical variants for which the true penetrance is known.

* Two distinct populations of simulated families are generated in these simulations. The first is based on the sibship distribution observed in a UK population 1974 birth cohort ([2](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable)). The second is based on the distribtion the 'Next Steps' cohort study of English families ([3](https://doi.org/10.1017/ehs.2020.54)).

* The accuracy of penetrance estimates made for hypothetical variants occuring in these two populations is compared according to the states modelled (see `states` argument [details]( https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance)) and whether the user supplies `adpenetrance` with information about the distribution of sibships in their dataset (see `define_sibstructure` argument [details]( https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance)).

* Full documentation for the adpenetrance function is provided [here]( https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __simulation_1_correct_parameter_specification.R__

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 20 ground truth penetrance values when input parameters for the function (`N` and `RX`) are correctly specified.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __simulation_2_N_parameter_incorrect.R__

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 5 ground truth penetrance values when the mean sibship size parameter (`N`) is incorrectly specified, testing different degrees of misspecification in `N`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __simulation_3_obsRX_parameter_incorrect.R__

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 5 ground truth penetrance values when the observed rate of state X parameter (`RX`) is incorrectly specified, testing different degrees of misspecification in `RX`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __Simulation_4_accountForG.R__

* This script calls `adpenetrance_function.R` to test the adjusted penetrance accuracy according to the degree of residual disease risk *g* among people not harbouring the variant and whether *g* is taken into account when running the approach or is *g* assumed to be 0.

___

#### References
1. Spargo, T. P., Opie-Martin, S., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2021). Calculating variant penetrance using family history of disease and population data. *medRxiv* 2021.03.16.21253691; doi: 10.1101/2021.03.16.21253691

2. Office for National Statistics. Childbearing for women born in different years. 2020. Available from: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable

3. Sheppard P, Monden C. When does family size matter? Sibship size, socioeconomic status and education in England. Evol Hum Sci. 2020; 2, e51:1-21. doi: 10.1017/ehs.2020.54