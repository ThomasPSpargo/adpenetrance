## ADPenetrance approach validation steps
___Updated 16/12/2022___

_The repository is maintained by Thomas Spargo (<thomas.spargo@kcl.ac.uk>) - please reach out with any questions._

___

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __adpenetrance_MLE_function.R__
* This script generates the `adpenetrance.MLE` function, which adapts the primary `adpenetrance` function to produce additional unadjusted penetrance estimates via a maximum likelihood approach. These can be compared to the unadjusted penetrance estimates generated using the lookup table approach otherwise employed.

* The `adpenetrance.MLE` function is to be called within the `MLE_validation.R` script.
    

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __MLE_validation.R__

* This script calls `adpenetrance_MLE_function.R` and `case_data.csv` to make penetrance estimates for those case studies presented in Table 2 of the associated manuscript ([1](https://doi.org/10.1186/s13073-022-01142-7)).

* This re-analysis compares unadjusted penetrance estimates obtained from the lookup table approach employed within the main `adpenetrance` function to those obtained via a maximum likelihood method.

* Negligable differences are observed between the results of the two approaches.

* Running this script automatically retrieves the `adpenetrance.MLE`,`adpenetrance`, and `getResidualRisk` functions from GitHub.


##### simulation_studies/

The following scripts each test the accuracy of adjusted penetrance estimates produced by `adpenetrance` in datasets of  simulated families harbouring hypothetical variants for which the true penetrance is known.

Two distinct populations of simulated families are generated in these simulations. The first is based on the sibship distribution observed in a UK population 1974 birth cohort ([2](https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable)). The second is based on the distribtion the 'Next Steps' cohort study of English families ([3](https://doi.org/10.1017/ehs.2020.54)).

The accuracy of penetrance estimates made for hypothetical variants occuring in these two populations is compared according to the states modelled (see `states` argument [details]( https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance)) and whether the user supplies `adpenetrance` with information about the distribution of sibships in their dataset (see `define_sibstructure` argument [details]( https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance)).

Full documentation for the adpenetrance function is provided [here]( https://github.com/ThomasPSpargo/adpenetrance/wiki/ADPenetrance).

The set of `TimeSimulation*.R` scripts, draw upon the `simADPenetrance` [function]( https://github.com/ThomasPSpargo/adpenetrance/wiki/simADPenetrance) to simulate the impact of age of sampling upon accuracy of lifetime penetrance estimation across several scenarios, described below.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __simulation_1_correct_parameter_specification.R__

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 20 ground truth penetrance values when input parameters for the function (`N` and `RX`) are correctly specified.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __simulation_2_N_parameter_incorrect.R__

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 5 ground truth penetrance values when the mean sibship size parameter (`N`) is incorrectly specified, testing different degrees of misspecification in `N`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __simulation_3_obsRX_parameter_incorrect.R__

* This script calls `adpenetrance_function.R` to test the accuracy of adjusted penetrance estimates across 5 ground truth penetrance values when the observed rate of state X parameter (`RX`) is incorrectly specified, testing different degrees of misspecification in `RX`.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __Simulation_4_accountForG.R__

* This script calls `adpenetrance_function.R` to test the adjusted penetrance accuracy according to the degree of residual disease risk *g* among people not harbouring the variant and whether *g* is taken into account when running the approach or is *g* assumed to be 0.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __TimeSimulation1_onesample.R__

* This script simulates a scenario where the only families harbouring a variant of penetrance 'f' have been sampled.

* Disease state rates are expected to vary over time, affecting accuracy in estimating lifetime penetrance according to age of sampling.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __TimeSimulation2_equalonset.R__

* This script simulates a scenario where families harbouring variant *M* which has penetrance f are sampled alongside other families where the variant doesnt occur, instead harbouring one of several other variants with varying penetrance.

* Equal variability of disease onset is observed between people with variant *M* and people with other variants. Therefore change in family disease state rates are expected to be more stable across ages of sampling. This permits reasonable lifetime penetrance estimates at younger ages.


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __TimeSimulation3_1pt3onset.R__

* This script simulates a scenario where families harbouring variant *M* which has penetrance f are sampled alongside other families where the variant doesnt occur, instead harbouring one of several other variants with varying penetrance.

* Unlike `TimeSimulation2_equalonset.R`, equal onset variability is not observed in this simulation, with a  shorter 'window' for disease onset among people harbouring variant *M*, scaled to be 1.3 times shorter than the onset window of people without *M*. Lifetime penetrance estimation is expected to more influenced by sampling age than in `TimeSimulation2_equalonset.R`, but still tolerable by this degree of unequal onset variability.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; __TimeSimulation4_0pt77onset.R__

* This script simulates a scenario where families harbouring variant *M* which has penetrance f are sampled alongside other families where the variant doesnt occur, instead harbouring one of several other variants with varying penetrance.

* Unlike `TimeSimulation2_equalonset.R`, equal onset variability is not observed in this simulation, with a  longer 'window' for disease onset among people harbouring variant *M*, scaled to be 0.77 shorter (i.e. $\ 1/0.77≈1.3 $ times longer) than the onset window for people without *M*. Lifetime penetrance estimation is expected to more influenced by sampling age than in `TimeSimulation2_equalonset.R`, but still tolerable by this degree of unequal onset variability.

___

#### References
1. Spargo, T.P., Opie-Martin, S., Bowles, H., Lewis, C. M., Iacoangeli, A., & Al-Chalabi, A. (2022). Calculating variant penetrance from family history of disease and average family size in population-scale data. *Genome Med* 14, 141. doi: 10.1186/s13073-022-01142-7

2. Office for National Statistics. Childbearing for women born in different years. 2020. Available from: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/conceptionandfertilityrates/datasets/childbearingforwomenbornindifferentyearsreferencetable

3. Sheppard P, Monden C. When does family size matter? Sibship size, socioeconomic status and education in England. Evol Hum Sci. 2020; 2, e51:1-21. doi: 10.1017/ehs.2020.54