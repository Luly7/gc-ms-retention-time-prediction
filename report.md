# GC-MS Retention Time Prediction Validation Report

## Executive Summary

This validation study evaluated retention time prediction performance on a HP-5MS column using 5 aromatic hydrocarbon compounds. The model consistently overestimated retention times with a mean relative error of 22.79%.

## Mean Error

The mean relative error across all compounds was 22.79% with a standard deviation of 3.69%. The absolute error averaged 3.37 minutes, suggesting the model systematically overestimates retention times for aromatic hydrocarbons.

## Agreement Distribution

The predictions showed the following distribution:
- Moderate agreement (10-20% error): 40.0% (2 compounds)
- Poor agreement (>20% error): 60.0% (3 compounds)

No predictions achieved excellent (<5%) or good (5-10%) agreement levels, indicating significant room for improvement in the model's performance for aromatic hydrocarbons.

## Top 3 Best Predictions

1. **Benzene**
   - Experimental RT: 13.27 min
   - Predicted RT: 15.59 min
   - Relative Error: 17.48%
   - Agreement: Moderate

2. **o-Xylene**
   - Experimental RT: 15.81 min
   - Predicted RT: 18.85 min
   - Relative Error: 19.23%
   - Agreement: Moderate

3. **p-Xylene**
   - Experimental RT: 15.27 min
   - Predicted RT: 19.14 min
   - Relative Error: 25.34%
   - Agreement: Poor

The performance analysis shows consistent overestimation across all tested compounds. The systematic bias suggests that model recalibration specifically for aromatic hydrocarbons could significantly improve accuracy. Further validation with a broader range of compound classes is recommended to determine whether this bias extends to other molecular structures.