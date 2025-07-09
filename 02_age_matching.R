#############################################################################################
### Summary:
### 1. 1:1 Matching using MatchIt
#############################################################################################


library(dplyr)
library(survival)
library(MatchIt)


### Case/Control matching --------------------------------------------------------------------------------------------------------
## Unusually there are more cases than controls.
## Matched pairs should have the maximum distance of 0.2 SD derived from age


## On the cluster MatchIt is only available on certain R versions:
# source /share/apps/source_files/R/R-4.3.2.source
# R


matchit_result <- matchit(
  subject_type_logical ~ screener_age,
  data = phenotype,
  method = "nearest",
  caliper = 0.2,      # 0.25 SD of age
  std.caliper = TRUE,  # SD units
  ratio = 1
)
matched_data <- match.data(matchit_result)
