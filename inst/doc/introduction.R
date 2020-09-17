## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fairsubset)

## -----------------------------------------------------------------------------
#Create normally distributed data
input_list <- list(a= stats::rnorm(100, mean = 3, sd = 2),
b = stats::rnorm(50, mean = 5, sd = 5),
c= stats::rnorm(75, mean = 2, sd = 0.5))

#Run fairsubset, using "mean" as the averaging choice
output <- fairSubset(input_list, subset_setting = "mean", manual_N = 10, random_subsets = 1000)

#Print the report, which describes best and worst subset features
output$report


## -----------------------------------------------------------------------------
#Create bimodal distributed data
input_list <- list(a= c( stats::rnorm(100, mean = 3, sd = 2)
                        ,stats::rnorm(100, mean = 5, sd = 4) ),
b = c( stats::rnorm(50, mean = 5, sd = 5)
      ,stats::rnorm(50, mean = 2, sd = 1)),
#Create skewed data
c= stats::rnbinom(75, 10, .5))

#Run fairsubset, using "ks" to enable ks.test metrics to decide the best subset
output <- fairSubset(input_list, subset_setting = "ks", manual_N = 10, random_subsets = 1000)

#Print the report, which describes best and worst subset features
output$report


## -----------------------------------------------------------------------------
#Run fairsubset, using manual_N = NULL to default to largest possible subsets consistent between all samples
output <- fairSubset(input_list, subset_setting = "ks", manual_N = NULL, random_subsets = 1000)

#Print the size of the largest possible consistent N subset
nrow(output$best_subset)

#Evidently, sample c the lowest number of observed events, so fairsubset chose N = 75 for all subsets. Note that the "best_subset" using this setting retains the original set of data for the smallest sample (ie, output$best_subset$c is the same as input_list$c, just shuffled)

