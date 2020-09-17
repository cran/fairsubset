#' @name fairsubset
#' @description Allows user to obtain subsets of columns of data or vectors within a list. These subsets will match the original data in terms of average and variation, but have a consistent length of data per column. It is intended for use on automated data generation which may not always output the same N per replicate or sample.
#' @param input_list A list, data frame, or matrix.  If matrix or data frame, columns should represent each sample's data.
#' @param subset_setting Choose from c("mean", "median", "ks"). Mean or median will use these averages to choose the best subset.  "ks" will use the Kolmogorov Smirnov test to choose the best subset.   Defaults to "mean".
#' @param manual_N To manually choose how many data points should be in each sample, enter an integer value here.  Otherwise, fairSubset chooses the length of the sample with the most data.  Defaults to NULL.
#' @param random_subsets To manually choose how many random subsets should be used to choose the best subset, enter an integer value here.  Defaults to 1000.
#' @title fairsubset
#' @author Joe Delaney
#' @keywords manip array
#'
#' @return Returns a list.
#' @return $best_subset is a data.frame containing data best representative of original data, given the parameters chosen for fairsubset
#' @return $worst_subset is a data.frame containing data as far from the original as observed in all randomly chosen subsets. It is used solely as a comparator for the worst case scenario from randomly choosing subsets
#' @return $report is a data.frame of averages and variation regarding original data, best subset, and worst subset
#' @return $warning is a character string. If != "", it represents known errors
#'
#' @importFrom matrixStats colMedians
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats ks.test
#' @importFrom stats rnorm
#' @import stats
#' @export
#' @examples
#' input_list <- list(a= stats::rnorm(100, mean = 3, sd = 2),
#' b = stats::rnorm(50, mean = 5, sd = 5),
#' c= stats::rnorm(75, mean = 2, sd = 0.5))
#' fairSubset(input_list, subset_setting = "mean", manual_N = 10, random_subsets = 1000)$report

fairSubset <- function(input_list, subset_setting = "mean", manual_N = NULL, random_subsets = 1000){

    input_data <- NULL
    list_names <- NULL
    warning <- ""

    if( is.list(input_list) == TRUE ){
    if(!is.null(names(input_list))){list_names <- names(input_list)}
    max_N <- max(sapply(1:length(input_list), function(list_item){length(input_list[[list_item]])}),na.rm=TRUE)
    input_list <- data.frame(matrix(unlist(
                              lapply(1:length(input_list), function(list_item){c(input_list[[list_item]], rep(NA, max_N - length(input_list[[list_item]])))})
                ), nrow = max_N, byrow = FALSE),stringsAsFactors=FALSE)
    if(!is.null(list_names)){names(input_list) <- list_names}
    input_data <- input_list
    } else {
    input_data <- as.data.frame(input_list, stringsAsFactors = FALSE)
    }

    pasted_data_header <- colnames(input_data)
    input_data <- lapply(1:ncol(input_data), function(column){unlist(input_data[,column])[!is.na(unlist(input_data[,column]))]})

    if(length(input_data[sapply(unlist(input_data), is.numeric)]) != length(unlist(input_data))){
      suppressWarnings(
        input_data <- lapply(1:length(input_data), function(list_item){
          as.numeric(unlist(input_data[[list_item]]))[!is.na(as.numeric(unlist(input_data[[list_item]])))]
        })
      )
      warning <- "Your data may contain non-numeric values.  These have been removed prior to calculations."
    }

    shortest_data <- min(sapply(1:length(input_data), function(list_item){length(input_data[[list_item]])}))
    if(!is.null(manual_N)){shortest_data <- as.integer(manual_N)}
    sampled_data_report <- as.data.frame(matrix(0, nrow=2, ncol = length(input_data)))
    colnames(sampled_data_report) <- pasted_data_header
    row.names(sampled_data_report) <- c("average","standard deviation")

    sampled_data_best_indices <- as.data.frame(matrix(0, nrow=shortest_data, ncol = length(input_data)))
    colnames(sampled_data_best_indices) <- pasted_data_header

    sampled_data <- as.data.frame(matrix(0, nrow=shortest_data, ncol = length(input_data)))
    colnames(sampled_data) <- pasted_data_header

    data_vector <- rep(0,shortest_data)
    average_value <- 0
    standard_deviation_value <- 0

    if(subset_setting == "mean"){
      average_values_all <- sapply(1:length(input_data), function(list_item){mean(input_data[[list_item]])})
    }

    if(subset_setting %in% c("median", "ks") ){
      average_values_all <- sapply(1:length(input_data), function(list_item){stats::median(input_data[[list_item]])})
    }

    standard_deviation_values_all <- sapply(1:length(input_data), function(list_item){stats::sd(input_data[[list_item]])})

    all_sampled_data <- lapply(1:random_subsets, function(iteration){
      for(column in 1:length(input_data)){
        sampled_data[,column] <- sample(input_data[[column]],size = shortest_data, replace = FALSE)
      }
      return(sampled_data)
    })

    if(subset_setting == "mean"){
      average_values_randomized <- lapply(1:random_subsets, function(iteration){
        abs(colMeans(all_sampled_data[[iteration]]) - average_values_all) #subtract to get distance vector from original
      })
    }
    if(subset_setting %in% c("median", "ks")){
      average_values_randomized <- lapply(1:random_subsets, function(iteration){
        abs(matrixStats::colMedians(as.matrix(all_sampled_data[[iteration]])) - average_values_all)
      })
    }

    standard_deviation_values_randomized <- lapply(1:random_subsets, function(iteration){
      abs(sapply(1:ncol(sampled_data), function(column){stats::sd(all_sampled_data[[iteration]][,column])}) - standard_deviation_values_all) #subtract to get distance vector from original
    })

    average_vector <- rep(0,random_subsets)
    standard_deviation_vector <- rep(0,random_subsets)
    sum_vector <- rep(0,random_subsets)

    if(subset_setting == "ks"){
      KS_pvals <- lapply(1:random_subsets, function(iteration){
        sapply(1:ncol(sampled_data), function(column){
          suppressWarnings({ #otherwise, R fills with tie warnings
            stats::ks.test(input_data[[column]],all_sampled_data[[iteration]][,column], alternative = "two.sided", exact = NULL)$p.value
          })
        })
      })

      best_and_worst_simulations <- lapply(1:length(input_data), function(list_item){

        p_val_vector <- sapply(1:random_subsets, function(iteration){return(KS_pvals[[iteration]][list_item])})
        max_pval <- max(p_val_vector)
        min_pval <- min(p_val_vector)

        best_simulations_for_further_testing  <- as.integer(which(p_val_vector == max(p_val_vector)))
        worst_simulations_for_further_testing <- as.integer(which(p_val_vector == min(p_val_vector)))

        if(length(best_simulations_for_further_testing) > 1 | length(worst_simulations_for_further_testing) > 1){

          average_vector <-
            sapply(1:random_subsets, function(iteration){
              average_values_randomized[[iteration]][list_item]
            }) /
            sum(sapply(1:random_subsets, function(iteration){ #dividing by sum results in equal weight of average and standard deviation
              average_values_randomized[[iteration]][list_item]
            })) * random_subsets

          standard_deviation_vector <-
            sapply(1:random_subsets, function(iteration){
              standard_deviation_values_randomized[[iteration]][list_item]
            }) /
            sum(sapply(1:random_subsets, function(iteration){ #dividing by sum results in equal weight of average and standard deviation
              standard_deviation_values_randomized[[iteration]][list_item]
            })) *random_subsets

          sum_vector <- sapply(1:random_subsets, function(iteration){sum( c(average_vector[iteration],standard_deviation_vector[iteration]), na.rm=TRUE)}) #determine best simulation for given column

          sum_vector_best  <- sum_vector[best_simulations_for_further_testing]
          sum_vector_worst <- sum_vector[worst_simulations_for_further_testing]

          return(list(
            best  = min(as.integer(which(sum_vector  == min(sum_vector_best))),na.rm=TRUE )
            ,worst = min(as.integer(which(sum_vector == min(sum_vector_worst))),na.rm=TRUE )
          ))

        } else {
          return(list(
            best  = as.integer(which(p_val_vector == max(p_val_vector)))
            ,worst = as.integer(which(p_val_vector == min(p_val_vector)))
          ))
        }

      })
      best_simulations <- as.numeric(unlist(sapply(1:length(input_data), function(list_item){best_and_worst_simulations[[list_item]]$best})))
      worst_simulations <-   as.numeric(unlist(sapply(1:length(input_data), function(list_item){best_and_worst_simulations[[list_item]]$worst})))
    }

    best_simulations <- unlist(sapply(1:length(input_data), function(list_item){

      if(sum(sapply(1:random_subsets, function(iteration){
        average_values_randomized[[iteration]][list_item]
      })) == 0
      ){return(1)} else {

        average_vector <-
          sapply(1:random_subsets, function(iteration){
            average_values_randomized[[iteration]][list_item]
          }) /
          sum(sapply(1:random_subsets, function(iteration){ #dividing by sum results in equal weight of average and standard deviation
            average_values_randomized[[iteration]][list_item]
          })) * random_subsets

        standard_deviation_vector <-
          sapply(1:random_subsets, function(iteration){
            standard_deviation_values_randomized[[iteration]][list_item]
          }) /
          sum(sapply(1:random_subsets, function(iteration){ #dividing by sum results in equal weight of average and standard deviation
            standard_deviation_values_randomized[[iteration]][list_item]
          })) *random_subsets

        sum_vector <- average_vector + standard_deviation_vector #determine best simulation for given column

        return(as.integer(which(sum_vector == min(sum_vector))))
      }

    }))

    worst_simulations <- unlist(sapply(1:length(input_data), function(list_item){

      if(sum(sapply(1:random_subsets, function(iteration){
        average_values_randomized[[iteration]][list_item]
      })) == 0
      ){return(1)} else {

        average_vector <-
          sapply(1:random_subsets, function(iteration){
            average_values_randomized[[iteration]][list_item]
          }) /
          sum(sapply(1:random_subsets, function(iteration){ #dividing by sum results in equal weight of average and standard deviation
            average_values_randomized[[iteration]][list_item]
          })) *random_subsets

        standard_deviation_vector <-
          sapply(1:random_subsets, function(iteration){
            standard_deviation_values_randomized[[iteration]][list_item]
          }) /
          sum(sapply(1:random_subsets, function(iteration){ #dividing by sum results in equal weight of average and standard deviation
            standard_deviation_values_randomized[[iteration]][list_item]
          })) *random_subsets

        sum_vector <- average_vector + standard_deviation_vector #determine best simulation for given column

        return(as.integer(which(sum_vector == max(sum_vector))))
      }

    }))

    for(column in 1:ncol(sampled_data)){
      sampled_data[,column] <- all_sampled_data[[best_simulations[column]]][,column]
    }

    worst_sampled_data <- sampled_data
    for(column in 1:ncol(worst_sampled_data)){
      worst_sampled_data[,column] <- all_sampled_data[[worst_simulations[column]]][,column]
    }

    report <- as.data.frame(matrix(0,nrow = 6, ncol=length(input_data)))
    colnames(report) <- pasted_data_header
    if(subset_setting == "mean"){
      row.names(report) <- c("Mean of original data", "Mean of best subset of data", "Mean of worst subset of data",
                             "Standard deviation of original data", "Standard deviation of best subset of data", "Standard deviation of worst subset of data")
      report["Mean of original data",] <- average_values_all
      report["Mean of best subset of data",] <- colMeans(sampled_data)
      report["Mean of worst subset of data",] <- colMeans(worst_sampled_data)
      report["Standard deviation of original data",] <- standard_deviation_values_all
      report["Standard deviation of best subset of data",] <- sapply(1:length(input_data), function(column){stats::sd(sampled_data[,column])})
      report["Standard deviation of worst subset of data",] <- sapply(1:length(input_data), function(column){stats::sd(worst_sampled_data[,column])})
    } else { #Median or KS
      row.names(report) <- c("Median of original data", "Median of best subset of data", "Median of worst subset of data",
                             "Standard deviation of original data", "Standard deviation of best subset of data", "Standard deviation of worst subset of data")
      report["Median of original data",] <- average_values_all
      report["Median of best subset of data",] <- matrixStats::colMedians(as.matrix(sampled_data))
      report["Median of worst subset of data",] <- matrixStats::colMedians(as.matrix(worst_sampled_data))
      report["Standard deviation of original data",] <- standard_deviation_values_all
      report["Standard deviation of best subset of data",] <- sapply(1:length(input_data), function(column){stats::sd(sampled_data[,column])})
      report["Standard deviation of worst subset of data",] <- sapply(1:length(input_data), function(column){stats::sd(worst_sampled_data[,column])})
    }

    if(warning != ""){warning(warning)}

    return(list(
       best_subset = sampled_data
      ,worst_subset = worst_sampled_data
      ,report = report
      ,warning = as.character(warning)

    ))

  }
