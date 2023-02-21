# A function to view matched observations in a "matchit" class object after
# genetic matching. Takes a "matchit" object where method = "genetic" or a list
# of them (e.g. after multiple imputation and matching). Returns either a vector
# of strings listing all matched IDs or a random sample of them. Removes
# duplicate matches across imputed data sets. Function does not check for
# correct usage.
#
# Barry Hashimoto | barryhashimoto@gmail.com
#
# Feb 20, 2023

library(dplyr)
library(randomNames)
library(Amelia)
library(MatchIt)
library(rgenoud)

# Define function
seeMatchesGenetic <- function(match_data,
                              ID,
                              all = TRUE,
                              n_sample = NULL) {
  
  # match_data: a data frame or list of data frames (i.e. imputed and matched
  # data frames) that must have been produced by matchit() using 'cem' method.
  #
  # ID: the string name of the vector in `source_data` identifying the
  # observations that will be displayed by the function's return value, the list
  # of matched units.
  #
  # all: returns every set of matches if TRUE. Defaults to TRUE.
  #
  # n_sample: number of matches to sample for returning, if all = FALSE.
  require(rgenoud)
  require(MatchIt)
  require(Amelia)
  require(dplyr)

  # Where user passes in a single data set, nest it as the sole element of a
  # list and prune the source_data file so that listwise deletion does not
  # result in differing row lengths.
  if (class(match_data) == "matchit") { 
    match_data <- list(match_data)
  }
  
  matches <- list()
  
  for(i in seq_along(match_data)){
    
    # Get data and matching model results
    model_data <- match_data[[i]]$model$data
  	match_matrix <- match_data[[i]]$match.matrix
  	
  	# Create a dictionary of matched rownames and ID names
  	nameKey <- 
  	  cbind(
  	    rownames(model_data), 
  	    as.character(model_data[[ID]])
  	    ) |> 
  	  as.data.frame(stringsAsFactors = FALSE)
	
  	# Create a dictionary of matched rownames
  	matchMat <- 
  	  cbind(
  	    rownames(match_matrix), 
  	    match_matrix
  	    ) |> 
  	  as.data.frame(stringsAsFactors = FALSE)
  	
  	# Use the empty-index trick to replace rownames in matchMat with their
  	# corresponding ID names from nameKey. The unlist function is used to convert
  	# the matchMat data frame into a vector for efficient matching. The base
  	# match function is used to find the indices of the matching rownames between
  	# matchMat and nameKey$V1. And nameKey$V2 is used to replace the values in
  	# matchMat with the corresponding ID names.
  	matchMat[] <- nameKey$V2[match(unlist(matchMat), nameKey$V1)]
  	
  	# Format and sort all the matches into a vector of strings. Split the
  	# matchMat data frame into a list of data frames, with one data frame for
  	# each row Remove any columns with missing values from each data frame.
  	# Convert each data frame to a single string by concatenating the non-NA
  	# values separated by "&". Combine the strings into a single vector. Then
  	# sort the vector alphabetically.
  	matches[[i]] <- split(matchMat, f = seq(nrow(matchMat))) |> 
  	  lapply(function(x) x[, !is.na(x)]) |> 
  	  lapply(function(x) paste(x, collapse = " & ")) |>
  	  unlist(use.names = FALSE) |> 
  	  sort()
  	}
  
  # Extract matched names into a sorted vector of strings
  matches <- unlist(matches) |> 
    unique() |> 
    sort()

  # Optionally sample from the matches to reduce console output
  if (all == FALSE) {
    matches <- sample(matches, n_sample)
  }
  
  return(matches)
  }

# Example: impute, match, and see matches in the "cem" package's LeLonde data.
# Set unit IDs to random names.

# Load the LeLonde data and set the unit IDs to random names.
data(LeLonde)

LeLonde <-
  LeLonde |> 
  mutate(
    name = randomNames::randomNames(
      n = nrow(LeLonde),
      name.order = "first.last",
      name.sep = " "
    )
  )

# Make a single matchit object, using listwise deletion for missing values.
listwise_deleted_data <- LeLonde[complete.cases(LeLonde), ]

# Declare a linear matching formula
theta <-
  as.formula(paste("treated ~ black + hispanic + married + nodegree + age", 
                   "+ education + u74 + u75 + q1"))

# Perform genetic matching
matched_data <-
  matchit(
    formula = theta,
    data = listwise_deleted_data,
    method = "genetic"
    )

# Impute two data sets via the multivariate normal model from Amelia II
imputed_list <-
  amelia(
    LeLonde,
    m = 2,
    p2s = F,
    idvars = "name",
    noms = c("black",
             "hispanic",
             "treated",
             "married",
             "nodegree",
             "u74",
             "u75",
             "q1"
    )
  )$imputations

# Print first element of the returned data
imputed_list[[1]] |> glimpse()

# Initialize a list to hold matched and imputed data sets
matched_list <- list()

# Populate the list with matched, imputed data sets
for (i in seq_along(imputed_list)) {
  matched_list[[i]] <-
    matchit(
      formula = theta,
      data = imputed_list[[i]],
      method = "genetic"
    )
}

# Reveal returns of the resulting object
matched_list[[1]] |> str()

# See matches for a single data set after listwise deletion
seeMatchesGenetic(matched_data, "name", all = FALSE, n_sample = 10)

# Output:  
# [1] "Maria Prutch & Joshua Loya"              
# [2] "Elizabeth Harris & Nicholas Branaman"    
# [3] "Robert Morquecho & Caitlin Domingo"      
# [4] "Korena Chandler & Carlo Hyden-Terry"     
# [5] "Ruqayya al-Mahmud & Jay Heredia"         
# [6] "Ian Clark & Irving Silva"                
# [7] "Hui Ouch & Olivia Giron Lovato"          
# [8] "Nickolas Roberts & Adrian Christopher Ha"
# [9] "James Flyinghawk & Brian Alvarez"        
# [10] "Shannon Larson & Chelsea Ducklow"     

# See a sample of matches for multiply imputed data sets, duplicates removed
seeMatchesGenetic(matched_list, "name", all = FALSE, n_sample = 10)

# Output: 
# [1] "Mahdeeya el-Jamil & Al Lyn Solomon"         
# [2] "Kendall Shibre & Tayler Jennings"           
# [3] "Corbin Curtin & Melissa Pham"               
# [4] "Shureetha Mendoza-Lagunas & Samantha Jarmon"
# [5] "Naseer el-Hamed & Breanna Coley"            
# [6] "Kristopher House & Steven Cathey"           
# [7] "Jonee Porter & Alicia Her"                  
# [8] "Auston Colaizzi & Brittney Lupinski"        
# [9] "Imtinaan al-Farrah & Earl Allen"            
# [10] "Hasana al-Rahman & Crystal Minnich"     

# See all matches for the above data set
if (FALSE){
  seeMatchesGenetic(matched_list, "name", all = TRUE)
}

# References
# Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart (2007). Matching as
# Nonparametric Preprocessing for Reducing Model Dependence in Parametric Causal
# Inference. Political Analysis 15(3): 199-236.
# http://gking.harvard.edu/files/abs/matchp-abs.shtml Honaker, J., King, G.,
# Blackwell, M. (2011). Amelia II: A Program for Missing Data. Journal of
# Statistical Software, 45(7), 1–47. URL http://www.jstatsoft.org/v45/i07/.
# Walter Mebane, Jr. and Jasjeet S. Sekhon. 2011. ``Genetic Optimization Using
# Derivatives: The rgenoud package for R.'' Journal of Statistical Software,
# 42(11): 1-26. Jasjeet S. Sekhon. 2011. ``Multivariate and Propensity Score
# Matching Software with Automated Balance Optimization: The Matching package
# for R.'' Journal of Statistical Software, 42(7): 1-52.
