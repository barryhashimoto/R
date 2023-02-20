# A function to view matched observations of a "matchit" class object after
# coarsened exact matching (CEM). Takes a "matchit" object after coarsened exact
# matching, or a list of such objects (e.g. after multiple imputation and
# matching). Returns either a vector of strings listing all matched IDs or a
# random sample of them. Removes duplicate matches across multiply imputed data
# sets that have each been separately processed by the matching algorithm.
# Function does not check for correct usage.
#
# Barry Hashimoto | barryhashimoto@gmail.com
#
# Feb 20, 2023

library(dplyr)
library(randomNames)
library(Amelia)
library(MatchIt)
library(cem)

# Define function
seeMatchesCEM <- function(match_data, 
                          source_data, 
                          ID, 
                          all = TRUE, 
                          n_sample = NULL) {
  
  # match_data: a data frame or list of data frames (i.e. imputed and matched
  # data frames) that must have been produced by matchit() using 'cem' method.
  #
  # source_data: the source data frame processed by matching and possibly also
  # by imputation. User passes this in so the ID column may be merged.
  #
  # ID: the string name of the vector in `source_data` identifying the
  # observations that will be displayed by the function's return value, the list
  # of matched units.
  #
  # all: returns every set of matches if TRUE. Defaults to TRUE.
  #
  # n_sample: number of matches to sample for returning, if all = FALSE.
  require(MatchIt)
  require(Amelia)
  require(cem)

  # Where user passes in a single data set, nest it as the sole element of a
  # list and prune the source_data file so that listwise deletion does not
  # result in differing row lengths.
  if (class(match_data) == "matchit") { 
    match_data <- list(match_data)
    source_data <- source_data[complete.cases(source_data), ]
    }
  
  # Initalize list to store results
  matches <- list()
  
  # Get matched IDs for each data frame within `match_data`
  for (i in seq_along(match_data)) {
    
    # Get data frame and matching subclass
    dfs <- data.frame(match_data[[i]]$X, 
                      subclass = match_data[[i]][["subclass"]])
    
    # Merge IDs from source to matched
    dfs[[ID]] <- source_data[[ID]]
    
    # Built nested list wherein each element holds the matched names as strings
    matches[[i]] <-
      aggregate(
        dfs[[ID]],
        by = list(dfs[["subclass"]]),
        FUN = function(x) {
          paste0(x, sep = "")
        }
      )[, 2] |>
      lapply(sort)
  }

  # Extract matched names into a sorted vector of strings
  matches <- unlist(matches, recursive = FALSE) |>
    unique() |>
    lapply(function(z) paste0(z, collapse = " & ")) |>
    unlist() |>
    sort()
  
  # Optionally sample from the matches to reduce console output
  if (all == FALSE) {
    matches <- sample(matches, n_sample)
  }
  
  return(matches)
}

# Example: 

# Impute, match, and see matches in the "cem" package's LeLonde data.

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

# Initialize a list to hold matchit objects from multiply imputed data frames
matched_list <- list()

# Declare a linear matching formula
theta <-
  as.formula(paste("treated ~ black + hispanic + married + nodegree + age", 
                   "+ education + u74 + u75 + q1"))

# Perform coarsened exact matching in each imputed data set
for (i in seq_along(imputed_list)) {
  matched_list[[i]] <- 
    matchit(
      formula = theta,
      data = imputed_list[[i]],
      method = "cem"
    )
  }

# Reveal returns of the resulting object
matched_list[[1]] |> str()

# Make a single matchit object using listwise deletion for missing values
matched_data <-
  matchit(
    formula = theta,
    data = LeLonde[complete.cases(LeLonde),],
    method = "cem"
    )

# Sample of matches using unimputed data
seeMatchesCEM(match_data = matched_data,
              source_data = LeLonde,
              ID = "name",
              all = FALSE,
              n_sample = 10) |> 
  print()

# Output:

# [1] "Adrianna Compos & Annnees el-Mustafa & Eden Lockett & Jacob Sobon &
# Michael Heideman"
# [2] "Avery Parks & Brandon Liddell & Dakota Maez & Elizabeth Abdella & Heather
# Parsons"
# [3] "Flint Tutag & Monica Burd"                                                           
# [4] "Aaron Thompson & Lacey Gorrono & Nikola Kim"                                         
# [5] "Melissa Pham & Robert Morquecho & Vy Deck & Ya'qoob al-Khawaja"                      
# [6] "Abraham Dejesus Morales & Rayyaan al-Baten & Samantha Washington"                    
# [7] "Haleema al-Salloum & Jayden Snowhawk"                                                
# [8] "Michaela Cooper & Shafaaa al-Shaban"                                                 
# [9] "Alvin Patterson & Cody Anderson"                                                     
# [10] "Breanna Reynolds & Levi Gutierrez-Solis"

# Sample of matches using imputed data
seeMatchesCEM(match_data = matched_list,
                source_data = LeLonde,
                ID = "name",
                all = FALSE,
                n_sample = 10) |> 
  print()

# Output:
# [1] "Chelsea Ducklow & Nazeeha al-Bacho & Richard Ortega & Shahla al-Fares &
# Shannon Larson"
# [2] "Avery Parks & Brandon Liddell & Dakota Maez & Destiny Geshick & Elizabeth
# Abdella & Heather Parsons"
# [3] "Anthony Debus & Auston Colaizzi & Zainab el-Wahba"                                                                                    
# [4] "Alexander Martinez & Christopher Butcher & Domonique Soto-Rodriguez"                                                                  
# [5] "Aayid al-Mahmood & Brianna Celestino & Marlenny Vogel & Nico Burnsed"                                                                 
# [6] "Alysia Price & Natsima Heideman"                                                                                                      
# [7] "Brandon Vega & Cody Reed & Craig Finney & Mattewos Johnson & Nicholas
# Huntley & Raakaan el-Jabara & Saaliha el-Sharaf & Xavier Powell"
# [8] "Daniel Zerna & Raymond Padilla & Selina Gonzalez-Carrazco"                                                                            
# [9] "Abraham Dejesus Morales & Rayyaan al-Baten & Samantha Washington"                                                                     
# [10] "Jocelynn Taylor & Sabeeka al-Abdulla"        

# All matches using imputed data
if (FALSE){
  seeMatchesCEM(match_data = matched_list,
                source_data = LeLonde,
                ID = "name",
                all = TRUE) |> 
    print()
}

# References: 
# Daniel Ho, Kosuke Imai, Gary King, and Elizabeth Stuart (2007).
# Matching as Nonparametric Preprocessing for Reducing Model Dependence in
# Parametric Causal Inference. Political Analysis 15(3): 199-236.
# http://gking.harvard.edu/files/abs/matchp-abs.shtml. Honaker, J., King, G.,
# Blackwell, M. (2011). Amelia II: A Program for Missing Data. Journal of
# Statistical Software, 45(7), 1–47. URL http://www.jstatsoft.org/v45/i07/.
# Stefano Iacus, Gary King, Giuseppe Porro, “Matching for Casual Inference
# Without Balance Checking: Coarsened Exact Matching,”
# http://gking.harvard.edu/files/abs/cem-abs.shtml
