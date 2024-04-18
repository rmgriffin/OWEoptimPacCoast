rm(list=ls())
PKG <- c("nsga2R","tidyverse","rgenoud","lpSolve","rPref","mco")

for (p in PKG) {
  if(!require(p,character.only = TRUE)) {  
    install.packages(p)
    require(p,character.only = TRUE)}
}

rm(p,PKG)

# mco nsga2 with repair function 
# from this (modified) example https://stackoverflow.com/questions/38009498/global-multi-optimization-function-specification-in-r
set.seed(0)
all.options<-data.frame(num.option=1:100,main.goal1 = abs(rnorm(100)),
                        main.goal2 = abs(rnorm(100)),
                        main.goal3 = abs(rnorm(100)),soil=c(rep("soilType1",50),rep("soilType2",50))) # all possible combinations of the 3 goals

repair_solution <- function(solution, n_select, n_items) { # Turns a draw of D real numbers into a draw of D unique integer numbers with a range of 1:n_items 
  unique_solution <- unique(round(solution),0)
  
  while(length(unique_solution) < n_select) {
    missing_items <- setdiff(1:n_items, unique_solution)
    unique_solution <- c(unique_solution, sample(missing_items, 1))
  }
  
  if(length(unique_solution) > n_select) {
    unique_solution <- sample(unique_solution, n_select)
  }
  
  return(sort(unique_solution))
}

optimized_repair_solution <- function(solution, n_select, n_items) {
  # Ensure the solution does not exceed the desired length, and remove duplicates
  unique_solution <- unique(round(solution),0)
  
  # If the unique solution has more items than needed, randomly select n_select items
  if(length(unique_solution) > n_select) {
    return(sample(unique_solution, n_select))
  }
  
  # If the unique solution has fewer items than needed, fill in with missing items
  if(length(unique_solution) < n_select) {
    # Pre-calculate all missing items once
    all_items <- 1:n_items
    missing_items <- all_items[!all_items %in% unique_solution]
    # Randomly select enough missing items to reach n_select
    to_add <- sample(missing_items, n_select - length(unique_solution))
    unique_solution <- c(unique_solution, to_add)
  }
  
  # Return the repaired solution, ensuring it's sorted for consistency
  return(sort(unique_solution))
}

main.goal1<-function(x){    # x - a vector of row indices
  x_repaired <- optimized_repair_solution(x, D, nrow(all.options))
  return(sum(all.options[x_repaired, 2]))}

main.goal2<-function(x){    # x - a vector of row indices
  x_repaired <- optimized_repair_solution(x, D, nrow(all.options))
  return(mean(all.options[x_repaired, 3]))}

eval<-function(x){
  return(c(main.goal1(x),main.goal2(x)))} # objective function

D<-50 # Number of values in the function must equal D (is replacement an issue, ie. could it use the same observation several times in a set?)

system.time(G<-nsga2(fn=eval,
         idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn 
         odim=2, # Output dimensions
         lower.bounds=rep(0,D),
         upper.bounds=rep(as.numeric(nrow(all.options)),D),
         popsize=200,generations=1000, cprob = 0.7, cdist = 5,
         mprob = 0.2, mdist = 10))

plot(G)

# mco nsga2 with penalty  
# from this (modified) example https://stackoverflow.com/questions/38009498/global-multi-optimization-function-specification-in-r
set.seed(0)
all.options<-data.frame(num.option=1:100,main.goal1 = abs(rnorm(100)),
                        main.goal2 = abs(rnorm(100)),
                        main.goal3 = abs(rnorm(100)),soil=c(rep("soilType1",50),rep("soilType2",50))) # all possible combinations of the 3 goals

psf<-1000 # Penalty scale factor for duplicates

main.goal1<-function(x){    # x - a vector of row indices
  main.goal1=sum(all.options[x,2]) # compute main.goal1 - subsetting like this, even if you get a decimal continuous value, it returns the index of the integer the value starts with
  return(main.goal1)}

main.goal2<-function(x){    # x - a vector 
  main.goal2=mean(all.options[x,3]) # compute main.goal2
  return(main.goal2)}

calculate_penalty <- function(x) {
  length(all.options[x,"num.option"]) - length(unique(all.options[x,"num.option"]))} # Having an index in the dataframe allows the penalty function to address decimal draws in x

eval<-function(x){
  penalty <- calculate_penalty(x) * psf
  
  obj1_val <- main.goal1(x) + penalty
  obj2_val <- main.goal2(x) + penalty
  
  return(c(obj1_val, obj2_val))} # objective function
  
D<-50 # Number of values in the function must equal D (is replacement an issue, ie. could it use the same observation several times in a set?)

G<-nsga2(fn=eval,
        idim=D, # Length of the subset of indices drawn by nsga2 through the eval fn 
        odim=2, # Output dimensions
        lower.bounds=rep(0,D),
        upper.bounds=rep(as.numeric(nrow(all.options)),D),
        popsize=200,generations=1000, cprob = 0.7, cdist = 5,
        mprob = 0.2, mdist = 10)

plot(G)



## rPref only marks pareto from dominated
p <- high(mtcars$mpg) * high(mtcars$hp)
psel(mtcars,p)

## Genoud

# Example dataframe
set.seed(42)
x <- data.frame(
  vector1 = runif(10, 1, 100),  # Example data
  vector2 = runif(10, 1, 100)
)

# Desired fixed sample size
n <- 3

# Modified Objective Function
objective_function <- function(x) {
  selected_indices <- sample(nrow(x), n)
  
  # Calculate objectives
  sum_v1 <- -sum(x$vector1[selected_indices])  # Objective 1: Minimize sum of v1
  avg_v2 <- -mean(x$vector2[selected_indices])  # Objective 2: Minimize average of v2
  
  return(c(sum_v1, avg_v2))
}

objective_function(x)

# Optimization setup (conceptual, as genoud might not directly apply without custom adaptation)
result <- genoud(
  fn = objective_function,
  nvars = 2,  # Number of decision variables, corresponding to dataframe rows
  pop.size = 1000,  # Population size
  max.generations = 100  # Max number of generations
)

# Please note: This setup assumes an adaptation for ranking-based selection,
# which may require additional implementation effort to align with genoud's capabilities.



## NSGA2 - doesn't seem to work because it can't respect paired vectors as an input

# Assuming a dataset or vector of values
set.seed(42)
dataset <- runif(100, min = 0, max = 100)  # A vector of 100 random values

# Define objective functions
# Objective 1: Minimize variance of a sample of size n from the dataset
objective1 <- function(x) {
  sample_indices <- sample(1:length(dataset), size = 5, replace = FALSE)
  sample_variance <- var(dataset[sample_indices])
  return(sample_variance)
}

# Objective 2: Maximize the mean of a different sample of size n from the dataset
objective2 <- function(x) {
  sample_indices <- sample(1:length(dataset), size = 5, replace = FALSE)
  sample_mean <- mean(dataset[sample_indices])
  return(-sample_mean)  # Negate to convert to a minimization problem
}

# Define the optimization problem parameters
num_variables <- 2  # Two decision variables influencing sample size or characteristics
lower_bounds <- c(1, 1)  # Minimum sample size of 1 (assuming indices or parameters start at 1)
upper_bounds <- c(50, 50)  # Arbitrarily chosen upper bound

# Run NSGA-II
result <- nsga2(
  fn = function(x) c(objective1(x), objective2(x)),  # Combined objective functions
  popsize = 100,  # Population size
  generations = 100,  # Number of generations
  idim = 2,
  odim = 2,
  lower.bounds = lower_bounds,
  upper.bounds = upper_bounds
)

plot(result, xlab="y1", ylab="y2", main="Objective space")

# Extract Pareto-optimal solutions
pareto_solutions <- as.data.frame(result) %>% as.data.frame() %>% filter(pareto.optimal == TRUE)

# Plot the Pareto frontier (using transformed objectives to original scale)
plot(-pareto_solutions[,1], -pareto_solutions[,2], 
     xlab = "Sample Variance", ylab = "Sample Mean", 
     main = "Pareto Frontier", col = "blue", pch = 19)
