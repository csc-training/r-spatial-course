library(tictoc)

# Just a demo slow function, that waits for 1 second
slow_function<-function(i) {
  Sys.sleep(1) 
  return(i)
}
# Input data vector, the slow function is run for each element.
input = 1:7

# SERIAL options
# Basic FOR loop
a <- 0
tic()
for(i in input) {  
  a[i] <- slow_function(i)
}
toc()

# Basic lapply
tic()
b <- lapply(input, slow_function)
toc()

# purrr, map
library(purrr) 
tic()
c <- map(input, slow_function)
toc()

# PARALLEL options with future
library(future.apply) 
plan(multisession)
#options(future.availableCores.methods = "Slurm")
future::availableCores()

# Parallel lapply
tic()
d <- future_lapply(input, slow_function)
toc()

# Parallel map
library(furrr) 
tic()
e <- future_map(input, slow_function)
toc()
