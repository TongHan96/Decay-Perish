## Optimal Policy Loss: A Quick Guide

### Introduction

In order to derive optimal policy loss in supply chain and inventory management, we can utilize computational approaches. In this guide, we'll delve into how this can be done using R to get the optimal policy loss under different scenarios: Lost-sales and Backlogging.

### Prerequisite

Ensure you have R and necessary libraries installed. You may need to install additional R packages and load them in your script, depending on your requirements.

### Code Example

```r
source('new.R')
# OR
source('optimal_n.R')
start_time = Sys.Date()
```

#### Lost-sales

In lost-sales scenarios, we assess situations where missed sales opportunities incur a cost or penalty. Example code snippet is provided below.

```r
suppressWarnings(
    test(
        p_vals = c(8, 20, 40), 
        theta_vals = c(3, 6, 8), 
        gamma_vals = c(1, 0.5, 0.1), 
        c_vals = c(0, 5), 
        type = 'poi', 
        PP = 5, 
        h = 1, 
        T_ = 100
    )
)

suppressWarnings(
    test(
        p_vals = c(8, 20, 40), 
        theta_vals = c(3, 6, 8), 
        gamma_vals = c(1, 0.5, 0.1), 
        c_vals = c(0, 5), 
        type = 'geom', 
        PP = 6, 
        h = 1, 
        T_ = 100, 
        case = 'LIFO'
    )
)

```

#### Backlogging

Backlogging involves scenarios where unfulfilled orders accumulate and are fulfilled at a later time. An example code snippet is provided below.

```r
suppressWarnings(
    test(
        p_vals = c(8, 20, 40), 
        theta_vals = c(3, 6, 8), 
        gamma_vals = c(1, 0.5, 0.1), 
        c_vals = c(0), 
        type = 'poi', 
        PP = 5, 
        h = 1, 
        T_ = 100,
        m = 1,
        L = 1,
        CASE = 'Back'
    )
)

```

### Optimal n
```r
# Assuming compute_cost returns a numeric value
results <- data.frame(m=integer(), n=integer(), type=character(), gamma=numeric(), cost=numeric(), stringsAsFactors = FALSE)

for(m in c(3, 4, 5)) {
  for(n in c(1, 3, 5)) {
    for(type in c('poi', 'geom')) {
      for(gamma in c(0.1, 0.5)){
        cost <- test_n(m=m, n=n, type=type, gamma=gamma)
        result <- data.frame(m=m, n=n, type=type, gamma=gamma, cost=cost)
        results <- rbind(results, result)
      }
    }
  }
}

write.csv(results, "gamma_simulation_results.csv", row.names = FALSE)

```

### Explanation of Parameters

- `p_vals`, `theta_vals`, `gamma_vals`, and `c_vals`: These are different parameter values for testing optimal policy loss under varying conditions.
  
- `type`: Indicates the type of distribution ('poi' for Poisson, 'geom' for Geometric).

- `PP`: Policy parameter.

- `h`: Holding cost per unit per unit of time.

- `T_`: Time.

- `case` (or `CASE`): The case being tested, e.g., 'LIFO' for Last-In-First-Out or 'Back' for Backlogging.

### Note

The above code snippets are assumed to be calling a function named `test()` which is not a built-in function in R and thus should be defined elsewhere in your code with the appropriate logic to handle the parameters and scenarios.

### Conclusion

With this guide, you can explore different scenarios under lost-sales and backlogging by modifying the parameter values as per your requirement. For a comprehensive evaluation, consider analyzing the results under various parameter combinations and summarizing the insights in a tabular or graphical format for better understanding and decision-making.
