

fun_pos <- function(x){return(max(x, 0))}

is_increasing <- function(vec) {
  all(diff(vec) >= 0)
}



get_cost <- function(S, T_=10000, c=0, p=8, theta=3, h=1, m=3, gamma=1, type='poi', pois=5, case='FIFO'){
  x_t = rep(0, m)
  C_all = 0
  for(t in 1:T_){
    x_t[m] = S
    if(type == 'poi')
      D_t = rpois(n = 1, lambda = pois)
    else
      D_t = rgeom(1, 1/pois)
    if(case=="FIFO"){
      o_t = gamma * fun_pos(x_t[1]-D_t) + (1-gamma) * fun_pos(x_t[m]-D_t)
      C_t = ((1-gamma) * theta+h) * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + gamma * theta * fun_pos(x_t[1]-D_t) + c * (S-x_t[m-1])
    }
    if(case=="LIFO"){
      o_t = gamma * min(fun_pos(x_t[m]-D_t), x_t[1]) + (1-gamma) * fun_pos(x_t[m]-D_t)
      C_t = ((1-gamma) * theta+h) * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + gamma * theta * min(fun_pos(x_t[m]-D_t), x_t[1]) + c * (S-x_t[m-1])
    }
    
    for(i in 1:(m-1)){
      if(i == 1){
        x_1 = x_t[i]
      }
      if(case=="FIFO")
        x_t[i] = gamma * fun_pos(x_t[i+1] - max(D_t, x_1))
      if(case=="LIFO")
        x_t[i] = gamma * (min(x_t[i+1], fun_pos(x_t[m]-D_t)) - gamma * min(fun_pos(x_t[m]-D_t), x_1))
    }
    C_all = C_all + C_t
  }
  # print(paste('S =',S, 'c =', c, 'p =', p, 'm =', m, 
  #             'gamma =', gamma))
  # print(paste0('C(', S, ') = ', C_all/T_))
  # flush.console()
  return(C_all/T_)
}



get_cost_backlog <- function(S, T_=10000, c=0, p=8, theta=3, h=1, m=1, L=1, gamma=1, type='poi', pois=5){
  x_t = rep(0, m+L)
  C_all = 0
  for(t in 1:T_){
    x_t[m+L] = S
    if(type == 'poi')
      D_t = rpois(n = 1, lambda = pois)
    else
      D_t = rgeom(1, 1/pois)
    o_t = gamma * fun_pos(x_t[1]-D_t) + (1-gamma) * fun_pos(x_t[m]-D_t)
    C_t = ((1-gamma) * theta+h) * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + gamma * theta * fun_pos(x_t[1]-D_t) + c * (S-x_t[m+L-1])
    for(i in 1:(m+L-1)){
      if(i == 1){
        x_1 = x_t[i]
      }
      x_t[i] = gamma * (x_t[i+1] - max(D_t, x_1))
    }
    C_all = C_all + C_t
  }
  # print(paste('S =',S, 'c =', c, 'p =', p, 'm =', m, 
  #             'gamma =', gamma))
  # print(paste0('C(', S, ') = ', C_all/T_))
  # flush.console()
  return(C_all/T_)
}



tune_range <- function(p, gamma, h, theta, PP, c=0, m=3, L=0, type='poi', case="Lost_sales") {
  if(case == "Lost_sales")
    p = p - c
  else
    p = p
  theta = theta + c
  th = (1 - gamma) * theta + h
  if(type == "poi"){
    lower_bound <- qpois(p / (p + theta + h), lower.tail = TRUE, lambda = PP)
    if(case == 'Lost_sales')
      upper_bound <- qpois(p / (p + th), lower.tail = TRUE, lambda = PP)
    else
      upper_bound <- qpois(p / ((gamma*theta - th * sum(gamma^{1:L}) )/(sum(gamma^{0:(m+L-1)})) + p + th), lower.tail = TRUE, lambda = PP)
  }
  else{
    lower_bound <- qgeom(p / (p + theta + h ), lower.tail = TRUE, prob = 1/PP)
    if(case == 'Lost_sales')
      upper_bound <- qgeom(p / (p + th), lower.tail = TRUE, prob = 1/PP)
    else
      upper_bound <- qgeom(p / ((gamma*theta - th * sum(gamma^{1:L}) )/(sum(gamma^{0:(m+L-1)})) + p + th), lower.tail = TRUE, prob = 1/PP)
  }
  return(c(lower_bound, upper_bound))
}



# We will compute the lower bound using monte method
compute_cost <- function(S, c=0, p=8, theta=3, h=1, m=3, gamma=1, pois=5, n_sims=1000, type='poi', case='FIFO') {
  if (type == 'poi') {
    # Computing E[S-D]^+
    expected_S_minus_D <- sum(sapply(0:(20*pois), function(d) {
      max(S - d, 0) * dpois(d, pois)
    }))
    
    # Computing E[D-S]^+
    expected_D_minus_S <- sum(sapply(0:(20*pois), function(d) {
      max(d - S, 0) * dpois(d, pois)
    }))
    
  } else if (type == 'geom') {
    
    # Computing E[S-D]^+
    expected_S_minus_D <- sum(sapply(0:(20*pois), function(d) {
      max(S - d, 0) * dgeom(d, prob=1/pois)
    }))
    
    # Computing E[D-S]^+
    expected_D_minus_S <- sum(sapply(0:(20*pois), function(d) {
      max(d - S, 0) * dgeom(d, prob=1/pois)
    }))
    
  } else {
    stop("Invalid type. Please use 'poi' or 'geom'.")
  }
  
  # Using Monte Carlo simulation for the third term
  simulated_term3 <- replicate(n_sims, {
    if (type == 'poi') {
      D_vector <- rpois(m, pois)
    } else {
      D_vector <- rgeom(m, prob=1/pois)
    }
    
    weights <- gamma^(0:(m-1))
    if(case == 'FIFO'){
      expected_Di <- sum(weights * D_vector)
      max(S*gamma^(m-1) - expected_Di, 0)
    }
    else if(case == 'LIFO') {
      # For LIFO, consider the most recent demand first
      shortages <- gamma^(0:(m-1)) * sapply(rev(1:m), function(j) max(S - D_vector[j], 0))
      min(shortages)
    } else {
      stop("Invalid case. Please use 'FIFO' or 'LIFO'.")
    }

  })
  
  expected_term3 <- mean(simulated_term3)
  expected_D <- ifelse(type=="poi", yes = pois, no=pois-1)
  
  cost <- (h + (theta+c)*(1-gamma))*expected_S_minus_D + (p-c)*expected_D_minus_S + 
          (gamma*(theta+c)/sum(gamma^(0:(m-1))))*expected_term3 + c*expected_D
  
  return(cost)
}





# We will compute the lower bound using monte method
compute_C_tilde_backlog <- function(S, b=8, h=1, theta=3, c=0, m=1, L=1, gamma=1, pois=5, n_sims=10000, type='poi') {
  
  theta = theta + c
  tilde_h = h + (1-gamma) * theta
  if(type == "poi"){
    d = rpois((L + m) * n_sims, lambda = pois)
    D_vectors = matrix(d, nrow = n_sims, ncol = L+m)
  } else {
    d = rgeom((L + m) * n_sims, prob=1/pois)
    D_vectors <-matrix(d, nrow = n_sims, ncol = L+m)
  }
  
  term1_samples <- apply(D_vectors[,1:(L+1)], 1, function(D_vec) {
    fun_pos(gamma^L*S - sum(gamma^(0:L) * D_vec))
  })
  expected_term1 <- mean(term1_samples)
  
  term2_samples <- apply(D_vectors[,1:(L+1)], 1, function(D_vec) {
    fun_pos(sum(gamma^(0:L) * D_vec) - gamma^L*S)
  })
  expected_term2 <- mean(term2_samples)
  
  term3_samples <- apply(D_vectors, 1, function(D_vec) {
    fun_pos(gamma^(m + L - 1)*S - sum(gamma^(0:(m + L - 1)) * D_vec))
  })
  expected_term3 <- mean(term3_samples)
  
  cost <- tilde_h * expected_term1 +
          b * expected_term2 +
          (gamma * (theta - tilde_h * sum(gamma^(0:(L-1))))) / sum(gamma^(0:(L + m - 1))) * expected_term3
  return(cost)
}



test <- function(p_vals, theta_vals, gamma_vals, c_vals, type, PP, h=1, T_=1000,m=3, L=0, case='FIFO', CASE='Lost_sales'){
  # Given parameter values
  
  # Initialize the CSV file with headers
  headers <- data.frame(p = numeric(), theta = numeric(), gamma = numeric(), c = numeric(),
                        optimal_S = numeric(), min_get_cost_value = numeric(),
                        min_compute_cost_value = numeric())
  write.csv(headers,paste0("results_", type, "_", case,"_", CASE,".csv"), row.names = FALSE)

  # Iterate through all combinations of p, theta, and gamma
  for(p in p_vals){
    for(theta in theta_vals){
      for(gamma in gamma_vals){
        for(c in c_vals){
          # Get range of possible S for current parameter combination
          if(CASE=='Lost_sales')
            S_range <- tune_range(p, gamma, h = h, theta, PP = PP, c = c, type=type)
          else
            S_range <- c(5,50)
          
          min_get_cost <- Inf
          min_compute_cost <- Inf
          optimal_S <- NA
          
          # Evaluate costs over the entire S range with a step of 0.1
          for (S in seq(S_range[1], S_range[2], by = 0.1)){
            # Compute cost using get_cost
            if(CASE=='Lost_sales')
              current_get_cost_value <- mean(replicate(T_, get_cost(S, c = c, p = p, 
                                                                m = m, 
                                                                theta = theta, 
                                                                gamma = gamma,
                                                                type = type, pois = PP,
                                                                case = case)))
            else
              current_get_cost_value <- mean(replicate(T_, get_cost_backlog(S, c = c, p = p,
                                                                m = m, L = L,
                                                                theta = theta, 
                                                                gamma = gamma,
                                                                type = type, pois = PP)))
  
            # Check if current S results in a smaller cost
            if(current_get_cost_value < min_get_cost){
              min_get_cost <- current_get_cost_value
              optimal_S <- S
            }
  
            # Compute cost using compute_cost
            if(CASE=='Lost_sales')
              current_compute_cost_value <- compute_cost(S, c = c, p = p, 
                                                      m = m,
                                                      theta = theta, gamma = gamma, 
                                                      pois = PP,
                                                      type = type, case = case)
            else
              current_compute_cost_value <- compute_C_tilde_backlog(S, c = c, b = p, 
                                                      m = m, L = L,
                                                      theta = theta, gamma = gamma, 
                                                      pois = PP,
                                                      type = type)
            
            # Update min_compute_cost if current value is smaller
            if(current_compute_cost_value < min_compute_cost){
              min_compute_cost <- current_compute_cost_value
            }
          }
  
          # Append to results CSV file
          new_row <- data.frame(p = p, theta = theta, gamma = gamma, c = c,
                                optimal_S = optimal_S, 
                                min_get_cost_value = min_get_cost, 
                                min_compute_cost_value = min_compute_cost)
          write.table(new_row, file = paste0("results_", type, "_", case,"_", CASE,".csv"), sep = ",", append = TRUE, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
      }
    }
  }
}

