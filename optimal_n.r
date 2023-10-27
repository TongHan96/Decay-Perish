fun_pos <- function(x){return(max(x, 0))}

is_increasing <- function(vec) {
  all(diff(vec) >= 0)
}

tune_range <- function(p, gamma, h, theta, PP, c=0, m=3, L=0, type='poi') {
  p = p - c
  theta = theta + c
  th = (1 - gamma) * theta + h
  if(type == "poi"){
    lower_bound <- qpois(p / (p + theta + h), lower.tail = TRUE, lambda = PP)
    upper_bound <- qpois(p / (p + th), lower.tail = TRUE, lambda = PP)
  }
  else{
    lower_bound <- qgeom(p / (p + theta + h ), lower.tail = TRUE, prob = 1/PP)
    upper_bound <- qgeom(p / (p + th), lower.tail = TRUE, prob = 1/PP)
  }
  return(c(lower_bound, upper_bound))
}


get_cost <- function(S, T_=10000, c=0, p=8, theta=3, h=1, m=3, gamma=1, type='poi', pois=5, case='FIFO', Bino=FALSE, n=0){
  x_t = rep(0, m)
  C_all = 0
  for(t in 1:T_){
    # generate N(n)
    if(n == 0)
      N_n = 1
    else
      N_n =  which.max(sapply(1:10*n, function(i) sum(rexp(10*n,1)[1:i]))<=n)
    x_t[m] = S
    # print(x_t)
    if(type == 'poi')
      D_t = rpois(n = N_n, lambda = pois)
    else
      D_t = rgeom(n = N_n, 1/pois)
    if(case=="FIFO"){
      if(Bino==FALSE){
        o_t = gamma * fun_pos(x_t[1]-D_t) + (1-gamma) * fun_pos(x_t[m]-D_t)
        C_t = ((1-gamma) * theta+h) * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + gamma * theta * fun_pos(x_t[1]-D_t) + c * (S-x_t[m-1])
      }
      else{
        o_vec = sapply(2:m, function(i) rbinom(1, fun_pos(x_t[i] - max(x_t[i-1], D_t)), 1-gamma))
        o_t = sum(o_vec) + fun_pos(x_t[1]-D_t)
        C_t = h * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + theta * o_t + c * (S-x_t[m-1])
      }
    }
    if(case=="LIFO"){
      if(Bino==FALSE){
        o_t = gamma * min(fun_pos(x_t[m]-D_t), x_t[1]) + (1-gamma) * fun_pos(x_t[m]-D_t)
        C_t = ((1-gamma) * theta+h) * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + gamma * theta * min(fun_pos(x_t[m]-D_t), x_t[1]) + c * (S-x_t[m-1])
      }
      else{
        o_vec = sapply(2:m, function(i) rbinom(1, fun_pos(min(x_t[i], fun_pos(x_t[m] - D_t)) - x_t[i-1]), 1-gamma))
        o_t = sum(o_vec) + min(fun_pos(x_t[m]-D_t), x_t[1])
        C_t = h * fun_pos(x_t[m]-D_t) + p * fun_pos(D_t-x_t[m]) + theta * o_t + c * (S-x_t[m-1])
      }
    }
    if(Bino==FALSE){
      for(i in 1:(m-1)){
        if(i == 1){
          x_1 = x_t[i]
        }
        if(case=="FIFO")
          x_t[i] = gamma * fun_pos(x_t[i+1] - max(D_t, x_1))
        if(case=="LIFO")
          x_t[i] = gamma * (min(x_t[i+1], fun_pos(x_t[m]-D_t)) - gamma * min(fun_pos(x_t[m]-D_t), x_1))
      }
    } else {
      x_old = x_t
      for(i in 1:(m-1)){
        if(case=='FIFO'){
          if(i ==1)
            x_t[i] = fun_pos(x_old[i+1] -max(x_old[i], D_t)) - o_vec[i]
          else
            x_t[i] = x_t[i-1] + fun_pos(x_old[i+1] -max(x_old[i], D_t)) - o_vec[i] 
          
        }
        else {
          if(i == 1)
            x_t[i] = fun_pos(min(x_old[i+1], fun_pos(x_old[m] - D_t)) - x_old[i]) - o_vec[i]
          else
            x_t[i] = x_t[i-1] + fun_pos(min(x_old[i+1], fun_pos(x_old[m] - D_t)) - x_old[i]) - o_vec[i]
        }
      }
    }
    C_all = C_all + C_t
  }
  return(C_all/T_)
}

compute_cost <- function(S, c=0, p=8, theta=3, h=1, m=3, gamma=1, pois=5, n_sims=5000, type='poi', case='FIFO', n=0, Ratio=FALSE) {
  
  # Function to perform Monte Carlo simulation
  simulate <- function(type, pois, n_sims, n) {
    if (n == 0) {
      if (type == 'poi') {
        return(rpois(n_sims, pois))
      } else if (type == 'geom') {
        return(rgeom(n_sims, prob = 1 / pois))
      } else {
        stop("Invalid type. Please use 'poi' or 'geom'.")
      }
    } else {
      N_n <- numeric(n_sims)
      for (i in 1:n_sims) {
        N_n[i] <- which.max(sapply(1:(10 * n), function(j) sum(rexp(10 * n, 1)[1:j])) <= n)
      }
      return(sapply(N_n, function(N) sum(rpois(N, pois))))
    }
  }
  
  # Monte Carlo for E[S-D]^+ and E[D-S]^+
  demand_samples <- simulate(type, pois, n_sims, n)
  expected_S_minus_D <- mean(pmax(S - demand_samples, 0))
  expected_D_minus_S <- mean(pmax(demand_samples - S, 0))
  
  # Monte Carlo for the third term
  simulated_term3 <- replicate(n_sims, {
    demand_vector <- simulate(type, pois, m, n)
    weights <- gamma^(0:(m-1))
    
    if (case == 'FIFO') {
      expected_Di <- sum(weights * demand_vector)
      max(S * gamma^(m-1) - expected_Di, 0)
    } else if (case == 'LIFO') {
      shortages <- weights * pmax(S - demand_vector, 0)
      min(shortages)
    } else {
      stop("Invalid case. Please use 'FIFO' or 'LIFO'.")
    }
  })
  
  expected_term3 <- mean(simulated_term3)
  expected_D <- mean(demand_samples)
  
  # Calculate cost
  cost <- (h + (theta + c) * (1 - gamma)) * expected_S_minus_D + (p - c) * expected_D_minus_S + 
    (gamma * (theta + c) / sum(gamma^(0:(m-1)))) * expected_term3 + c * expected_D
  
  diff_part = (gamma * (theta + c) / sum(gamma^(0:(m-1)))) * expected_term3  * (sum(gamma^{1:(m-1)}))
  if(Ratio)
    return(diff_part / cost)
  else
    return(cost)
}


test_n <- function(c=0, p=6, theta=3, h=1, m=3, gamma=1, pois=5, n_sims=5000, type='poi', case='FIFO', n=0) {
  S_range = tune_range(p=p, gamma=gamma, h = h, theta = theta, PP = pois, c = c, type=type)
  min_COST = Inf
  min_S = NA
  for(S in seq(S_range[1],S_range[2],0.1)){
    COST = compute_cost(S, c=c, p=p, theta=theta, h=h, m=m, gamma=gamma, pois=pois, n_sims=n_sims, type=type, case=case, n=n)
    # print(paste0('S = ',S, " COST = ", COST))
    if(COST<min_COST){
      min_COST = COST
      min_S = S
    }
  }
  print(compute_cost(min_S, c=c, p=p, theta=theta, h=h, m=m, gamma=gamma, pois=pois, n_sims=n_sims, type=type, case=case, n=n, Ratio = TRUE))
}

