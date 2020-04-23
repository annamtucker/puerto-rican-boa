library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

A = matrix(c(0.1, 0, 0.6, 1.35,
             0.2, 0.41, 0, 0,
             0, 0.5, 0.54, 0,
             0, 0, 0.18, 0.9), nrow = 4, byrow = T)

N = matrix(0, nrow = 4, ncol = 50)
N[,1] = rep(100, 4)

for(i in 2:ncol(N)){
  N[,i] <- A %*% N[,i-1]
}

matplot(t(N), type ="l")

toplot = as_tibble(expand.grid(stage = c("Y", "J", "S", "A"),
                               year= c(1:50)))

toplot$N = c(N)

ggplot(toplot, aes(x = year, y = N, col = stage)) +
  geom_line(lwd = 2) +
  theme(legend.position = "top",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))+
  xlab("Year") 

## add K
A = matrix(c(0.1, 0, 0.6, 1.35,
             0.2, 0.41, 0, 0,
             0, 0.5, 0.54, 0,
             0, 0, 0.18, 0.9), nrow = 4, byrow = T)

N = matrix(0, nrow = 4, ncol = 50)
N[,1] = rep(100, 4)

K = 450

for(i in 2:ncol(N)){
  mat <- A
  if(sum(N[,i-1]) > K){
    mat[1,3:4] <- c(0, 0)
  }
  N[,i] <- mat %*% N[,i-1]
}

matplot(t(N), type ="l")

toplot = as_tibble(expand.grid(stage = c("Y", "J", "S", "A"),
                               year= c(1:50)))

toplot$N = c(N)

ggplot(toplot, aes(x = year, y = N, col = stage)) +
  geom_line(lwd = 2) +
  theme(legend.position = "top",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))+
  xlab("Year") 

toplot %>% 
  group_by(year) %>% 
  summarize(N = sum(N)) %>% 
  ggplot(aes(x = year, y = N)) +
  geom_line(lwd = 2)+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22)) +
  xlab("Year")


K = 500
N = c(1:1000)
F1 = ifelse(N < K, 5, 0)
F2 = ifelse(N < K, 5*(1-N/K), 0)

plot2 = data.frame(N = N,
                   F1 = F1, 
                   F2 = F2)

plot2 %>% 
  gather(type, val, 2:3) %>% 
  ggplot(aes(x = N, y = val)) +
  geom_line(lwd = 2) +
  facet_wrap(~type, scales = "free") +
  ylab("Fecundity") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))


data.frame(x = rnorm(10000, 19.4, 2.91)) %>% 
  filter(x > 0) %>% 
  ggplot(aes(x = x)) +
  geom_density(fill = "gray80") +
  xlab("Possible clutch size") +
  ylab("Density")+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))


data.frame(x = rnorm(50, 19.4, 2.91),
           yr = c(1:50)) %>% 
  ggplot(aes(x = yr, y = x)) +
  geom_point(size = 4) +
  geom_line(lwd = 1.5) +
  geom_hline(yintercept = 19.4, lty = 2, lwd = 2) +
  xlab("Year") + 
  ylab("Clutch size")+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))


# add para uncertainty, env. stoch, demo stoch.

beta_shape_parms = function(mean, sd){
  avg = mean
  var = sd*sd
  
  #a = ifelse(var < avg*(1-avg), avg*(((avg*(1-avg))/var)-1), 100*avg)
  #b = ifelse(var < avg*(1-avg), (1-avg)*(((avg*(1-avg))/var)-1), 100*(1-avg))
  
  a = 100*avg
  b = 100*(1-avg)
  
  return(list(a = a, b = b))
}

get_beta_vals = function(n, mean, sd){
  x = beta_shape_parms(mean, sd)
  rbeta(n, x$a, x$b)
}

reps = 1000

A = matrix(c(0.1, 0, 0.6, 1.35,
             0.2, 0.41, 0, 0,
             0, 0.5, 0.54, 0,
             0, 0, 0.18, 0.9), nrow = 4, byrow = T)

mat = array(0, dim = c(reps, 4, 4))
for(r in 1:reps){
  mat[r,,] <- matrix(get_beta_vals(16, A, A*0.1))
  mat[r,1,4] <- rlnorm(1, 1.35, 0.1)
}

yr.effect = matrix(runif(reps*50, 0.5, 1.5), nrow = reps, ncol = 50)


N = array(0, dim = c(4, 50, reps))
N[,1,] = round(runif(reps*4, 50, 200))

K = round(runif(reps, 350, 550))

for(r in 1:reps){
  for(i in 2:ncol(N[,,r])){
    
    yrmat <- mat[r,,]*yr.effect[r,i]
    yrmat[1:4,1:2] <- ifelse(yrmat[1:4,1:2] > 1, 1, yrmat[1:4,1:2])  
    yrmat[3:4,3:4] <- ifelse(yrmat[3:4,3:4] > 1, 1, yrmat[3:4,3:4])  
    
    if(sum(N[,i-1,r]) > K[r]){
      yrmat[1,3:4] <- c(0, 0)
    }
    
    N[1,i,r] <- rbinom(1, N[1,i-1,r], yrmat[1,1]) + rpois(1, N[3,i-1,r]*yrmat[1,3]) + rpois(1, N[4,i-1,r]*yrmat[1,4])
    N[2,i,r] <- rbinom(1, N[1,i-1,r], yrmat[2,1]) + rbinom(1, N[2,i-1,r], yrmat[2,2])
    N[3,i,r] <- rbinom(1, N[2,i-1,r], yrmat[3,2]) + rbinom(1, N[3,i-1,r], yrmat[3,3])
    N[4,i,r] <- rbinom(1, N[3,i-1,r], yrmat[4,3]) + rbinom(1, N[4,i-1,r], yrmat[4,4])
    
  }
}



toplot = as_tibble(expand.grid(stage = c("Y", "J", "S", "A"),
                               year= c(1:50),
                               rep = c(1:reps)))

toplot$N = c(N)



toplot %>% 
  group_by(year, rep) %>% 
  summarize(N = sum(N)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  summarize(medN = median(N),
            lcl = quantile(N, 0.025),
            ucl = quantile(N, 0.975)) %>% 
  ggplot(aes(x = year, y = medN)) +
  geom_line(lwd = 2)+
  geom_line(lwd = 1, lty = 2, aes(y = lcl)) +
  geom_line(lwd = 1, lty = 2, aes(y = ucl)) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22)) +
  xlab("Year") +
  ylab("N")





