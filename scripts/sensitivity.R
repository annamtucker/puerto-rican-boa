# simple projection for sensitivity analysis

source("scripts/packages-functions.R")

# inputs ----

n.iters = 100    # number of random parameter draws
n.reps = 1000     # reps per parameter draw
n.years = 30
qe.threshold = 5000

filename = "output/sensitivity_Ninit_new.rds"

# demographic rate ranges
#as.range = c(0.01, 0.99)
#ss.range = c(0.01, 0.99)
#js.range = c(0.01, 0.99)
#ys.range = c(0.01, 0.99)
#fec.range = c(0.1, 10)
 
as.range = c(0.9, 0.9)
ss.range = c(0.72, 0.72)
js.range = c(0.9, 0.9)
ys.range = c(0.3, 0.3)
fec.range = c(4.5, 4.5)


# initial population size
#Ninit.min = 113709
#Ninit.max = 113709

#Ninit.min = 37903
#Ninit.max =  189515
#Ninit.max = 379029

Ninit.min = 10000
Ninit.max = 1000000


# maximum population size
K.min = 758058
K.max = 758058 

#K.min = 379029
#K.max = 1137087
#K.max = 2274074

# growth estimates
yj = 0.67
js = 0.55
sa = 0.25

# fecundity estimates
subad.bp = 0.15
fec.sd = 1.3

# demographic rate CVs
para.cv = 0.15
temp.cv = 0.15

# simulation ----

surv.mask = fec.mask = array(dim = c(n.iters, 4,4))
for(i in 1:n.iters){
  surv.mask[i,,] <- matrix(c(1, 0, 0, 0,
                       1, 1, 0, 0,
                       0, 1, 1, 0,
                       0, 0, 1, 1),
                     nrow = 4, ncol = 4, byrow = T)

  fec.mask[i,,] <- matrix(c(0, 0, 1, 1, 
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0),
                    nrow = 4, ncol = 4, byrow = T)
}


# draw parameter values
f = runif(n.iters, fec.range[1], fec.range[2])
young.surv = runif(n.iters, ys.range[1], ys.range[2])
juv.surv = runif(n.iters, js.range[1], js.range[2])
subad.surv = runif(n.iters, ss.range[1], ss.range[2])
ad.surv = runif(n.iters, as.range[1], as.range[2])

# matrix parameters
T.YY = young.surv * (1-yj)
T.JY = young.surv * yj
T.JJ = juv.surv * (1-js)
T.SJ = juv.surv * js
T.SS = subad.surv * (1-sa)
T.AS = subad.surv * sa
T.AA = ad.surv
F.S = f * subad.bp
F.A = f

avg.mat = array(dim = c(n.iters, 4, 4))
for(i in 1:n.iters){
  avg.mat[i,,] <- matrix(c(T.YY[i], 0, F.S[i], F.A[i],
                           T.JY[i], T.JJ[i], 0, 0,
                           0, T.SJ[i], T.SS[i], 0,
                           0, 0, T.AS[i], T.AA[i]),
                         nrow = 4, ncol = 4, byrow = T)
}

# parametric uncertainty SD
para.sd.mat = avg.mat*para.cv


# replicate-level average and SD for each rates

# survival and transition rates
surv = avg.mat*surv.mask
surv.sd = para.sd.mat*surv.mask

# fecundity
fec = avg.mat*fec.mask
fec.sd = para.sd.mat*fec.mask


avg.surv = sd.surv = array(dim = c(n.iters, n.reps, 4,4))
avg.fec = sd.fec = array(dim = c(n.iters, n.reps, 4, 4))
for(i in 1:n.iters){
  for(r in 1:n.reps){
    avg.surv[i,r,,] <-  matrix(get_beta_vals(16, surv[i,,], surv.sd[i,,]),
                                nrow = 4, ncol = 4)
    sd.surv[i,r,,] <- matrix(rinvgauss(16, surv.sd[i,,], 1),
                             nrow = 4, ncol = 4)
    
    avg.fec[i,r,,] <- matrix(rlnorm(16, log(fec[i,,]), fec.sd[i,,]),
                     nrow = 4, ncol = 4)
    sd.fec[i,r,,] = matrix(rinvgauss(16, fec.sd[i,,], 1),
                    nrow = 4, ncol = 4)
  }
}

sd.surv[is.na(sd.surv)] <- 0
sd.fec[is.na(sd.fec)] <- 0

# average matrices for each rep
avg = avg.surv + avg.fec

# stable stage distribution for each rep
dist = array(dim = c(n.iters, n.reps, 4))
for(i in 1:n.iters){
  for(r in 1:n.reps){
    dist[i,r,] <-  eigen.analysis(avg[i,r,,])$stable.stage
  }
}

# initial population size
Ninit = matrix(round(runif(n.iters*n.reps, Ninit.min, Ninit.max)),
               ncol = n.reps, nrow = n.iters)

# carrying capacity
K = matrix(round(runif(n.iters*n.reps, K.min, K.max)),
               ncol = n.reps, nrow = n.iters)


# year-specific parameters
surv.sd.yr = avg.surv*temp.cv

yr.surv = array(dim = c(n.iters, n.reps, n.years, 4, 4))
yr.fec = array(dim = c(n.iters, n.reps, n.years, 4, 4))
for(i in 1:n.iters){
  for(r in 1:n.reps){
    for(t in 1:n.years){
      yr.surv[i,r,t,,] <- matrix(get_beta_vals(16, avg.surv[i,r,,], surv.sd.yr[i,r,,]),
                               nrow = 4, ncol = 4)
      yr.fec[i,r,t,,] <-  matrix(rpois(16, avg.fec[i,r,,]),
                               nrow = 4, ncol = 4)
    }
  }
}

# iter-, rep- and year-specific matrix
mat <- yr.surv + yr.fec

# stage-specific N
N = array(dim = c(n.iters, n.reps, n.years, 4))
for(i in 1:n.iters){
  for(r in 1:n.reps){
    N[i,r,1,] <- round(Ninit[i,r]*dist[i,r,])
  }
}

# place to store total N
Ntot = array(dim = c(n.iters, n.reps, n.years))


for(i in 1:n.iters){
  for(r in 1:n.reps){
    for(t in 2:n.years){
      
      # set fecundity = 0 if pop exceeds K
      if(sum(N[i,r,t-1,]) > K[i,r]){
        mat[i,r,t-1,1,3:4] <- 0
      }
      
      N[i,r,t,] <- round(mat[i,r,t-1,,] %*% N[i,r,t-1,])
    }
    
    Ntot[i,r,] = apply(N[i,r,,], 1, sum)
      
  }
}


saveRDS(list(f, young.surv, juv.surv, subad.surv, ad.surv, Ntot), filename)

pushover("done")

#########################



all.res <- readRDS("output/sensitivity_results.rds")
names(all.res) <- c("f", "young.surv", "juv.surv", "subad.surv", "ad.surv", "Ntot")
list2env(all.res, globalenv())
n.iters = dim(Ntot)[1]
n.reps = dim(Ntot)[2]
n.years = 30

# summarize results ----

res <- tibble(iter = 1:n.iters,
             fec = f,
             young.surv = young.surv,
             juv.surv = juv.surv,
             subad.surv = subad.surv,
             ad.surv = ad.surv) %>% 
  full_join(as_tibble(expand.grid(iter = 1:n.iters,
                                  rep = 1:n.reps)))

Nres <- as_tibble(expand.grid(iter = 1:n.iters,
                              rep = 1:n.reps,
                              year = 1:n.years)) %>% 
  mutate(N = c(Ntot))

res <- full_join(res, Nres)

Nfec <- res %>% 
  group_by(iter, fec, year) %>% 
  summarize(medN = median(N)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = medN)) +
  geom_line(alpha = 0.5, aes(col = fec, group = iter)) +
  scale_color_viridis_c() +
  xlab("Year") +
  ylab("Median N")

Nys <-  res %>% 
  group_by(iter, young.surv, year) %>% 
  summarize(medN = median(N)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = medN)) +
  geom_line(alpha = 0.5, aes(col = young.surv, group = iter)) +
  scale_color_viridis_c()  +
  xlab("Year") +
  ylab("Median N")

Njs <-  res %>% 
  group_by(iter, juv.surv, year) %>% 
  summarize(medN = median(N)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = medN)) +
  geom_line(alpha = 0.5, aes(col = juv.surv, group = iter)) +
  scale_color_viridis_c()  +
  xlab("Year") +
  ylab("Median N")

Nss <-  res %>% 
  group_by(iter, subad.surv, year) %>% 
  summarize(medN = median(N)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = medN)) +
  geom_line(alpha = 0.5, aes(col = subad.surv, group = iter)) +
  scale_color_viridis_c()  +
  xlab("Year") +
  ylab("Median N")

Nas <-  res %>% 
  group_by(iter, ad.surv, year) %>% 
  summarize(medN = median(N)) %>% 
  ungroup() %>% 
  ggplot(aes(x = year, y = medN)) +
  geom_line(alpha = 0.5, aes(col = ad.surv, group = iter)) +
  scale_color_viridis_c()  +
  xlab("Year") +
  ylab("Median N")

Nplot <- plot_grid(Nfec, Nys, Njs, Nss, Nas,
                   labels = c("Fecundity", "Young Survival", "Juvenile Survival",
                              "Subadult Survival", "Adult Survival"),
                   nrow = 2, align = "hv")
Nplot

# quasi-extinction probability

# input rates

rates <- tibble(parm = c("fec", "young.surv", "juv.surv", "subad.surv", "ad.surv"),
                input = c(9, 0.3, 0.9, 0.72, 0.9))

# survival rates
qe.threshold = 5000

res.ext <- res %>% 
  filter(year == 30) %>% 
  mutate(ext = ifelse(N < qe.threshold, 1, 0)) %>% 
  group_by(iter, fec, young.surv, juv.surv, subad.surv, ad.surv) %>%
  summarize(pext = mean(ext)) %>% 
  ungroup()

pext_plot <- res.ext %>% 
  gather(parm, val, 2:6) %>% 
  full_join(rates) %>% 
  mutate(parmlab = fct_relevel(c("fec" = "Fecundity", "young.surv" = "Young Survival",
                     "juv.surv" = "Juvenile Survival", "subad.surv" = "Subadult Survival",
                     "ad.surv" = "Adult Survival")[parm],
         c("Fecundity", "Young Survival", "Juvenile Survival", "Subadult Survival", "Adult Survival"))) %>% 
  ggplot(aes(x = val, y = pext)) +
  geom_point(alpha = 0.1) +
  stat_smooth(se = F, lwd = 1, col = "black") +
  geom_vline(aes(xintercept = input), lty = 2) +
  facet_wrap(~parmlab, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(2,2,2,2))) +
  ylim(0,1) +
  xlab("Value") +
  ylab("Quasi-extinction probability \n (Initial N = 215,000)")
pext_plot

save_plot(pext_plot, file = "output/sens_pext.jpg",
          base_width = 10, base_height = 6)


# what is the minimum adult survival above which all pext = 0
res.ext %>% 
  filter(ad.surv > 0.83) %>% 
  pull(pext) %>% 
  sum()
  


res.lam <- res %>% 
  group_by(iter, rep) %>% 
  nest() %>% 
  mutate(lambda = map(data, calc_lambda)) %>% 
  unnest() 

lam_plot <- res.lam %>% 
  group_by(iter, rep, fec, young.surv, juv.surv, subad.surv, ad.surv) %>% 
  summarize(mean.lambda = geo_mean(lambda, N)) %>% 
  ungroup() %>% 
  group_by(iter, fec, young.surv, juv.surv, subad.surv, ad.surv) %>% 
  summarize(lcl = quantile(mean.lambda, 0.025),
            med = median(mean.lambda),
            ucl = quantile(mean.lambda, 0.975)) %>% 
  ungroup() %>% 
  gather(parm, val, 2:6) %>% 
  full_join(rates) %>% 
  mutate(parmlab = fct_relevel(c("fec" = "Fecundity", "young.surv" = "Young Survival",
                                 "juv.surv" = "Juvenile Survival", "subad.surv" = "Subadult Survival",
                                 "ad.surv" = "Adult Survival")[parm],
                               c("Fecundity", "Young Survival", "Juvenile Survival", "Subadult Survival", "Adult Survival"))) %>% 
  ggplot(aes(x = val)) +
  geom_linerange(aes(ymin = lcl, ymax = ucl), alpha = 0.1) +
  geom_point(aes(y = med), alpha = 0.3) +
  geom_vline(aes(xintercept = input), lty = 2) +
  facet_wrap(~parmlab, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(2,2,2,2))) +
  xlab("Value") +
  ylab("Average population growth rate (\u03bb)")
lam_plot  

save_plot(lam_plot, file = "output/sens_lambda.jpg",
          base_width = 10, base_height = 6)


# regression ----
res.mod <- res.lam %>% 
  group_by(iter, rep, fec, young.surv, juv.surv, subad.surv, ad.surv) %>% 
  summarize(mean.lambda = geo_mean(lambda, N)) %>% 
  ungroup() 
  

mod <- lm(mean.lambda ~ fec + young.surv + juv.surv + subad.surv + ad.surv - 1,
          data = res.mod)
summary(mod)
plot(resid(mod))
sens.mod <- tidy(mod)

sens.mod %>% 
  write_csv("output/sens_ests.csv")


# probability of increase/decline ----
lam.prob <- res.lam %>%   
  group_by(iter, rep, fec, young.surv, juv.surv, subad.surv, ad.surv) %>% 
  summarize(mean.lambda = geo_mean(lambda, N)) %>% 
  ungroup() %>% 
  mutate(incr = ifelse(mean.lambda > 1, 1, 0),
         stable = ifelse(mean.lambda == 1, 1, 0),
         decr = ifelse(mean.lambda < 1, 1, 0)) %>% 
  group_by(iter, fec, young.surv, juv.surv, subad.surv, ad.surv) %>% 
  summarize(p_incr = mean(incr),
            p_stable = mean(stable),
            p_decr = mean(decr)) %>% 
  ungroup()
  
lam_incr <- lam.prob %>% 
  gather(parm, val, 2:6) %>% 
  full_join(rates) %>% 
  mutate(parmlab = fct_relevel(c("fec" = "Fecundity", "young.surv" = "Young Survival",
                                 "juv.surv" = "Juvenile Survival", "subad.surv" = "Subadult Survival",
                                 "ad.surv" = "Adult Survival")[parm],
                               c("Fecundity", "Young Survival", "Juvenile Survival", "Subadult Survival", "Adult Survival"))) %>% 
  ggplot(aes(x = val, y = p_incr)) +
  geom_point(alpha = 0.1) +
  stat_smooth(se = F, lwd = 1, col = "black") +
  geom_vline(aes(xintercept = input), lty = 2) +
  facet_wrap(~parmlab, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(2,2,2,2))) +
  xlab("Value") +
  ylab("Probability of population increase \n (\u03bb > 1)")
lam_incr

save_plot(lam_incr, file = "output/sens_pincr.jpg",
          base_width = 10, base_height = 6)
  

lam_decr <- lam.prob %>% 
  gather(parm, val, 2:6) %>% 
  full_join(rates) %>% 
  mutate(parmlab = fct_relevel(c("fec" = "Fecundity", "young.surv" = "Young Survival",
                                 "juv.surv" = "Juvenile Survival", "subad.surv" = "Subadult Survival",
                                 "ad.surv" = "Adult Survival")[parm],
                               c("Fecundity", "Young Survival", "Juvenile Survival", "Subadult Survival", "Adult Survival"))) %>% 
  ggplot(aes(x = val, y = p_decr)) +
  geom_point(alpha = 0.1) +
  stat_smooth(se = F, lwd = 1, col = "black") +
  geom_vline(aes(xintercept = input), lty = 2) +
  facet_wrap(~parmlab, scales = "free") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(margin = margin(2,2,2,2))) +
  xlab("Value") +
  ylab("Probability of population decline \n (\u03bb < 1)")
lam_decr

save_plot(lam_decr, file = "output/sens_pdecr.jpg",
          base_width = 10, base_height = 6)
