# demographic projection without loops
# 16 april 2019

source("scripts/packages-functions.R")
source("scripts/projection-functions-2.R")

set.seed(419)

# simulation parameters ----
n.iters = 10000               # for parametric uncertainty
n.reps = 1                    # replicates of each parameter draw
n.years = 30                  # number of years to project 
n.habs = 2

# scenarios 
prop.urban = 0.43
habitat.area = 379029
urb.rates = c(0, 0.008, 0.016, 0.024)  # rate of urban growth per year


# demographic parameters ----

# initial population size
Ninit.min = habitat.area*0.1
Ninit.max =  habitat.area*0.5

# maximum population size
K.min = habitat.area*1
K.max = habitat.area*3

# demographic rate CVs
para.cv = 0.15
temp.cv = 0.15

# adjustment on average demographic rates for urban habitats
urban = c(0.5, 1)

# survival rates
young.surv = 0.3
juv.surv = 0.9
subad.surv = 0.72
ad.surv = 0.9

# growth rates
yj = 0.67
js = 0.55
sa = 0.25

# fecundity and breeding propensity of subadults
f = 4.5
subad.bp = 2/f
max.fec = 35

# matrix parameters
T.YY = young.surv * (1-yj)
T.JY = young.surv * yj
T.JJ = juv.surv * (1-js)
T.SJ = juv.surv * js
T.SS = subad.surv * (1-sa)
T.AS = subad.surv * sa
T.AA = ad.surv
F.SY = f * subad.bp * young.surv
F.AY = f * young.surv

avg.mat = matrix(c(T.YY, 0, F.SY, F.AY,
               T.JY, T.JJ, 0, 0,
               0, T.SJ, T.SS, 0,
               0, 0, T.AS, T.AA),
             nrow = 4, ncol = 4, byrow = T)

# parametric uncertainty SD
para.sd.mat = avg.mat*para.cv

# manually change sd for subad fecundity (less than 1 otherwise and leads to error)
#para.sd.mat[1,c(3,4)] <- c(log(F.SY)*0.15, log(F.AY)*0.15)

# to isolate survival rates only
surv.mask = matrix(c(1, 0, 0, 0,
                     1, 1, 0, 0,
                     0, 1, 1, 0,
                     0, 0, 1, 1),
                   nrow = 4, ncol = 4, byrow = T)

# to isolate fecundities only
fec.mask = matrix(c(0, 0, 1, 1, 
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0),
                  nrow = 4, ncol = 4, byrow = T)


# run simulation ----

iters <- as_tibble(expand.grid(rate = urb.rates,
                               prop.hab = prop.urban,
                               iter = c(1:n.iters))) %>% 
  mutate(parms = purrr::map(iter, draw_parms)) %>% 
  unnest()

sim <- iters %>% 
  mutate(projection = pmap(list(avg, prop.hab, rate, Ninit, K), project_pop)) 

saveRDS(sim, file = "output/final-projection.rds")


pushover("pr boa projection done")

#sim <- readRDS("output/final-projection.rds")


Ntot <- sim %>% 
  select(iter, rate, projection, Ninit, K) %>% 
  unnest() 




# plot outcomes ----

percent = gsub(" ", "", paste(urb.rates*1000, "%"))

scenario.labs = sapply(percent, FUN = function(x) paste(x, "urban growth\nper decade"))

labels <- tibble(rate = urb.rates,
                scenario = scenario.labs) %>% 
  mutate(scenario = fct_reorder(scenario, rate))

# population size over time
Nplot <- Ntot %>% 
  full_join(labels) %>% 
  group_by(scenario, iter, year, Ninit) %>% 
  summarize(Ntot = sum(Ntot)) %>% 
  ungroup() %>% 
  mutate(deltaN = Ntot-Ninit) %>% 
  group_by(scenario, year) %>% 
  summarize(medN = median(deltaN, na.rm = T),
            lciN = quantile(deltaN, 0.025, na.rm = T),
            uciN = quantile(deltaN, 0.975, na.rm = T)) %>% 
  ungroup() %>% 
  #mutate(scenariolab = scenario.labs[scenario]) %>% 
  ggplot(aes(x = year, y = medN/100000)) +
  geom_ribbon(aes(ymin = lciN/100000, ymax = uciN/100000), fill = "gray80", alpha = 0.7) +
  geom_line(lwd = 1.5) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Change in population size\n(x100,000)") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(1, 30)) +
  facet_wrap(~scenario, scales = "free") +
  scale_y_continuous(limits = c(NA, 6), breaks = c(-1, 0, 1, 2, 3, 4, 5)) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, margin = margin(2,2,2,2)),
        panel.spacing = unit(2, "lines"))

Nplot

save_plot(Nplot, file = "output/projection-final.jpg", base_width = 8, base_height = 6)



# overall probabilities
probs <- Ntot %>%
  rename(scenario = rate) %>% 
  group_by(scenario, iter, year) %>% 
  summarize(N = sum(Ntot)) %>% 
  ungroup() %>% 
  group_by(scenario, iter) %>% 
  nest() %>% 
  mutate(lambda = map(data, calc_lambda)) %>% 
  unnest() %>% 
  group_by(scenario, iter) %>% 
  summarize(inc = ifelse(geo_mean(lambda, N) >= 1, 1, 0),
            dec = ifelse(geo_mean(lambda, N) < 1, 1, 0),
            qe_50 = ifelse(any(N < 50), 1, 0),
            qe_500 = ifelse(any(N < 500), 1, 0),
            qe_1000 = ifelse(any(N < 1000), 1, 0),
            qe_5000 = ifelse(any(N < 5000), 1, 0)) %>% 
  ungroup() %>% 
  #mutate(Scenario = scenario.labs[scenario]) %>% 
  group_by(scenario) %>% 
  summarize(pext_50 = mean(qe_50),
            pext_500 = mean(qe_500),
            pext_1000 = mean(qe_1000),
            pext_5000 = mean(qe_5000),
            pinc = mean(inc),
            pdec = mean(dec))
probs

write_csv(probs, "output/scenario_probabilities-final.csv")
   

# quasi-extinction probability
qe_plot <- probs %>% 
  gather(prob, val, 2:5) %>% 
  mutate(prob = fct_relevel(prob, "pext_50", "pext_500", "pext_1000", "pext_5000")) %>% 
  ggplot(aes(x = scenario*1000, y = val, color = prob)) +
  geom_path(lty = 2, lwd = 1, position = position_dodge(width = 0.75), alpha = 0.8) +
  geom_point(size = 3, position = position_dodge(width = 0.75)) +
  xlab("Urbanization rate (% growth per decade)") +
  ylab("Quasi-extinction probability") +
  scale_x_continuous(breaks = c(0, 8, 16, 24)) +
  ylim(0, 0.02) +
  scale_color_viridis_d(name = "Quasi-extinction threshold",
                        labels = c(50, 500, 1000, 5000),
                        option = "D", 
                        end = 0.8) +
  theme(legend.position = "top")
qe_plot

save_plot(qe_plot, file = "output/qe_probs-final.jpg",
          base_width = 8, base_height = 5)



### average lambda and percent change in pop size for each scenario ----

pop_change <- Ntot %>% 
  group_by(iter, rate) %>% 
  summarize(mean_lam = geo_mean(lambda, Ntot),
         percent_change = calc_p_change(Ntot)) %>% 
  ungroup()

pop_change %>% 
  gather(var, val, 3:4) %>% 
  group_by(rate, var) %>% 
  summarize(med = median(val, na.rm = T),
            lci = quantile(val, 0.025, na.rm = T),
            uci = quantile(val, 0.975, na.rm = T)) %>% 
  filter(var == "mean_lam")

  
ggplot(pop_change, aes(x = mean_lam)) +
  geom_histogram(col = "gray") +
  facet_wrap(~rate)

ggplot(pop_change, aes(x = percent_change)) +
  geom_histogram(col = "gray") +
  facet_wrap(~rate)


Ntot %>% 
  group_by(rate, iter, year, Ninit, Ntot) %>% 
  ungroup() %>% 
  mutate(deltaN = Ntot-Ninit) %>% 
  group_by(rate, year) %>% 
  summarize(medN = median(deltaN),
            lciN = quantile(deltaN, 0.025),
            uciN = quantile(deltaN, 0.975)) 

  ungroup() %>% 
  mutate(percent_change = medN - Ninit)
