
### compare effects of reducing maximum K and/or maximum initial N


sim1 <- readRDS("output/Ninit-max-1_K-max-3.rds")
sim2 <- readRDS("output/Ninit-max-0.5_K-max-3.rds")
sim3 <- readRDS("output/Ninit-max-1_K-max-6.rds")
sim4 <- readRDS("output/Ninit-max-0.5_K-max-6.rds")


Ntot1 <- sim1 %>% 
  select(iter, rate, projection, Ninit, K) %>% 
  unnest() %>% 
  mutate(maxInitDens = 1,
         maxKDens = 3,
         sim = "Max. initial N = 379,029\nMax. K = 1,137,087")

Ntot2 <- sim2 %>% 
  select(iter, rate, projection, Ninit, K) %>% 
  unnest() %>% 
  mutate(maxInitDens = 0.5,
         maxKDens = 6,
         sim = "Max. initial N = 189,515\nMax. K = 2,274,174")

Ntot3 <- sim3 %>% 
  select(iter, rate, projection, Ninit, K) %>% 
  unnest() %>% 
  mutate(maxInitDens = 0.5,
         maxKDens = 3,
         sim = "Max. initial N = 189,515\nMax. K = 1,137,087")

Ntot4 <- sim4 %>% 
  select(iter, rate, projection, Ninit, K) %>% 
  unnest() %>% 
  mutate(maxInitDens = 1,
         maxKDens = 6,
         sim = "Max. initial N = 379,029\nMax. K = 2,274,174")

Ntot <- Ntot1 %>% 
  bind_rows(Ntot2) %>% 
  bind_rows(Ntot3) %>% 
  bind_rows(Ntot4)


# population size over time
Nplot <- Ntot %>% 
  filter(rate == 0) %>% 
  mutate(deltaN = Ntot-Ninit) %>% 
  group_by(year, sim) %>% 
  summarize(medN = median(deltaN, na.rm = T),
            lciN = quantile(deltaN, 0.025, na.rm = T),
            uciN = quantile(deltaN, 0.975, na.rm = T)) %>% 
  ungroup() %>% 
  filter(!is.na(sim)) %>% 
  ggplot(aes(x = year, y = medN/1000000, col = sim, fill = sim)) +
  geom_ribbon(aes(ymin = lciN/1000000, ymax = uciN/1000000), alpha = 0.5) +
  geom_line(lwd = 1.5, alpha = 0.9) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Change in population size\n(millions)") +
  xlab("Year") +
  facet_wrap(~sim, scales = "free") +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(1, 30)) +
  scale_y_continuous(limits = c(NA, 2)) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, margin = margin(2,2,2,2)),
        strip.text.y = element_text(angle = 0),
        panel.spacing = unit(2, "lines"),
        legend.position = "none")
Nplot

save_plot(Nplot, file = "output/reduce_max_K_Ninit.jpg", base_width = 8, base_height = 6)

 
 Nplot_med <- Ntot %>% 
   full_join(labels) %>% 
   group_by(scenario, iter, year, Ninit, sim) %>% 
   summarize(Ntot = sum(Ntot)) %>% 
   ungroup() %>% 
   mutate(deltaN = Ntot-Ninit) %>% 
   group_by(scenario, year, sim) %>% 
   summarize(medN = median(deltaN),
             lciN = quantile(deltaN, 0.025),
             uciN = quantile(deltaN, 0.975)) %>% 
   ungroup() %>% 
   #mutate(scenariolab = scenario.labs[scenario]) %>% 
   ggplot(aes(x = year, y = medN/1000000, col = sim, fill = sim)) +
   geom_line(lwd = 1.5, alpha = 0.9) +
   geom_hline(yintercept = 0, lty = 2) +
   ylab("Median change in population size\n(millions)") +
   xlab("Year") +
   facet_wrap(~scenario, scales = "free") +
   scale_x_continuous(breaks = seq(0, 30, 5), limits = c(1, 30)) +
   scale_y_continuous(limits = c(NA, 1)) +
   theme(strip.background = element_rect(fill = "white"),
         strip.text = element_text(size = 14, margin = margin(2,2,2,2)),
         panel.spacing = unit(2, "lines"),
         legend.position = "top") +
   guides(col = guide_legend(override.aes = list(alpha = 1, lwd = 3))) +
   scale_color_discrete(name = "")
 Nplot_med

 save_plot(Nplot_med, file = "output/reduce_max_K_Ninit_medians.jpg", base_width = 10, base_height = 8)
 
 
 probs <- Ntot %>%
   filter(rate == 0) %>% 
   group_by(iter, year, sim, maxInitDens, maxKDens) %>% 
   summarize(N = sum(Ntot)) %>% 
   ungroup() %>% 
   group_by(iter, sim, maxInitDens, maxKDens) %>% 
   nest() %>% 
   mutate(lambda = map(data, calc_lambda)) %>% 
   unnest() %>% 
   group_by(iter, sim, maxInitDens, maxKDens) %>% 
   summarize(inc = ifelse(geo_mean(lambda, N) >= 1, 1, 0),
             dec = ifelse(geo_mean(lambda, N) < 1, 1, 0),
             qe_50 = ifelse(any(N < 50), 1, 0),
             qe_500 = ifelse(any(N < 500), 1, 0),
             qe_1000 = ifelse(any(N < 1000), 1, 0),
             qe_5000 = ifelse(any(N < 5000), 1, 0)) %>% 
   ungroup() %>% 
   group_by(sim, maxInitDens, maxKDens) %>% 
   summarize(pext_50 = mean(qe_50),
             pext_500 = mean(qe_500),
             pext_1000 = mean(qe_1000),
             pext_5000 = mean(qe_5000),
             pinc = mean(inc),
             pdec = mean(dec))
 probs
 
 write.csv(probs, "output/probs_reduce_max_K_Ninit.csv", row.names = F)
 
 
 probs = read_csv("output/probs_reduce_max_K_Ninit.csv")
 
 
 qe_plot <- probs %>% 
   gather(prob, val, 3:6) %>% 
   mutate(prob = fct_relevel(prob, "pext_50", "pext_500", "pext_1000", "pext_5000")) %>% 
   ggplot(aes(x = scenario*1000, y = val, color = sim)) +
   geom_path(lty = 2, lwd = 1, position = position_dodge(width = 0.75), alpha = 0.8) +
   geom_point(size = 3, position = position_dodge(width = 0.75)) +
   xlab("Urbanization rate (% growth per decade)") +
   ylab("Quasi-extinction probability") +
   scale_x_continuous(breaks = c(0, 8, 16, 24)) +
   # scale_color_viridis_d(name = "Quasi-extinction threshold",
   #                       labels = c(50, 500, 1000, 5000),
   #                       option = "D", 
   #                       end = 0.8) +
   theme(legend.position = "top",
         strip.background = element_rect(fill = "white")) +
   facet_wrap(~prob)
 qe_plot
 
 
 
 
 
 
 Ntot %>% 
   filter(rate == 0 & year == 30) %>% 
   mutate(qext = ifelse(Ntot < 5000, 1, 0),
          inc = ifelse(Ntot >= Ninit, 1, 0)) %>% 
   select(iter, year, Ninit, K, maxInitDens, maxKDens, Ntot, qext, inc) %>% 
   mutate(Kint = cut_interval(K, length = 100000),
          Kint = fct_reorder(Kint, K)) %>%
   group_by(Kint, maxKDens) %>% 
   summarize(midpointK = mean(K),
             pext = mean(qext),
             pinc = mean(inc)) %>% 
   ungroup() %>% 
   gather(type, val, pext, pinc) %>% 
   mutate(typelab = c("pext" = "Quasi-extinction (<5000 individuals)",
                      "pinc" = "Stability or growth")[type]) %>% 
   ggplot(aes(x = midpointK, y = val, col = as.character(maxKDens))) +
   geom_smooth(se = F) +
   facet_wrap(~typelab)
 
 
 Ntot %>% 
   filter(rate == 0 & year == 30) %>% 
   mutate(qext = ifelse(Ntot < 5000, 1, 0),
          inc = ifelse(Ntot >= Ninit, 1, 0)) %>% 
   select(iter, year, Ninit, K, maxInitDens, maxKDens, Ntot, qext, inc) %>% 
   mutate(minNint = cut_interval(Ninit, length = 100000)) %>%
   group_by(minNint, maxInitDens) %>% 
   summarize(midpoint = mean(Ninit),
             pext = mean(qext),
             pinc = mean(inc)) %>% 
   ungroup() %>% 
   gather(type, val, pext, pinc) %>% 
   mutate(typelab = c("pext" = "Quasi-extinction (<5000 individuals)",
                      "pinc" = "Stability or growth")[type]) %>% 
   ggplot(aes(x = midpoint, y = val, col = as.character(maxInitDens))) +
   geom_smooth(se = F) +
   facet_wrap(~typelab)
 