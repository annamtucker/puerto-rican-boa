

out <- readRDS("output/eval-Ninit.rds")
head(out)

Ntot <- out %>% 
  select(iter, rate, projection, Ninit, K) %>% 
  unnest() 

qe.probs <- Ntot %>% 
  filter(year == 30) %>% 
  group_by(Ninit) %>% 
  summarize(`50` = mean(Ntot < 50),
            `100` = mean(Ntot < 100),
            `500` = mean(Ntot < 500),
            `1000` = mean(Ntot < 1000),
            `5000` = mean(Ntot < 5000)) %>% 
  ungroup() %>% 
  gather(threshold, prob, 2:6) %>% 
  mutate(threshold = as.numeric(threshold)) 

all.probs <- as_tibble(expand.grid(prob = seq(0, 1, 0.05),
                                   threshold = c(50, 100, 500, 1000, 5000)))

all.probs$Ninit = NA

for(i in 1:nrow(all.probs)){
  all.probs$Ninit[i] <- min(qe.probs$Ninit[qe.probs$threshold == all.probs$threshold[i] &
                                          qe.probs$prob <= all.probs$prob[i]])
}

all.probs %>% 
  mutate(threshold.fct = as.factor(threshold), levels = threshold) %>% 
  ggplot(aes(x = prob, y = threshold.fct)) +
  geom_tile(aes(fill = Ninit)) +
  geom_label(aes(label = Ninit))


qe.probs %>% 
  mutate(threshold.fct = as.factor(threshold), levels = threshold) %>% 
  ggplot(aes(x = Ninit, y = threshold.fct)) +
  geom_tile(aes(fill = prob)) +
  scale_fill_viridis_c()

write.csv(all.probs, file = "output/eval_Ninit_qe.csv", row.names = F)
write.csv(qe.probs, file = "output/all_qe-prob_Ninit.csv", row.names = F)
