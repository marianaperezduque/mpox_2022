library(ggplot2)
library(rstan)
library(cmdstanr)
library(janitor)
library(tidyverse)
library(readxl)

# read in data

mpx <- read_delim("mpx_clean.csv")

count_data <- mpx %>% 
  group_by(data_do_inicio_dos_sintomas, travel) %>% 
  summarize(n_cases = dplyr::n()) %>% 
  ungroup()

count_data$travel[is.na(count_data$travel)] = 0


# pop 2021 PT census
pop <- read_csv("pop.csv", col_names = FALSE)
colnames(pop) <- c("age_group", "n")

# plot epicurve
df <- count_data %>% drop_na()
df$travel <- as.factor(df$travel)
ggplot(df, aes(data_do_inicio_dos_sintomas, n_cases, fill = travel))+ geom_col() + theme_classic()

df$travel <- as.numeric(df$travel)

# dates
dates <- seq.Date(from=min(df$data_do_inicio_dos_sintomas),
                  to=max(df$data_do_inicio_dos_sintomas), by='day')

# reformat data for seir model (incidence object)
ndf <- data.frame(date=dates, cases=NA)
for(i in 1:nrow(ndf)){
  ndf$cases[i] <- sum(df$n_cases[df$data_do_inicio_dos_sintomas==ndf$date[i]])
  ndf$imports[i] <- sum(df$travel[df$data_do_inicio_dos_sintomas==ndf$date[i]])
}


# model inputs
data <- list()
data$Nt <- nrow(ndf) # N time points
data$cases <- ndf$cases # cases time series
data$pop <- 107328 #sum(pop$n) # population size https://static-content.springer.com/esm/art%3A10.1186%2F1471-2458-13-919/MediaObjects/12889_2012_5864_MOESM1_ESM.pdf
data$seed <- 15/data$pop # seeded infections
data$alpha <- 1/5.6 # 5.6 day incubation period
data$sigma <- 1/21 # 21 day infectious period
data$imports <- ndf$imports/data$pop


# fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path('/Users/marianaperezduque/.cmdstan/cmdstan-2.32.1')
setwd('~/Documents/DGS/monkeypox/mpx')
mod <- cmdstan_model('SimpleSEIR.stan', pedantic=T)
fit <- mod$sample(data=data, chains=2, parallel_chains=2, iter_sampling=2000, refresh=100, iter_warmup=1000)
stanfit <- rstan::read_stan_csv(fit$output_files())

print(get_elapsed_time(stanfit))
summary_stanfit <- summary(stanfit)
summary(summary_stanfit$summary[,"Rhat"])

# check convergence
chains <- rstan::extract(stanfit)
traceplot(stanfit, pars=c('rho','phi'))

# reporting rate estimate
quantile(chains$rho, c(0.5,0.025,0.975)) 

#--- plot model fit
fit <- data.frame(t=dates, obs=ndf$cases, pred=NA, ciL=NA, ciU=NA, imports= ndf$imports)
for(t in 1:data$Nt) fit[t,3:5] <- quantile(chains$ecases[,t], c(0.5,0.025,0.975))


B <-  ggplot(df, aes(data_do_inicio_dos_sintomas, n_cases))+ 
  geom_col(aes(fill=factor(travel)))+
  geom_line(data=fit,aes(t,pred),col='plum', size=1)+ ylab('Reported cases')+ xlab('')+
  geom_ribbon(data=fit, aes(x=t, y=pred, ymin=ciL, ymax=ciU),fill='plum',alpha=0.5)+
  scale_fill_manual(values=c("snow4", "yellow3"), 
                    name="",
                    labels=c("No travel history", "Travel history"))+
  theme_classic()+
  theme(text=element_text(size=15), legend.position = c(0.8, 0.8)) 


sum(fit$pred) # total predicted cases
sum(fit$obs) # total observed cases


#--- extract Rt, Susc, Infected, Recovered estimates
rt <- data.frame(date=dates, med=NA, ciL=NA, ciU=NA)
susc <- data.frame(date=dates, med=NA, ciL=NA, ciU=NA)
ii <- data.frame(date=dates, med=NA, ciL=NA, ciU=NA)
rec <- data.frame(date=dates, med=NA, ciL=NA, ciU=NA)
for(i in 1:nrow(rt)){
  rt[i,2:4] <- quantile(chains$Rt[,i], c(0.5,0.025,0.975))
  susc[i,2:4] <- quantile(chains$S[,i], c(0.5,0.025,0.975))
  ii[i,2:4] <- quantile(chains$inc[,i], c(0.5,0.025,0.975))
  rec[i,2:4] <- quantile(chains$R[,i], c(0.5,0.025,0.975))
  
}

#--- plot Rt, Susc, Infected, Recovered estimates

rt_14
rt_28

# Rt
A <- ggplot(rt, aes(date, med))+ geom_line(col='indianred')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='indianred', alpha=0.3)+
  theme_classic()+ ylab('R[t]')+ xlab('')+ theme(text=element_text(size=15))+
  geom_hline(yintercept=1, linetype='dashed') 

# Susc prop
C <- ggplot(susc, aes(date, med))+ geom_line(col='CornFlowerblue')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='CornFlowerblue', alpha=0.3)+
  theme_classic()+ ylab('Susceptible')+ xlab('')+ theme(text=element_text(size=15)) 

# Susc N
ggplot(susc, aes(date, med*data$pop))+ geom_line(col='CornFlowerblue')+
  geom_ribbon(aes(ymin=ciL*data$pop,ymax=ciU*data$pop), fill='CornFlowerblue', alpha=0.3)+
  theme_classic()+ ylab('Susceptible')+ xlab('')+ theme(text=element_text(size=15)) 

# Infected prop
D <- ggplot(ii, aes(date, med))+ geom_line(col='plum')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='plum', alpha=0.3)+
  theme_classic()+ ylab('Infected')+ xlab('')+ theme(text=element_text(size=15))  

# Infected N
ggplot(ii, aes(date, med*data$pop))+ geom_line(col='plum')+
  geom_ribbon(aes(ymin=ciL*data$pop,ymax=ciU*data$pop), fill='plum', alpha=0.3)+
  theme_classic()+ ylab('Infected')+ xlab('')+ theme(text=element_text(size=15))  

# Recovered prop
E <- ggplot(rec, aes(date, med))+ geom_line(col='chartreuse4')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='chartreuse4', alpha=0.3)+
  theme_classic()+ ylab('Recovered')+ xlab('')+ theme(text=element_text(size=15)) 

# Recovered N
ggplot(rec, aes(date, med*data$pop))+ geom_line(col='chartreuse4')+
  geom_ribbon(aes(ymin=ciL*data$pop,ymax=ciU*data$pop), fill='chartreuse4', alpha=0.3)+
  theme_classic()+ ylab('Recovered')+ xlab('')+ theme(text=element_text(size=15))  

library(cowplot)
top_row <- plot_grid(A, B, nrow=1, labels = "AUTO")
bottom_row <- plot_grid(C, D, E, nrow=1, labels = c("C", "D", "E"))

plot_grid(top_row, bottom_row, label_size = 12, ncol = 1)

# SIR numbers
sum(ii$med*data$pop)
sum(ii$ciL*data$pop)
sum(ii$ciU*data$pop)

# Number max infected
ii %>%
filter(med == max(med)) %>%
pull(date, med)

ii[48, "med"]*data$pop
ii[48, "ciL"]*data$pop
ii[48, "ciU"]*data$pop



### Supplementary results

## Rt 14, 21 and 28 days 
rt_14
rt_28

combined_data <- rbind(
  data.frame(dataset = "21 days", rt),
  data.frame(dataset = "14 days", rt_14),
  data.frame(dataset = "28 days", rt_28)
)

ggplot(combined_data, aes(date, med)) +
  geom_line(aes (col= dataset)) +
  geom_ribbon(aes(ymin = ciL, ymax = ciU, fill = dataset), alpha = 0.3) +
  scale_fill_manual(values = c('14 days' = '#7fcdbb', '21 days' = '#2c7fb8', '28 days' = '#253494'), ) +
  scale_color_manual(values = c('14 days' = '#7fcdbb', '21 days' = '#2c7fb8', '28 days' = '#253494')) +
  labs(fill = "Infectious period", col ="") + 
  guides(col = FALSE)+
  ylab('R[t]') +
  xlab('') +
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(0.8, 0.8)) +
  geom_hline(yintercept = 1, linetype = 'dashed') 

# Different numbers of seeds
fit_2 # 2 seeds
rt_2 

fit_15
rt_15


combined_data_2 <- rbind(
  data.frame(dataset = "2 seeds", rt_2),
  data.frame(dataset = "10 seeds", rt),
  data.frame(dataset = "15 seeds", rt_15)
)

A_S <- ggplot(combined_data_2, aes(date, med)) +
  geom_line(aes (col= dataset)) +
  geom_ribbon(aes(ymin = ciL, ymax = ciU, fill = dataset), alpha = 0.5) +
  labs(fill = "Number of infection seeds", col ="") + 
  scale_fill_manual(values = c('2 seeds' = '#7fcdbb', '10 seeds' = '#2c7fb8', '15 seeds' = '#253494'),breaks=c('2 seeds', '10 seeds', '15 seeds') ) +
  scale_color_manual(values = c('2 seeds' = '#7fcdbb', '10 seeds' = '#2c7fb8', '15 seeds' = '#253494')) +
  guides(col = FALSE) +
  ylab('R[t]') +
  xlab('') +
  theme_classic()+
  theme(text = element_text(size = 15), legend.position = c(0.8, 0.8)) +
  geom_hline(yintercept = 1, linetype = 'dashed') 


combined_data_3 <- rbind(
  data.frame(dataset = "2 seeds", fit_2),
  data.frame(dataset = "10 seeds", fit),
  data.frame(dataset = "15 seeds", fit_15)
)


  
B_S <-  ggplot(data=combined_data_3,aes(t,pred))+
  geom_line(aes(col=dataset))+ 
  geom_ribbon(aes(ymin=ciL, ymax=ciU, fill = dataset),alpha=0.3)+
  scale_fill_manual(values = c('2 seeds' = '#7fcdbb', '10 seeds' = '#2c7fb8', '15 seeds' = '#253494'), ) +
  scale_color_manual(values = c('2 seeds' = '#7fcdbb', '10 seeds' = '#2c7fb8', '15 seeds' = '#253494')) +
  ylab('Reported cases')+ 
  xlab('')+
  labs(fill = "Number of infection seeds", col ="") + 
  guides(col = FALSE, fill = FALSE)+
  theme_classic()+
  theme(text=element_text(size=15), legend.position = c(0.8, 0.8)) 


plot_grid(A_S, B_S, labels = "AUTO")



