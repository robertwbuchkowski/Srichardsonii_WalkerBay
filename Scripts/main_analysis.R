# Analysis of the Data

library(tidyverse)
library(bootRes)
require(dplR)
source("Scripts/custom_dendro_functions.R")

# Load in the data ----

d1 = read.csv("Data/shrubs_full.csv") %>%
  gather(-Year, key = Shrub, value = ringwidth) %>%
  filter(!is.na(ringwidth))

d2 = read.csv("Data/list_to_include.csv")

all(d1$Shrub %in% d2$Shrub)

# Compare age and year ----

ay1 = read.csv("Data/shrubs_full.csv")
ay2 = ay1

for(i in 2:71){
  N = length(ay2[,i][!is.na(ay2[,i])])
  ay2[,i][!is.na(ay2[,i])] = seq(1, N, 1)
}

ay2 = ay2 %>%
  gather(-Year, key = Shrub, value = age) %>%
  filter(!is.na(age)) %>%
  left_join(d2) %>%
  filter(Include == "Oui")

plot(age~Year, data = ay2)

cor(ay2$age,ay2$Year, method = "pearson")

# Return to analysis ----

d1 = d1 %>% 
  left_join(d2) %>%
  filter(Include == "Oui") %>%
  group_by(Year) %>%
  summarize(N = n(), sd = sd(ringwidth), chrono = mean(ringwidth))

shrubs = read.csv("Data/shrubs.csv")

# QCQA: Verify that the two calculated chronologies are the same.
plot(shrubs$stdringwidth, d1$chrono); abline(0,1)

# Load the climate data
CBclimate = read.csv("Data/en_climate_monthly_NU_2400600_1929-2015_P1M.csv") %>%
  select(Year, Month, `Total.Precip..mm.`, `Mean.Temp...C.`) %>%
  rename(TP = `Total.Precip..mm.`, MT = `Mean.Temp...C.`,
         year = Year, month = Month) %>%
  filter(year > 1947 & year < 2014)


# Response function analysis ----

dresp = data.frame(d1)
rownames(dresp) = dresp$Year
dresp$Year = NULL
dresp$sd = NULL
dresp$N = NULL

# This is the custom dcc function with higher replication
respCB = dcc_RWB(chrono = dresp, clim = CBclimate, method = "response", start = -6, end = 9)

respCB$ID = rownames(respCB)
rownames(respCB) = NULL

p1 = respCB %>% separate(ID, into = c("variable", "year", "month"), sep =c(2,8)) %>%
  mutate(Month = rep(c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se"), 2)) %>%
  select(-year, -month) %>%
  filter(variable == "MT") %>%
  ggplot(aes(x = Month, y = coef, color = significant)) + 
  geom_hline(yintercept = 0, linetype = 2) + geom_point() + theme_classic() + geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper)) + 
  scale_x_discrete(limits = c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se")) +
  scale_color_manual(values = c("grey", "black")) + ylab("Response coefficient") + theme(legend.position = "none") + ggtitle("Average temperature  (\u00B0C)")

p2 = respCB %>% separate(ID, into = c("variable", "year", "month"), sep =c(2,8)) %>%
  mutate(Month = rep(c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se"), 2)) %>%
  select(-year, -month) %>%
  filter(variable == "TP") %>%
  ggplot(aes(x = Month, y = coef, color = significant)) + 
  geom_hline(yintercept = 0, linetype = 2) + geom_point() + theme_classic() + geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper)) + 
  scale_x_discrete(limits = c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se")) +
  scale_color_manual(values = c("grey", "black")) + ylab("Response coefficient") + theme(legend.position = "none") + ggtitle("Total precipitation (mm)")

png("Plots/Figure3.png", width = 5, height = 7, units = "in", res = 600)
ggpubr::ggarrange(p1,p2, labels = "AUTO", ncol = 1, nrow = 2)
dev.off()

# Figure 1 ----
shrubs2 = subset(shrubs, Year > 1948)
m1 = lm(MeantempAnnual~Year, data=shrubs2)
pred.m1 = predict(m1, newdata = data.frame(Year = seq(1948, 2013, 1)))
shrubs2[,"Pred"] = pred.m1

GT = tibble(X = c(1960,1960),
            Y = c(-11.5,-12),
            L = c("italic(R)^2==0.32","italic(p)<0.0001"))

p1 = shrubs2 %>% 
  ggplot(aes(x = Year, y = MeantempAnnual)) + 
  geom_point(shape = 1) + theme_classic() + 
  scale_x_continuous(breaks = seq(1940, 2020, by = 10)) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T) +
  geom_line(aes(y = Pred), col = "blue") + ylab("Mean Temperature (\u00B0C)")

m1 = lm(MeantempJuly~Year, data=shrubs2)
pred.m1 = predict(m1, newdata = data.frame(Year = seq(1948, 2013, 1)))
shrubs2[,"Pred"] = pred.m1

GT = tibble(X = c(1960,1960),
            Y = c(12,11.5),
            L = c("italic(R)^2==0.15","italic(p)<0.001"))

p2 = shrubs2 %>% 
  ggplot(aes(x = Year, y = MeantempJuly)) + 
  geom_point(shape = 1) + theme_classic() + 
  scale_x_continuous(breaks = seq(1940, 2020, by = 10)) + 
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T) +
  geom_line(aes(y = Pred), col = "blue") + ylab("Mean July Temperature (\u00B0C)")

png("Plots/Figure1.png", width = 5, height = 7, units = "in", res = 600)
ggpubr::ggarrange(p1,p2, labels = "AUTO", ncol = 1, nrow = 2)
dev.off()     

png("Plots/Figure2.png", width = 6, height = 5, units = "in", res = 600)
d1 %>% mutate(upper = chrono + sd, lower = chrono -sd) %>%
  mutate(N = N/100) %>%
  ggplot(aes(x = Year, y = chrono)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey") + 
  theme_classic() + geom_line(lwd = 1.1) + ylab("Standardized Ring Width") + scale_y_continuous(sec.axis = sec_axis(~ . *100, name = "Sample depth")) +
  geom_ribbon(aes(x = Year, ymin = 0, ymax = N), fill = "snow4")
dev.off()

# Ring area analysis ----

# load the the data summarized in Excel
area = read_csv("Data/ringarea.csv")

# Verify the same July temperature vector
area$JulyTemp[!is.na(area$JulyTemp)]-shrubs2$MeantempJuly

p1 = area %>% 
  ggplot(aes(x = Year, y = ringarea)) + 
  geom_line() + theme_classic() + 
  scale_x_continuous(breaks = seq(1920, 2020, by = 10)) + ylab(expression(Ring~area~(mm^2)))

# Load in the raw data 
area2 = read_csv("Data/ringarea2.csv") %>% gather(-Year, key = ID, value = ringarea) %>% filter(!is.na(ringarea)) %>% mutate(ringarea = ringarea/1e6) %>% group_by(Year) %>%
  summarize(sd = sd(ringarea), ringarea = mean(ringarea), N = n()) 

m1 = lm(ringarea~Year, data = area2); summary(m1)
pred.m1 = predict(m1, newdata = data.frame(Year = seq(1922, 2014, 1)))
area2[,"Pred"] = pred.m1

GT = tibble(X = c(1930,1930),
       Y = c(2,1.8),
       L = c("italic(R)^2==0.74","italic(p)<0.001"))

p1 = area2 %>%
  mutate(upper = ringarea + sd, lower = ringarea -sd) %>%
  mutate(N = N/100) %>%
  mutate(lower = ifelse(lower < 0, 0, lower)) %>%
  ggplot(aes(x = Year, y = ringarea)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey") + 
  theme_classic() + geom_line(lwd = 1.1) + 
  geom_line(aes(y = Pred), col = "blue") + 
  scale_x_continuous(breaks = seq(1920, 2020, by = 10)) + ylab(expression(Ring~area~(mm^2))) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T)

area3 = area %>% select(Year, JulyTemp) %>%
  filter(!is.na(JulyTemp)) %>%
  left_join(
    area2 %>% select(Year, ringarea,sd)
  )

m1 = lm(ringarea~JulyTemp, data = area3); summary(m1)
pred.m1 = predict(m1, newdata = data.frame(JulyTemp = area3$JulyTemp))
area3[,"Pred"] = pred.m1

GT = tibble(X = c(5,5),
            Y = c(2,1.8),
            L = c("italic(R)^2==0.27","italic(p)<0.001"))

p2 = area3 %>%
  mutate(upper = ringarea + sd, lower = ringarea -sd) %>%
  mutate(lower = ifelse(lower < 0, 0, lower)) %>%
  ggplot(aes(x = JulyTemp, y = ringarea)) + 
  theme_classic() + geom_pointrange(aes(ymin = lower, ymax = upper),shape = 1) + 
  geom_line(aes(y = Pred), col = "blue") + 
  ylab(expression(Ring~area~(mm^2))) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T) + 
  xlab("Mean July Temperature (\u00B0C)") 

png("Plots/Figure4.png", width = 5, height = 7, units = "in", res = 600)
ggpubr::ggarrange(p1,p2, labels = "AUTO", ncol = 1, nrow = 2)
dev.off()  

# Verify with LMM

area4 = read_csv("Data/ringarea2.csv") %>% gather(-Year, key = ID, value = ringarea) %>% filter(!is.na(ringarea)) %>% mutate(ringarea = ringarea/1e6)

lmm1 <- nlme::lme(ringarea~Year, random=~1|ID, data = area4)
summary(lmm1) # Looks good!

png("Plots/INDringarea.png", width = 7, height = 7, units = "in", res = 600)
read_csv("Data/ringarea2.csv") %>% gather(-Year, key = ID, value = ringarea) %>% filter(!is.na(ringarea)) %>% mutate(ringarea = ringarea/1e6) %>%
  ggplot(aes(x=Year, y=ringarea)) +
  geom_line() + 
  facet_wrap(.~ID) + 
  ylab(expression(Ring~area~(mm^2))) +
  theme_classic() + 
  scale_x_continuous(breaks = c(1920, 1960, 2000),
                                       minor_breaks = NULL)
dev.off()
