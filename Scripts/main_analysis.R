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

# Compare age and year

ay1 = read.csv("Data/shrubs.csv")
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

# Return to analysis

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

respCB = dcc_RWB(chrono = dresp, clim = CBclimate, method = "response", start = -6, end = 9)

respCB$ID = rownames(respCB)
rownames(respCB) = NULL

mutate(variable = ifelse(variable == "MT", "Average temperature  (\u00B0C)", "Total precipitation (mm)"))

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

p1 = shrubs2 %>% 
  ggplot(aes(x = Year, y = MeantempAnnual)) + 
  geom_point(shape = 1, size = 3) + theme_classic() + 
  scale_x_continuous(breaks = seq(1940, 2020, by = 10)) + 
  geom_line(aes(y = Pred)) + ylab("Mean Temperature (\u00B0C)")

m1 = lm(MeantempJuly~Year, data=shrubs2)
pred.m1 = predict(m1, newdata = data.frame(Year = seq(1948, 2013, 1)))
shrubs2[,"Pred"] = pred.m1

p2 = shrubs2 %>% 
  ggplot(aes(x = Year, y = MeantempJuly)) + 
  geom_point(shape = 1, size = 3) + theme_classic() + 
  scale_x_continuous(breaks = seq(1940, 2020, by = 10)) + 
  geom_line(aes(y = Pred)) + ylab("Mean July Temperature (\u00B0C)")

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

