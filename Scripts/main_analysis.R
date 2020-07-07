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

# Check whether the July temperature has increased from 1996 to 2010:
m1_short = lm(MeantempJuly~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2010))
summary(m1_short)
plot(m1_short)

# Check whether total precipiration has changed from 1996 to 2010:
m1_precip = lm(TotalPrecipAnnual~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2010))
summary(m1_precip)
plot(m1_precip)

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


# Breakpoint analysis:
m1_breakpoint = segmented::segmented(m1)
summary(m1_breakpoint)
plot(m1_breakpoint)

out = data.frame(SY = 1924:2012)
out$BP = NA

for(i in 2:dim(out)[1]){
  out$BP[i] = segmented::segmented(m1, psi = out$SY[i])$psi[2]
}
plot(BP~SY, data = out)
# No robust breakpoint in these data: could be ~2007 or ~1948, but this is not capturing large sections of the data or dependent systematically on where the analysis starts
rm(out, m1_breakpoint)

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

# Test shrub growth over the same period as shrub cover (1996-2010)
m1_mod = lm(ringarea~Year, data = area2 %>% filter(Year > 1995 & Year < 2011)); summary(m1_mod) # Not significant!
plot(m1_mod)

# Test over all 14 year periods in the data set

outmove = data.frame(start = 1922:2000,
           end = 1936:2014, 
           coeff = NA,
           p = NA)

for(i in 1:dim(outmove)[1]){
  m1_mod = lm(ringarea~Year, data = area2 %>% filter(Year >= outmove$start[i] & Year <= outmove$end[i]))
  outmove$coeff[i] = summary(m1_mod)$coefficients[2,1]
  outmove$p[i] = summary(m1_mod)$coefficients[2,4]
}

outmove$sig = ifelse(outmove$p < 0.05, "blue", "orange")

plot(coeff~start, outmove, col = sig, type = "b", ylab = "Year Coefficient", xlab = "Start of 14-yr window")
abline(v = 1996, lty = 2)
legend("bottomleft", legend = c("p < 0.05", "p > 0.05"), col = c("blue", "orange"), pch = 1)
# Verify with LMM

area4 = read_csv("Data/ringarea2.csv") %>% gather(-Year, key = ID, value = ringarea) %>% filter(!is.na(ringarea)) %>% mutate(ringarea = ringarea/1e6)

lmm1 <- nlme::lme(ringarea~Year, random=~1|ID, data = area4)
summary(lmm1) # Looks good!

lmm2 <- nlme::lme(ringarea~Year, random=~1|ID, data = area4 %>% filter(Year > 1995 & Year < 2011))
summary(lmm2) # Same small effect

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

# Figure 0: The map of our study location

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggspatial)
library(rgeos)

world <- ne_countries(scale = "medium", returnclass = "sf")

png("Plots/Figure0.png", width = 5, height = 5, units = "in",res = 600)
ggplot() +
  geom_sf(data = world) + xlab("Longtitude") + ylab("Latitude") + 
  coord_sf(xlim = c(-125, -95), ylim = c(65, 75), expand = FALSE) +
  theme_minimal()+
  annotation_north_arrow(location = "bl", which_north = "true", 
                         # pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(aes(x = c(-105.053056,-115.095,-108.083333,-108.416667),
                 y = c(69.117222,67.825556,68.350000,68.916667)), 
             pch = 19, size = 2) + 
  annotate(geom = "text", x = c(-105.053056-0.5,-115.095-3,-108.083333-3.2,-108.416667-2.5), 
           y =c(69.117222+0.75,67.825556,68.350000-0.1,68.916667+0.4), 
           label = c("Cambridge \n Bay", "Kugluktuk", "Walker Bay", "Byron Bay"), 
           color = 'black', size = 3)
dev.off()

# Moving analysis 11-year averages for the supplemental ----

# Years to consider
YEARS = unique(CBclimate$year)
YEARS2 = data.frame(start = YEARS,
                    end = YEARS + 11,
                    middle = YEARS + 11/2) %>% filter(end <= max(YEARS) & start >= 1949)
YEARS2$`Ring width` = NA
YEARS2$`Ring area` = NA
selchrono = as.numeric(rownames(dresp))


for(i in 1:dim(YEARS2)[1]){
  CBclimatecur = CBclimate %>% filter(year >= YEARS2$start[i] & year <= YEARS2$end[i] & month == 7)
  
  drespcur = dresp[selchrono >= YEARS2$start[i] & selchrono <= YEARS2$end[i],]
  
  areacur = area3 %>% filter(Year >= YEARS2$start[i] & Year <= YEARS2$end[i]) %>% as.data.frame()
  
  YEARS2$`Ring width`[i] = cor(CBclimatecur$MT, drespcur, method = 'pearson')
  YEARS2$`Ring area`[i] = cor(CBclimatecur$MT, areacur$ringarea, method = 'pearson')
}

YEARS2 %>% write_csv("DataOut/movingcorr.csv")

YEARS2 %>% select(-end, -start) %>% gather(-middle,key = `Data Type`, value = value) %>% ggplot(aes(x = middle, y = value, linetype = `Data Type`)) + geom_line() + theme_classic() + ylab("Correlation: July temperature and growth") + xlab("Center of 11-yr window")

YEARS2 %>% filter(middle < 1970) %>% summarize(mean(`Ring area`), mean(`Ring width`))

YEARS2 %>% filter(middle > 1980) %>% summarize(mean(`Ring area`), mean(`Ring width`))