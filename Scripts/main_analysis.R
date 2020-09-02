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
cor.test(ay2$age,ay2$Year, method = "pearson")

ay3 = ay2 %>% group_by(Shrub) %>% 
  summarize(age = min(age)) %>%
  left_join(ay2)

cor.test(ay3$age, ay3$Year, method = "pearson")
plot(age~Year, data = ay3)

ay4 = ay2 %>% group_by(Year) %>% summarize(age = mean(age))

cor.test(ay4$age, ay4$Year, method = "pearson")
plot(age~Year, data = ay4)

# Calculate generic ring width statistics -----

d1 = d1 %>% 
  left_join(d2) %>%
  filter(Include == "Oui") %>%
  group_by(Year) %>%
  summarize(N = n(), sd = sd(ringwidth), chrono = mean(ringwidth))

shrubs = read.csv("Data/shrubs.csv")

# Mean sensitivity
sens1(shrubs$stdringwidth)

# Load in raw data 
rawdata<- read_csv("Data/raw_ring_data_June2019.csv")
rwdu <- read_csv("Data/final_chron_used_mod.csv")

# Prepare raw data for statistics that rely on individual cores
for_rs = rawdata %>% gather(-Year, key=ID, value=raw) %>%
  separate(ID, into =c("ID1", "ID2"), sep="-") %>%
  mutate(ID2 = str_pad(ID2, width = 3, side="right", pad= "_")) %>%
  mutate(ID2 = str_pad(ID2, width = 4, side="right", pad= "2")) %>%
  filter(!is.na(raw)) %>% 
  separate(ID2, into =c("ID2", "radius"), sep="_") %>%
  mutate(ID = paste0(ID1, "-", ID2)) %>%
  right_join(rwdu) %>%
  mutate(ID = paste0(ID1, "-", ID2, "-", radius)) %>%
  select(-ID1, -ID2, -radius) %>%
  pivot_wider(names_from = ID, values_from = raw) %>% data.frame()

rownames(for_rs) = for_rs$Year
for_rs$Year = NULL

rwi.stats(for_rs, ids = autoread.ids(for_rs))

# Prepare raw data for statistics calculated on individual shrubs
for_rwl = rawdata %>% gather(-Year, key=ID, value=raw) %>%
  separate(ID, into =c("ID1", "ID2"), sep="-") %>%
  mutate(ID2 = str_pad(ID2, width = 3, side="right", pad= "_")) %>%
  mutate(ID2 = str_pad(ID2, width = 4, side="right", pad= "2")) %>%
  filter(!is.na(raw)) %>% 
  separate(ID2, into =c("ID2", "radius"), sep="_") %>%
  mutate(ID = paste0(ID1, "-", ID2)) %>%
  right_join(rwdu) %>%
  mutate(ID = paste0(ID1, "-", ID2)) %>%
  select(-ID1, -ID2, -radius) %>%
  group_by(Year, ID) %>%
  summarize(raw = mean(raw)) %>%
  pivot_wider(names_from = ID, values_from = raw) %>% data.frame()

rownames(for_rwl) = for_rwl$Year
for_rwl$Year = NULL

rwl.report(for_rwl)

# QCQA: Verify that the two calculated chronologies are the same.
plot(shrubs$stdringwidth, d1$chrono); abline(0,1)

# Return to main analysis ----

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

dcplot(dcc(chrono = dresp, clim = CBclimate, method = "correlation", start = -6, end = 9))

respCB$ID = rownames(respCB)
rownames(respCB) = NULL

# Analysis of 1996 to 2010
respCB2 = dcc_RWB(chrono = data.frame(chrono = dresp[rownames(dresp) %in% seq(1996,2010,by=1),], row.names = seq(1996,2010,1))
                  , clim = CBclimate, method = "response", start = 3, end = 9)

respCB2$ID = rownames(respCB2)
rownames(respCB2) = NULL

respCB_plot = respCB %>% separate(ID, into = c("variable", "year", "month"), sep =c(2,8)) %>%
  mutate(Month = rep(c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se"), 2)) %>%
  select(-year, -month) %>%
  mutate(ID = "1948-2013") %>%
  bind_rows(
    respCB2 %>% separate(ID, into = c("variable", "year", "month"), sep =c(2,8)) %>%
      mutate(Month = rep(c("Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se"), 2)) %>%
      select(-year, -month) %>%
      mutate(ID = "1996-2010")
  )


p1 = respCB_plot  %>%
  filter(variable == "MT") %>%
  ggplot(aes(x = Month, y = coef, color = significant)) + 
  geom_hline(yintercept = 0, linetype = 2) + geom_point(aes(shape = ID)) + theme_classic() + geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper, linetype = ID)) + 
  scale_x_discrete(limits = c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se")) +
  scale_color_manual(values = c("grey", "black")) + ylab("Response coefficient") + theme(legend.position = "none") + ggtitle("Average temperature  (\u00B0C)") + 
  scale_shape_manual(values = c(19,1))

p2 = respCB_plot  %>%
  filter(variable == "TP") %>%
  ggplot(aes(x = Month, y = coef, color = significant)) + 
  geom_hline(yintercept = 0, linetype = 2) + geom_point(aes(shape = ID)) + theme_classic() + geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper, linetype = ID)) + 
  scale_x_discrete(limits = c("jn", "ju", "au", "se", "oc", "no", "de", "Ja", "Fe", "Mr", "Ap", "Ma", "Jn", "Ju", "Au", "Se")) +
  scale_color_manual(values = c("grey", "black")) + ylab("Response coefficient") + theme(legend.position = "none") +
  scale_shape_manual(values = c(19,1)) + ggtitle("Total precipitation (mm)")

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
            L = c("italic(R)^2==0.32","italic(p)<0.001"))

m2 = lm(MeantempAnnual~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2011))
pred.m2 = predict(m2, newdata = data.frame(Year = seq(1996, 2010, 1)))
shrubs2[,"Pred2"] = NA 
shrubs2[shrubs2$Year > 1995 & shrubs2$Year <2011, "Pred2"] = pred.m2

GT2 = tibble(X = c(1980,1980),
            Y = c(-11.5,-12),
            L = c("italic(R)^2==0","italic(p)==0.65"))

p1 = shrubs2 %>% 
  mutate(DataSet = ifelse(is.na(Pred2), "a", "b")) %>%
  ggplot(aes(x = Year, y = MeantempAnnual)) + 
  geom_point(aes(color = DataSet),shape = 1) + theme_classic() + 
  scale_x_continuous(breaks = seq(1940, 2020, by = 10)) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T, color = "blue") +
  # geom_text(aes(x = X, y = Y, label = L), data = GT2, parse = T, color = "orange") +
  geom_line(aes(y = Pred), col = "blue") + 
  # geom_line(aes(y = Pred2), col = "orange") + 
  ylab("Mean Temperature (\u00B0C)") +
  scale_color_manual(values = c("black", "orange"), guide = F)

m1 = lm(MeantempJuly~Year, data=shrubs2)
pred.m1 = predict(m1, newdata = data.frame(Year = seq(1948, 2013, 1)))
shrubs2[,"Pred"] = pred.m1

m2 = lm(MeantempJuly~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2011))
pred.m2 = predict(m2, newdata = data.frame(Year = seq(1996, 2010, 1)))
shrubs2[,"Pred2"] = NA 
shrubs2[shrubs2$Year > 1995 & shrubs2$Year <2011, "Pred2"] = pred.m2



GT = tibble(X = c(1960,1960),
            Y = c(12,11.5),
            L = c("italic(R)^2==0.15","italic(p)<0.001"))

GT2 = tibble(X = c(1980,1980),
             Y = c(12,11.5),
             L = c("italic(R)^2==0","italic(p)==0.39"))

p2 = shrubs2 %>% 
  mutate(DataSet = ifelse(is.na(Pred2), "a", "b")) %>%
  ggplot(aes(x = Year, y = MeantempJuly)) + 
  geom_point(aes(color = DataSet),shape = 1) + theme_classic() + 
  scale_x_continuous(breaks = seq(1940, 2020, by = 10)) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T, color = "blue") +
  # geom_text(aes(x = X, y = Y, label = L), data = GT2, parse = T, color = "orange") +
  geom_line(aes(y = Pred), col = "blue") + 
  # geom_line(aes(y = Pred2), col = "orange") + 
  ylab("Mean July Temperature (\u00B0C)") + 
  scale_color_manual(values = c("black", "orange"), guide = F)
  
# Check whether total precipiration has changed from 1996 to 2010:
m1_precip = lm(MeantempJuly~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2011)); summary(m1_precip)
m1_precip = lm(MeantempAnnual~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2011)); summary(m1_precip)
m1_precip = lm(TotalPrecipAnnual~Year, data=shrubs2 %>% filter(Year > 1995 & Year < 2011)); summary(m1_precip)

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

# check the response function with the ring area vector
dresparea = data.frame(chrono = area$ringarea)
rownames(dresparea) = area$Year

respCBarea = dcc_RWB(chrono = dresparea, clim = CBclimate, method = "response", start = -6, end = 9)

rm(dresparea)
dcplot(respCBarea) # yes the July temperature variable is still significant, now so is the past June, but not confident in that result

# Load in the raw data 
area2 = read_csv("Data/ringarea2.csv") %>% gather(-Year, key = ID, value = ringarea) %>% filter(!is.na(ringarea)) %>% mutate(ringarea = ringarea/1e6) %>% group_by(Year) %>%
  summarize(sd = sd(ringarea), ringarea = mean(ringarea), N = n()) 

m1 = lm(ringarea~Year, data = area2); summary(m1)
pred.m1 = predict(m1, newdata = data.frame(Year = seq(1922, 2014, 1)))
area2[,"Pred"] = pred.m1

GT = tibble(X = c(1930,1930),
       Y = c(2,1.8),
       L = c("italic(R)^2==0.74","italic(p)<0.001"))


m2 = lm(ringarea~Year, data = area2 %>% filter(Year > 1995 & Year < 2011)); summary(m2)
pred.m2 = predict(m2, newdata = data.frame(Year = seq(1996, 2010, 1)))
area2[,"Pred2"] = NA
area2[area2$Year > 1995 & area2$Year < 2011,"Pred2"] = pred.m2

GT2 = tibble(X = c(1950,1950),
            Y = c(2,1.8),
            L = c("italic(R)^2==0","italic(p)==0.74"))

p1 = area2 %>%
  mutate(upper = ringarea + sd, lower = ringarea -sd) %>%
  mutate(N = N/100) %>%
  mutate(lower = ifelse(lower < 0, 0, lower)) %>%
  ggplot(aes(x = Year, y = ringarea)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey") + 
  theme_classic() + 
  geom_line(lwd = 1.1) + 
  geom_line(data = area2 %>%
              filter(!is.na(Pred2)), color = "orange", lwd = 1.2) + 
  geom_line(aes(y = Pred), col = "blue") +
  # geom_line(aes(y = Pred2), col = "orange") +
  scale_x_continuous(breaks = seq(1920, 2020, by = 10)) + ylab(expression(Ring~area~(mm^2))) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T, col = "blue")
  # geom_text(aes(x = X, y = Y, label = L), data = GT2, parse = T, col = "orange")

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

m2 = lm(ringarea~JulyTemp, data = area3 %>% filter(Year > 1995 & Year < 2011)); summary(m2)
pred.m2 = predict(m2, newdata = data.frame(JulyTemp = area3[,"JulyTemp"]))
area3[,"Pred2"] = pred.m2
area3[area3$JulyTemp < min(area3[area3$Year > 1995 & area3$Year < 2011, "JulyTemp"]),"Pred2"] = NA

GT2 = tibble(X = c(6.5,6.5),
            Y = c(2,1.8),
            L = c("italic(R)^2==0.22","italic(p)<0.04"))

p2 = area3 %>%
  mutate(DataSet = ifelse(Year > 1995 & Year < 2011, "a", "b")) %>%
  mutate(upper = ringarea + sd, lower = ringarea -sd) %>%
  mutate(lower = ifelse(lower < 0, 0, lower)) %>%
  ggplot(aes(x = JulyTemp, y = ringarea)) + 
  theme_classic() + geom_pointrange(aes(ymin = lower, ymax = upper, color = DataSet),shape = 1) + 
  geom_line(aes(x = JulyTemp, y = Pred), col = "blue") + 
  geom_line(aes(x = JulyTemp, y = Pred2), col = "orange") + 
  ylab(expression(Ring~area~(mm^2))) +
  geom_text(aes(x = X, y = Y, label = L), data = GT, parse = T, col="blue") + 
  geom_text(aes(x = X, y = Y, label = L), data = GT2, parse = T, col = "orange") + 
  xlab("Mean July Temperature (\u00B0C)") +
  scale_color_manual(values = c("orange", "black"), guide = F)

png("Plots/Figure4.png", width = 5, height = 7, units = "in", res = 600)
ggpubr::ggarrange(p1,p2, labels = "AUTO", ncol = 1, nrow = 2)
dev.off()  

# Test shrub growth over the same period as shrub cover (1996-2010)
m1_mod = lm(ringarea~Year, data = area2 %>% filter(Year > 1995 & Year < 2011)); summary(m1_mod) # Not significant!
plot(m1_mod)

# Verify with LMM
area4 = read_csv("Data/ringarea2.csv") %>% gather(-Year, key = ID, value = ringarea) %>% filter(!is.na(ringarea)) %>% mutate(ringarea = ringarea/1e6)

lmm1 <- nlme::lme(ringarea~Year, random=~1|ID, data = area4)
summary(lmm1) # Looks good!

lmm2 <- nlme::lme(ringarea~Year, random=~1|ID, data = area4 %>% filter(Year > 1995 & Year < 2011))
summary(lmm2) # Same small effect

# Run a loop testing how many shrubs would be needed:
listshrubs = unique(area4$ID)
Nshrubs = length(listshrubs) 

out1 = array(NA, dim = c(53,2,100))

for(j in 1:100){
  for(i in 1:53){
    
    SS = base::sample(listshrubs, i)
    
    lmm1 <- nlme::lme(ringarea~Year, random=~1|ID, data = area4 %>% filter(ID %in% SS))
    
    out1[i,1,j] = summary(lmm1)$coefficients$fixed[2]
    out1[i,2,j] = summary(lmm1)$tTable[2,5]
  }
}

out2 = apply(out1, c(1,2), FUN = mean)
out2sd = apply(out1, c(1,2), FUN = sd)

plot(out2[,1]~seq(1:53))
plot(out2sd[,1]~seq(1:53))
plot(out2[,2]~seq(1:53))
# Conclusion: need over 10 shrubs from 1 site to get a robust relationship over time.

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

# Look at the ring series statistics for the area series
area5 = area4 %>% spread(key = ID, value = ringarea) %>% as.data.frame()

rownames(area5) = area5$Year
area5$Year = NULL

rwi.stats(area5)

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
           color = 'black', size = 3) +
  annotate(geom = "text", x = c(-110, -102, -122, -110), 
           y =c(71.3,67.5, 73, 65.5), 
           label = c("Victoria Island", "Kent Peninsula", "Banks \nIsland", "Mainland Nunavut"), 
           color = 'black', size = 4) + 
  annotate(geom = "segment", x = -106, y = 67.8, xend = -107, yend = 68.6, arrow = arrow(length = unit(0.10,"cm")))
dev.off()

# Moving analysis 14-year averages ----

# .....July temperature -----

# Years to consider
YEARS = unique(CBclimate$year)
YEARS2 = data.frame(start = YEARS,
                    end = YEARS + 14) %>% filter(end <= max(YEARS) & start >= 1949)
YEARS2$`Ring width_value` = NA
YEARS2$`Ring width_sig` = NA
YEARS2$`Ring area_value` = NA
YEARS2$`Ring area_sig` = NA
selchrono = as.numeric(rownames(dresp))


for(i in 1:dim(YEARS2)[1]){
  CBclimatecur = CBclimate %>% filter(year >= YEARS2$start[i] & year <= YEARS2$end[i] & month == 7)
  
  drespcur = dresp[selchrono >= YEARS2$start[i] & selchrono <= YEARS2$end[i],]
  
  areacur = area3 %>% filter(Year >= YEARS2$start[i] & Year <= YEARS2$end[i]) %>% as.data.frame()
  
  m1 = lm(drespcur~CBclimatecur$MT)
  YEARS2$`Ring width_value`[i] = coefficients(m1)[2]
  YEARS2$`Ring width_sig`[i] = coefficients(summary(m1))[2,4] < 0.05
  
  m1 = lm(areacur$ringarea~CBclimatecur$MT)
  YEARS2$`Ring area_value`[i] = coefficients(m1)[2]
  YEARS2$`Ring area_sig`[i] = coefficients(summary(m1))[2,4] < 0.05
}

YEARS2 %>% write_csv("DataOut/movingcorr_July.csv")

p1 = read_csv("DataOut/movingcorr_July.csv") %>% select(-end) %>% gather(-start,key = `Data Type`, value = value) %>% 
  separate(`Data Type`, into = c('Data Type', "temp"), sep = "_") %>%
  pivot_wider(names_from = temp, values_from = value) %>%
  mutate(sig = ifelse(sig == 0, "No", "Yes")) %>%
  ggplot(aes(x = start, y = value)) + geom_line(aes(linetype = `Data Type`), color = "grey") + theme_classic() + geom_point(aes(color = sig, shape = `Data Type`)) + ylab("Coefficient (Growth~July Temperature)") + xlab("Start of 14-yr window") + xlim(range(read_csv("DataOut/outmove.csv")$start)) +
  theme(legend.position = c(0.15, 0.7)) + geom_hline(yintercept = 0, lty = 2) +
  scale_color_manual(values = c("#E69F00", "purple4"), name = "Significant?", guide = F) +
  scale_shape_manual(values = c(19,1))  +
  annotate(geom= "segment",
           x = 1996, xend = 1996,
           y = 0.045, yend = 0.054,
           arrow = arrow(length = unit(0.1, "cm")))

# ......Test over all 14 year periods in the data set -----
outmove = data.frame(start = 1922:2000,
                     end = 1936:2014, 
                     coeff = NA,
                     p = NA)

for(i in 1:dim(outmove)[1]){
  m1_mod = lm(ringarea~Year, data = area2 %>% filter(Year >= outmove$start[i] & Year <= outmove$end[i]))
  outmove$coeff[i] = summary(m1_mod)$coefficients[2,1]
  outmove$p[i] = summary(m1_mod)$coefficients[2,4]
}

outmove$sig = ifelse(outmove$p < 0.05, "Yes", "No")
outmove$focus = ifelse(outmove$start ==1996, "Yes", "No")

write_csv(outmove, "DataOut/outmove.csv")

p2 = read_csv("DataOut/outmove.csv") %>% ggplot(aes(x = start, y = coeff))+ geom_hline(yintercept = 0, lty = 2) + geom_line(color = "grey") + geom_point(aes(col = sig)) + theme_classic() + xlab("") + ylab("Coefficient (Growth~Year)") +
  scale_color_manual(values = c("#E69F00", "purple4"), name = "Significant?") +
  theme(legend.position = c(0.15, 0.2)) +
  annotate(geom= "segment",
           x = 1996, xend = 1996,
           y = -0.01, yend = -0.006,
           arrow = arrow(length = unit(0.1, "cm")))


png("Plots/Figure6.png", width = 5, height = 7, units = "in", res = 600)
ggpubr::ggarrange(p2,p1, labels = "AUTO", ncol = 1, nrow = 2)
dev.off()

# Correlate chronology with 1996 to 2010 significant response function ----

selll = which(rownames(dresp) %in% c("1996", "2010"))

selll = dresp[seq(selll[1], selll[2],1),]

selll = CBclimate %>% 
  filter(year > 1995 & year < 2011) %>%
  select(-MT) %>%
  filter(month %in% c(5,7,8)) %>%
  pivot_wider(names_from = month, values_from = TP) %>%
  mutate(chron = selll) %>%
  rename(May = `5`,
         July = `7`,
         August = `8`)

m1 <- lm(chron~May, data = selll); summary(m1)
m1 <- lm(chron~July, data = selll); summary(m1)
m1 <- lm(chron~August, data = selll); summary(m1)
plot(chron~May, data= selll)
plot(chron~July, data= selll)
plot(chron~August, data= selll)

selll %>% write_csv("Data/TP_1996to2010.csv")

read_csv("Data/TP_1996to2010.csv") %>% pivot_longer(May:August) %>%
  rename(Year = year) %>%
  ggplot(aes(x = value, y = chron)) + 
  geom_point(aes(color = Year)) + facet_wrap(.~name, ncol = 1, nrow = 3) + 
  theme_classic() + 
  xlab("Total Precipitation (mm)") + 
  ylab("RWI") +
  stat_smooth(method = "lm") +
  scale_color_viridis_c()
