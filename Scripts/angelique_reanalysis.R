###############################################
library(bootRes)
library(treeclim)
require(dplR)

#to have the standardized and residual chronologies
data1 <- read.csv("Data/responsefunctionanalysis/raw_ringwidth.csv",sep=";")
row.names(data1) = data1$Year

#to do the crossdating (note that it took me several steps to conclude that I needed to remove the last 13 shrubs in the database)
sercor=corr.rwl.seg(data1[,-c(1,59:71)],seg.length=10,bin.floor = 8) # I removed the last 13 shrubs to have an eps > 0.85
rwl.report(data1[,-c(1,59:71)]) #autocorrelation1 = 0.42 and interseries correlation = 0.41

#to create the new database with an EPS > 0.85
data1 = subset(data1[,-c(1,59:71)], data1$Year >= 1949 & data1$Year  <= 2013) # we have data from only 2 shrubs in 1949 but three and more after, so I started the analysis in 1949 ond not in 1948 as indicated by Stephane (and stop in 2013 as Clara did)
datadetrend <- detrend(data1, method = "Spline",nyrs=9) # cubic spline with degree 9, as explained by Clara.
rwi.stats(datadetrend) # EPS>0.85

# to create the mean chronology
chronostd=chron(datadetrend, prefix = "CAM") # based on 57 shrubs

# Load the climate data
CBclimate = read.csv("Data/responsefunctionanalysis/climate_data.csv",sep=";")
CBclimate2=subset(CBclimate, CBclimate$Year >= 1948 & CBclimate$Year  <= 2013)

# to do the response function
respCB1=dcc(chronostd, CBclimate2, method = "response",selection = -6:9,boot = "stationary",ci = 0.05)
# note that the analysis is performed on the standardized chronology, not on the residual one as specified in the MS.
respCB1

# to do the response function: ROB
CBclimate = read.csv("Data/en_climate_monthly_NU_2400600_1929-2015_P1M.csv") %>%
  select(Year, Month, `Total.Precip..mm.`, `Mean.Temp...C.`) %>%
  rename(TP = `Total.Precip..mm.`, MT = `Mean.Temp...C.`,
         year = Year, month = Month) %>%
  filter(year > 1947 & year < 2014)

respCB1=dcc_RWB(chronostd, CBclimate, method = "response", start = -6, end = 9)
# note that the analysis is performed on the standardized chronology, not on the residual one as specified in the MS.
respCB1



# Raw reanalysis of trends over time

# rawdata<- read_csv("Data/raw_ring_data_June2019.csv")

# rwdu <- read_csv("Data/final_chron_used_mod.csv")

# rawsum <- rawdata %>% gather(-Year, key=ID, value=raw) %>%
#   separate(ID, into =c("ID1", "ID2"), sep="-") %>%
#   mutate(ID2 = str_pad(ID2, width = 3, side="right", pad= "_")) %>%
#   mutate(ID2 = str_pad(ID2, width = 4, side="right", pad= "2")) %>%
#   filter(!is.na(raw)) %>%
#   separate(ID2, into =c("ID2", "radius"), sep="_") %>%
#   mutate(ID = paste0(ID1, "-", ID2)) %>%
#   select(-ID1, -ID2) %>%
#   group_by(Year, ID) %>%
#   summarize(raw = mean(raw))

# data1$Year = rownames(data1)
# rownames(data1) = NULL

rawdata = cbind(data.frame(Year = as.numeric(rownames(data1))), data1)
rownames(rawdata) = NULL

stddata = finaldata = rawdata

number = dim(rawdata)[2]

selvec = rep(NA, number)
slopevec = rep(NA, number)
sqrvec = rep(NA, number)
sqrvec2 = rep(NA, number)

for(i in 2:number){
  rawdatacur = rawdata[,c(1,i)]
  cur =rawdatacur[!is.na(rawdatacur[,2]),] 
  colnames(cur) = c("Year", "Chron")
  cur["age"] = seq(1, dim(cur)[1], 1)
  cur["age2"] = cur["age"]^2
  m1 = lm(Chron~0, data=cur); summary(m1)
  m2 = lm(Chron~age, data=cur); summary(m2)
  m3 = lm(Chron~age2 + age, data=cur); summary(m3)
  sel = which(AIC(m1, m2, m3)$AIC == min(AIC(m1, m2, m3)$AIC))
  modelsel = rownames(AIC(m1, m2, m3))[sel]
  selvec[i] = modelsel
  slopevec[i] = coefficients(m2)[2]
  sqrvec[i] = coefficients(m3)[2]
  sqrvec2[i] = coefficients(m3)[3]
  resids = resid(get(modelsel))
  finaldata[as.numeric(names(resids)),i] = resids
}

for(i in 2:number){
  rawdatacur = rawdata[,c(1,i)]
  cur =rawdatacur[!is.na(rawdatacur[,2]),] 
  colnames(cur) = c("Year", "Chron")
  cur["age"] = seq(1, dim(cur)[1], 1)
  cur["age2"] = cur["age"]^2
  m1 = lm(Chron~0, data=cur); summary(m1)
  resids = resid(m1)
  stddata[as.numeric(names(resids)),i] = resids
}

# Get positive negative trends
ddd = data.frame(ID = colnames(rawdata)[-1], Slope = slopevec[-1])
ddd[,"PN"] = ifelse(ddd$Slope > 0, "+", "-")
ddd[,"Year"] = 1930
ddd[,"raw"] = 600

finaldataMEAN = finaldata %>% as_tibble() %>% 
  mutate(Year = rownames(data1)) %>%
  mutate(Year = as.numeric(Year)) %>%
  gather(-Year, key=ID, value=ringres) %>%
  filter(!is.na(ringres)) %>%
  group_by(Year) %>%
  summarize(mring = mean(ringres))

stddataMEAN = stddata %>% as_tibble() %>% 
  mutate(Year = rownames(data1)) %>%
  mutate(Year = as.numeric(Year)) %>%
  gather(-Year, key=ID, value=ringres) %>%
  filter(!is.na(ringres)) %>%
  group_by(Year) %>%
  summarize(mring = mean(ringres))

stdataMEAN2 = stddata %>%
  mutate(Year = as.character(Year)) %>%
  gather(-Year, key=ID, value=ringres) %>%
  filter(!is.na(ringres)) %>%
  group_by(ID) %>%
  summarize(N = n()) %>%
  filter(N > 10) %>% select(-N) %>% 
  left_join(
    stddata %>%
      mutate(Year = as.numeric(Year)) %>%
      gather(-Year, key=ID, value=ringres) %>%
      filter(!is.na(ringres))
  ) %>%
  group_by(Year) %>%
  summarize(mring = mean(ringres))

m1 <- lm(mring~Year, data = finaldataMEAN); summary(m1)
m2 <- lm(mring~Year, data = stddataMEAN); summary(m2)
m3 <- lm(mring~Year, data = stdataMEAN2); summary(m3)


# Subset shrubs
stddata %>%
  mutate(Year = as.character(Year)) %>%
  gather(-Year, key=ID, value=ringres) %>%
  filter(!is.na(ringres)) %>%
  group_by(ID) %>%
  summarize(N = n()) %>% ggplot(aes(x = N)) + geom_histogram()


output = matrix(NA, nrow = 10, ncol = 3)

for(i in 1:10){
  tc = seq(-1,-i,-1)
  
  temp = stddata %>%
    mutate(Year = as.character(Year)) %>%
    gather(-Year, key=ID, value=ringres) %>%
    filter(!is.na(ringres)) %>%
    group_by(ID) %>%
    summarize(N = n()) %>%
    arrange(N) %>%
    slice(tc) %>% select(-N) %>% 
    left_join(
      stddata %>%
        mutate(Year = as.numeric(Year)) %>%
        gather(-Year, key=ID, value=ringres) %>%
        filter(!is.na(ringres))
    ) %>%
    group_by(Year) %>%
    summarize(mring = mean(ringres))
  
  m4 <- lm(mring~Year, data = temp)
  output[i,] = summary(m4)$coefficients[2,c(1,2,4)]
}

output = data.frame(rbind(summary(m2)$coefficients[2,c(1,2,4)], output))
output$RM = seq(0,10,1)
colnames(output) = c("Coef", "StdEr","p", "RM")
output$Sig = ifelse(output$p < 0.05, "Yes","No")
output$lower = output$Coef - output$StdEr
output$upper = output$Coef + output$StdEr

# Create plots

p1 = finaldataMEAN %>%
  ggplot(aes(x=Year, y=mring)) +
  geom_line() +
  theme_classic() + ylab("Chronology")

p2 = stddataMEAN %>%
  ggplot(aes(x=Year, y=mring)) +
  geom_line() +
  stat_smooth(method="lm") + theme_classic() + ylab("Chronology") +
  geom_text(x=1960, y=120, label=paste0("p < ",signif(summary(m2)$coefficients[2,4],2)), size=3)

p3 = output %>% ggplot(aes(x = RM, y = Coef, color = Sig)) + geom_point() + theme_classic() + geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.1)) + 
  xlab("Number of shrubs removed") + ylab("Year Coefficient") + scale_color_discrete(name = "Significant")

ggpubr::ggarrange(p2,p3, labels = "AUTO")
