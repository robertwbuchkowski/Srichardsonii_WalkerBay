CBclimate %>%
  rename(Year = year, Month = month) %>%
  inner_join(
    bryonbay2
  ) %>%
  filter(Month == 4) %>%
  ggplot(aes(x = MeanT, y = MT)) + geom_point() + theme_classic() + geom_path(color = "grey")

CBclimate %>%
  rename(Year = year, Month = month) %>%
  inner_join(
    bryonbay2
  ) %>%
  filter(Month == 4) %>%
  ggplot(aes(x = Year, y = MT)) + geom_point() + theme_classic() + geom_path(color = "grey")


CBclimate2 = read.csv("Data/en_climate_monthly_NU_2400600_1929-2015_P1M.csv") %>%
  select(Year, Month, `Total.Precip..mm.`, `Mean.Temp...C.`) %>%
  rename(TP = `Total.Precip..mm.`, MT = `Mean.Temp...C.`,
         year = Year, month = Month) %>%
  filter(year > 1960 & year < 1993)


# Response function analysis ----

respCB2 = dcc_RWB(chrono = dresp, clim = CBclimate2, method = "response", start = -6, end = 9)

dcplot(respCB2)
