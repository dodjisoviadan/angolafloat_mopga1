################################################################################
# Clean ctd profiles and Clines computation on bgc data #
              # by SOVIADAN Y. Dodji adapted from Alexandre Accardo et al 2025 #
################################################################################
# install.packages("devtools")
# devtools::install_github("jiho/castr")
library(readr)
library(castr)
library(tidyverse)
library(dplyr)
library(gsw)
library(oce)
library(slider)

# csv file with BGC-ARGO data and several supplementary variables like density 
bgc_angola <- read_csv("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_angola.csv")

# select only variable of interest : 
bgc_angola <- bgc_angola %>% select(CYCLE_NUMBER, TEMP, PRES, POT_TEMP, PSAL, ABS_SAL, CHLA_ADJUSTED, DOXY_ADJUSTED, 
                                        BBP700, LATITUDE) # I don't select the density column because I will re-compute it later when data will be clean
# convert pressures values into depth
bgc_angola$depth <- swDepth(bgc_angola$PRES, bgc_angola$LATITUDE)

# First step : clean ctd cast for temperature and salinity variables
temp_sal <- bgc_angola %>% select(CYCLE_NUMBER, TEMP, POT_TEMP, PSAL, ABS_SAL, PRES, LATITUDE,depth)
# remove outliers with a global despiking :

# Create a function to detect outliers
outlier_detection <- function(x) {
  moving_median <- slide(x, k = 20, fun = median, n=1, na.rm = TRUE)
  mad <- slide(x, k = 20, fun = mad, n=1, na.rm = TRUE)
  abs_ano <- abs(x - moving_median)
  outliers <- abs_ano > 5.2*mad
  return(outliers)
}

# Apply the function to each group of data by CYCLE-NUMBER
temp_sal_outliers <- temp_sal %>%
  group_by(CYCLE_NUMBER) %>%
  mutate(outliers_temp = outlier_detection(TEMP),
         outliers_pot_temp = outlier_detection(POT_TEMP),
         outliers_sal = outlier_detection(PSAL),
         outliers_abs_sal = outlier_detection(ABS_SAL)) %>%
  ungroup()

# View what station contains outliers 

outlier_temp <- filter(temp_sal_outliers, outliers_temp == 'TRUE') # just to check 
outlier_sal <- filter(temp_sal_outliers, outliers_sal == 'TRUE')

# it seems to be ok ! 

# compare for some selected profiles if this method applied to the all data set give the same result than if it is applied on only one profile  
# select only the profile 3 and temperature ----

station3 <- filter(temp_sal, CYCLE_NUMBER == 3) %>% select(CYCLE_NUMBER, POT_TEMP, depth) # %>% na.omit()

station3 <-station3 %>% 
  mutate(moving_median = slide(POT_TEMP, k = 10, fun = median, n=10, na.rm = TRUE), 
         mad = slide(POT_TEMP, k = 10, fun = mad, n=1, na.rm = TRUE),
         abs_ano_moving_median = abs(POT_TEMP-moving_median),
         outliers = abs_ano_moving_median > 5.2*mad)

# check for salinity in station 2 
station2 <- filter(temp_sal, CYCLE_NUMBER == 2) %>% select(CYCLE_NUMBER, PSAL, depth) # %>% na.omit()

station2 <-station2 %>% 
  mutate(moving_median = slide(PSAL, k = 10, fun = median, n=10, na.rm = TRUE), 
         mad = slide(PSAL, k = 10, fun = mad, n=1, na.rm = TRUE),
         abs_ano_moving_median = abs(PSAL-moving_median),
         outliers = abs_ano_moving_median > 5.2*mad)

# It gives the same result ! ----

# now replace all outliers by NA 

# Replace values in the POT_TEMP column
temp_sal_outliers$POT_TEMP <- ifelse(temp_sal_outliers$outliers_pot_temp == TRUE, NA, temp_sal_outliers$POT_TEMP)
# Replace values in the TEMP column
temp_sal_outliers$TEMP <- ifelse(temp_sal_outliers$outliers_temp == TRUE, NA, temp_sal_outliers$TEMP)

# Replace values in the PSAL column
temp_sal_outliers$PSAL <- ifelse(temp_sal_outliers$outliers_sal == TRUE, NA, temp_sal_outliers$PSAL)

# Replace values in the ABS_SAL column
temp_sal_outliers$ABS_SAL <- ifelse(temp_sal_outliers$outliers_abs_sal == TRUE, NA, temp_sal_outliers$ABS_SAL)


# create a new data set with cleaned values 
bgc_angola_clean <- tibble(CYCLE_NUMBER = temp_sal_outliers$CYCLE_NUMBER,
                             TEMP= temp_sal_outliers$TEMP,
                             POT_TEMP = temp_sal_outliers$POT_TEMP,
                             ABS_SAL = temp_sal_outliers$ABS_SAL,
                             PSAL = temp_sal_outliers$PSAL, 
                             PRES = temp_sal_outliers$PRES, 
                             CHLA_ADJUSTED = bgc_angola$CHLA_ADJUSTED, 
                             BBP700 = bgc_angola$BBP700, 
                             DOXY_ADJUSTED = bgc_angola$DOXY_ADJUSTED, 
                             Latitude = bgc_angola$LATITUDE)

# convert pressures values into depth
bgc_angola_clean$depth <- swDepth(bgc_angola_clean$PRES, bgc_angola_clean$Latitude)
bgc_angola_clean$PSAL[35418:35438]<-NA
bgc_angola_clean$PSAL<-ifelse((bgc_angola_clean$PSAL<33)==TRUE,NA,bgc_angola_clean$PSAL)
# plot profiles to see if it's clean 
ggplot(bgc_angola_clean) + geom_path(aes(x=TEMP, y=-PRES, group=CYCLE_NUMBER), alpha=0.6)# good
ggplot(bgc_angola_clean) + geom_path(aes(x=PSAL, y=-PRES, group=CYCLE_NUMBER), alpha=0.6)# good
#bgc_angola_clean<-bgc_angola_clean%>%  filter(PSAL > 33)

# now compute potential density based on ABS_SAL and TEMP

bgc_angola_clean$dens <- gsw_pot_rho_t_exact(bgc_angola_clean$ABS_SAL, bgc_angola_clean$TEMP, bgc_angola_clean$PRES, p_ref = 1)

# plot it 
ggplot(bgc_angola_clean) + geom_path(aes(x=dens, y=-PRES, group=CYCLE_NUMBER), alpha=0.6) # good

# Despike each profile
detach("package:oce", unload = TRUE) # oce package mask despike function 
bgc_angola_clean_despike <- bgc_angola_clean %>% group_by(CYCLE_NUMBER) %>%
  mutate(TEMP=despike(TEMP, mult=3),
         POT_TEMP=despike(POT_TEMP, mult=3),
         PSAL=despike(PSAL, mult=3),
         ABS_SAL=despike(ABS_SAL, mult=3),
         DOXY_ADJUSTED=despike(DOXY_ADJUSTED, mult=3),
         CHLA_ADJUSTED=despike(CHLA_ADJUSTED, k=2, mult=4),
         BBP700=despike(BBP700, k=2, mult=4),
    # NB: chlorophyll and bbp are more spiky so be less stringent
         dens =despike(dens, mult=3))

compare <- tibble(before_despiking = bgc_angola_clean$TEMP,
                      after_despiking = bgc_angola_clean_despike$TEMP)

summary(compare)
compare <- compare[is.na(compare$after_despiking),]

# check if these methods applied to the all dataset give the same result than if it is applied to one profile and one variable ----

station1 <- bgc_angola_clean %>% filter(CYCLE_NUMBER == 1) %>% select(CYCLE_NUMBER, DOXY_ADJUSTED, depth)
station1 <- station1 %>% mutate(DOXY_ADJUSTED = despike(DOXY_ADJUSTED, mult=3))

station1_despike <- bgc_angola_clean_despike %>% filter(CYCLE_NUMBER == 1) %>% select(CYCLE_NUMBER, DOXY_ADJUSTED, depth)

identical(station1$DOXY_ADJUSTED, station1_despike$DOXY_ADJUSTED)

# it gives the same results ! ----

ggplot(bgc_angola_clean_despike) + geom_path(aes(x=dens, y=-depth, group=CYCLE_NUMBER), alpha=0.6, na.rm = TRUE)


# now binning : 5 meters bins and compute the mean of each variable inside
bgc_angola_binned <-  bgc_angola_clean_despike %>% mutate(binned_depth = round(depth/5)*5 + 2.5)
bgc_angola_binned <-  bgc_angola_binned %>% group_by(CYCLE_NUMBER, binned_depth) %>%
  summarise(TEMP = mean(TEMP, na.rm = TRUE), 
            POT_TEMP = mean(POT_TEMP, na.rm = TRUE),
            ABS_SAL= mean(ABS_SAL, na.rm = TRUE),
            PSAL = mean(PSAL, na.rm = TRUE),
            CHLA_ADJUSTED = mean(CHLA_ADJUSTED, na.rm = TRUE),
            DOXY_ADJUSTED = mean(DOXY_ADJUSTED, na.rm = TRUE), 
            BBP700 = mean(BBP700, na.rm = TRUE),
            dens = mean(dens, na.rm = TRUE))

# save it 

write.csv(bgc_angola_binned, "C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_angola_clean.csv", row.names = FALSE)


# now compute some variables for each station MLD, DCM, stratification index, sotck of chla (within 10m of the DCM) 

ctd_clines <- bgc_angola_binned %>% group_by(CYCLE_NUMBER) %>%
  summarise(thermocline = clined(POT_TEMP, binned_depth, n.smooth=0), #k = 2
            pycnocline = clined(dens, binned_depth),
            strat_index = stratif(dens, binned_depth, min.depths=0:5, max.depth=200:250),
            DCM = maxd(CHLA_ADJUSTED, binned_depth, n.smooth=2, k=3),
            MLD = mld(dens, binned_depth, ref.depths=0:5, default.depth=80),
            # average temperature in the mixed layer only
            temp_mld = integrate(TEMP, binned_depth, from=0, to=MLD, fun=mean),
            # stock of Chl a within 10 m of the DCM
            chla_dcm_stock = integrate(CHLA_ADJUSTED, binned_depth, from=DCM-10, to=DCM+10))

# check if these methods applied to the all data set give the same result than if it is applied to one profile and one variable ----
station1 <- bgc_angola_binned %>% filter(CYCLE_NUMBER == 1) %>% select(binned_depth, CYCLE_NUMBER, dens) %>% group_by(CYCLE_NUMBER) %>% summarise(strat_index = stratif(dens, binned_depth, min.depths=0:5, max.depth=200:250))
# It gives the same results ----



# Try to modify a little the clined function in order to set up a threshold for the maximum standard deviation. If the max sd is bellow this threshold the return depth is NA 
check_input <- function(x, depth=NULL) {
  ok <- TRUE
  # check the input
  if (all(is.na(x))) {
    ok <- FALSE
  }
  if (!is.null(depth)) {
    if (length(depth) != length(x)) {
      ok <- FALSE
      stop("The vector of data (n=", length(x), ") should be as long as the vector of depths (n=", length(depth), ")")
    }
  }
  return(ok)
}
get_depth <- function(i, depth) {
  if (length(i) > 0) {
    if (!is.null(depth)) {
      i <- depth[i]
    }
  } else {
    i <- NA
  }
  return(i)
}
modified_clined <- function(x, depth=NULL, n.smooth=0, k=2, threshold=NULL) {
  # check input
  ok <- check_input(x, depth)
  if (!ok) { return(NA) }
  
  # smooth the profile (if requested)
  x <- smooth(x, k=k, n=n.smooth)
  
  # compute the standard deviation
  s <- slide(x, k=k, stats::sd, na.rm=TRUE)
  # get its maximum
  max_s <- max(s)
  
  # check if the maximum standard deviation is below the threshold
  if (!is.null(threshold) && max_s < threshold) {
    return(NA)
  }
  
  # get the index of the maximum standard deviation
  i <- which.max(s)
  
  # if the depth is provided, extract the corresponding depth
  i <- get_depth(i, depth)
  
  return(i)
}

# try on temperature to see if it's work for thermocline 
new_ctd_clines <- bgc_angola_binned %>% group_by(CYCLE_NUMBER) %>%
  summarise(thermocline = modified_clined(POT_TEMP, binned_depth, n.smooth=2, k = 2, threshold = 0.27)) 
          
# plot time serie of each variable 
# merge with original data set 
bgc_angola <- read_csv("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/bgc_angola.csv")  %>% select(CYCLE_NUMBER, Date_Time) %>% group_by(CYCLE_NUMBER) %>% summarise(Date_Time = mean(Date_Time)) %>% merge(bgc_angola_binned, by = "CYCLE_NUMBER")



bgc_angola <- merge(bgc_angola, ctd_clines, by = "CYCLE_NUMBER")

ggplot(bgc_angola) + geom_path(aes(x = Date_Time, y = -thermocline))
# Filter data frame based on interval value in 'Date_Time' column in order to plot only specific profiles and check if there are valid
filtered_profiles <- subset(bgc_angola, CYCLE_NUMBER >= 178 & CYCLE_NUMBER <= 185)
ggplot(filtered_profiles) + geom_path(aes(x=POT_TEMP, y=-binned_depth))+
  geom_hline(aes(yintercept = -thermocline), color = "red")+
  facet_wrap(~CYCLE_NUMBER, ncol = 4)+
  theme_bw()

# dubious profiles 
dubious_p <- c(82:84, 114, 121, 124, 125, 128, 132:141, 143, 145, 147, 148, 151:155, 158, 159, 162, 164, 167:169, 171, 173, 175:180, 182, 184, 185)
# select them 
dubious_df <- bgc_angola[bgc_angola$CYCLE_NUMBER %in% dubious_p, ]

# plot them 
filtered_profiles <- subset(dubious_df, CYCLE_NUMBER >= 179 & CYCLE_NUMBER <= 185)
ggplot(filtered_profiles) + geom_path(aes(x=POT_TEMP, y=-binned_depth))+
  geom_hline(aes(yintercept = -thermocline), color = "red")+
  facet_wrap(~CYCLE_NUMBER, ncol = 3)+
  theme_bw()

# time serie without these profiles 
good_profiles <- bgc_angola[!(bgc_angola$CYCLE_NUMBER %in% dubious_p), ]
ggplot(good_profiles) + geom_path(aes(x = Date_Time, y = -thermocline))+ geom_point(aes(x = Date_Time, y = -thermocline), size = 0.3) + geom_smooth(aes(x = Date_Time, y = -thermocline))
# better without 

# There is a problem in the method that is used to compute the clines 

# plot MLD and stratification index in the same plot to see ----

ylim.prim <- c(0, 350)
ylim.sec <- c(0, 3)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1]-ylim.sec[1])

ggplot(bgc_angola)+
  geom_path(mapping = aes(x = Date_Time, y = MLD, color = "red"))+
  geom_smooth(mapping = aes(x = Date_Time, y = MLD, color = "red"))+
  geom_path(aes(x=Date_Time, y=a + strat_index*b, color = "blue"))+
  geom_smooth(aes(x=Date_Time, y=a + strat_index*b, color = "blue"))+
  scale_y_continuous(sec.axis = sec_axis(~(. - a)/b, name = "Stratification index"))+
  scale_color_manual(values = c("red" = "red", "blue" = "blue"), 
                     labels = c("MLD", "Stratification index"))+
  theme_bw()
# ----

write.csv(ctd_clines, 'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/station_clines.csv', row.names = FALSE)
