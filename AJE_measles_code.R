## Load necessary packages
my_packages<-c("tidyverse","magrittr","lubridate","ggplot2","dplyr","RColorBrewer",
               "mgcv","logbin","cowplot","DescTools", "ggpubr", "EnvStats",
               "hrbrthemes","viridis", "maps", "maptools", "ggmap", "rgdal", "scales", "sf")
lapply(my_packages, require, character.only = TRUE)   
theme_set(theme_classic())
remove(my_packages)

## Set working directory 
setwd("~/Documents/Data/")

#### Load and format data ####
full_df<-read.csv("AJE_measles_data.csv", stringsAsFactors=FALSE)


#### Section 1: Figs: 3,S2,S3 ####

## Subest data to just look at children under 10
## We are only interested in children over over 1 year old and under 10 for the main analysis 
df<- full_df %>% filter(full_df$FormalAge>=1 & full_df$FormalAge<10 ) 

## Create measles IgG seroprevalence and GMT figure for samples collected 2005-2015
## Add seronegative and seropos per year 
am<-table(df$SampleCollectionYear,df$MeaslesIgGSerostatus)

## Calculate measles IgG seroprevalence per year
df_prop<- data.frame(
  SampleCollectionYear=seq(2005,2015,1),
  seropos=am[,3],
  seroneg=am[,2],
  equiv=am[,1]
)
df_prop$Seropositive=(df_prop$seropos/(df_prop$seroneg+df_prop$seropos+df_prop$equiv))*100
df_prop$Seronegative=(df_prop$seroneg/(df_prop$seroneg+df_prop$seropos+df_prop$equiv))*100
df_prop$Equivocal=(df_prop$equiv/(df_prop$seroneg+df_prop$seropos+df_prop$equiv))*100

## Create functions to calculate geometric mean titer confidence intervals 
lower <- function(x, conf){
  p <- Gmean(x, conf.level = conf)
  return(p[2])
}
upper <- function(x, conf){
  p <- Gmean(x, conf.level = conf)
  return(p[3])
}

## Create a separate dataset with no seronegative individuals
df_titer<-df %>% filter(MeaslesIgGSerostatus!="negative")
df_titer_pos_over1<-df %>% filter(MeaslesIgGSerostatus!="negative") %>% filter(AgeCat>=1)
## Calculate GMT per year 
FullGMT<-df_titer_pos_over1 %>% 
  group_by(SampleCollectionYear) %>% 
  dplyr::summarize(g_mean = log10(Gmean(MeaslesIgGQuant)), lower=log10(lower(MeaslesIgGQuant,0.95)), upper=log10(upper(MeaslesIgGQuant,0.95)))
df_prop<-left_join(df_prop,FullGMT,by="SampleCollectionYear")
df_prop<- df_prop %>% mutate(samplesize=seropos+seroneg+equiv)
## Pivot longer
df_prop_long <- pivot_longer(df_prop, cols = c(Seropositive, Seronegative, Equivocal),
                             names_to = "Serostatus", values_to = "Proportion")

df_prop_long$Serostatus <- factor(df_prop_long$Serostatus, levels = c("Seronegative", "Equivocal", "Seropositive"))

## Create Figure 3 Yearly distribution of measles IgG antibody seroprevalence and geometric mean titer (GMT)
fig3 <- ggplot(df_prop_long) +geom_bar(aes(x=as.factor(SampleCollectionYear),y=Proportion, fill=Serostatus), stat="identity", width = 0.7)  + 
  geom_text(aes(x=as.factor(SampleCollectionYear),y=104, label = samplesize),  size=5, vjust=1) +
  labs(title="",
       x="Year of sample collection",
       y="Percentage (%)")   + 
  theme(plot.title = element_text(hjust = 0.5, size=25)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  scale_fill_manual(values = c("grey80","grey40", "black")) +
  theme(legend.position = "bottom",  
        legend.justification = "center")

fig3<-fig3 + geom_line(aes(x=as.factor(SampleCollectionYear),y=g_mean*25),stat="identity",color="red", group = 1, lwd=1.5) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,110,by=10),sec.axis=sec_axis(~.+0,name="Geometric Mean Titer (mIU/mL)",breaks=seq(0,110,by=12.5),labels = c(10^0,"",10^1,"",10^2,"",10^3,"",10^4))) 
## Export as landscape 7x10

## Plot age distribution by sampling year
figS2 <- ggplot(df, aes(x = as.factor(AgeCat))) +
  geom_bar(fill = "black", color = "black", width = 0.7) +
  labs(title = "",
       x = "Age Category",
       y = "Number of Samples") +
  facet_wrap(~ SampleCollectionYear) +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 0, size = 13)) +
  theme(axis.text.y = element_text(vjust = 0.6, size = 14)) +
  theme(axis.title.y = element_text(vjust = 0.6, size = 17)) +
  theme(axis.title.x = element_text(vjust = 0.6, size = 17)) +
  theme(legend.position = "none")
## Export as 7x5


## Create a plot to look at the relationship between serostatus and year
## FigS3 Measles serostatus among all individuals 2005-2015.

## Remove equivocal samples and create a binary serostatus variable
df_no_equivocal<- df %>% filter(MeaslesIgGSerostatus!="equivocal")
df_no_equivocal$Serostatus_binary<-ifelse(df_no_equivocal$MeaslesIgGSerostatus=="positive",1,0)

figS3<- ggplot(df_no_equivocal, aes(SampleCollectionYear, Serostatus_binary)) +
  geom_jitter(aes(), height=0.03) +
  geom_smooth(method=mgcv::gam,method.args=list(family=binomial)) +
  scale_x_continuous(breaks=seq(2005,2015, by=1)) +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Year",
       y="Measles serostatus \n (0:seronegative; 1:seropositive)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) 
## Export as landscape 7x5

## Remove any dataframes which aren't needed for additional sections for the sake of cleaning the global environment 
remove(am,FullGMT,df_no_equivocal,df_prop, df_prop_long)


#### Section 2. Pull in DHS Data (Figure S4a) ####

## Import 2008/2009 DHS data 
## Set working directory where data is stored 
setwd("~/Documents/Data/DHS")
dhs <- read.table("Madagascar2008.csv", sep = ",", header = TRUE)


## Calculate vaccine coverage using region 
DHS2008<- dhs %>% 
  filter(age.in.months>=9) %>% 
  dplyr::group_by(region.residence) %>%
  dplyr::summarize(m.vacc=sum(measles.y), total=n(), vacc.cov=m.vacc/total)

## Fix region names
DHS2008$region.residence[DHS2008$region.residence=="alaotra mangoro"]<-"Alaotra Mangoro"
DHS2008$region.residence[DHS2008$region.residence=="anamoroni'i mania"]<- "Amoron I Mania"
DHS2008$region.residence[DHS2008$region.residence=="analamanga"]<-"Analamanga"
DHS2008$region.residence[DHS2008$region.residence=="analanjirofo"]<-"Analanjirofo"
DHS2008$region.residence[DHS2008$region.residence=="androy"]<-"Androy"
DHS2008$region.residence[DHS2008$region.residence=="anosy"]<-"Anosy"
DHS2008$region.residence[DHS2008$region.residence=="atsimo andrefana"]<-"Atsimo Andrefana"
DHS2008$region.residence[DHS2008$region.residence=="atsimo atsinanana"]<-"Atsimo Atsinanana"
DHS2008$region.residence[DHS2008$region.residence=="atsinanana"]<-"Atsinanana"
DHS2008$region.residence[DHS2008$region.residence=="betsiboka"]<-"Betsiboka"
DHS2008$region.residence[DHS2008$region.residence=="boeny"]<-"Boeny"
DHS2008$region.residence[DHS2008$region.residence=="bongolava"]<-"Bongolava"
DHS2008$region.residence[DHS2008$region.residence=="diana"]<-"Diana"
DHS2008$region.residence[DHS2008$region.residence=="haute matsiatra"]<-"Haute Matsiatra"
DHS2008$region.residence[DHS2008$region.residence=="ihorombe"]<-"Ihorombe"
DHS2008$region.residence[DHS2008$region.residence=="itasy"]<-"Itasy"
DHS2008$region.residence[DHS2008$region.residence=="melaky"]<-"Melaky"
DHS2008$region.residence[DHS2008$region.residence=="menabe"]<-"Menabe"
DHS2008$region.residence[DHS2008$region.residence=="sava"]<-"Sava"
DHS2008$region.residence[DHS2008$region.residence=="sofia"]<-"Sofia"
DHS2008$region.residence[DHS2008$region.residence=="vakinankaratra"]<-"Vakinankaratra"
DHS2008$region.residence[DHS2008$region.residence=="vatovavy fitovinany"]<-"Vatovavy Fitovinany"
names(DHS2008)[names(DHS2008) ==  "region.residence"] <- "Regions"

## Load necessary package
library(sf)

## Troubleshooting option
# if (!require(gpclib)) install.packages("gpclib", type="source")
# gpclibPermit()


## Load in shape file data for Madagascar plots (these can be downloaded online)
setwd("/Users/arthur/Documents/Data")
sf2 <- read_sf(dsn = "./MadaMaps/data/ShpFiles/", layer = "region_mada")

# Convert to ggplot-compatible data frame
region.map <- st_as_sf(sf2) %>%
  st_sf() %>%
  st_transform(crs = 4326)  # Adjust the CRS as needed

## Load in DHS data
df.mapping <- DHS2008

## Merge region map with DHS vacc data
dff.reg <- df.mapping %>%
  mutate(region = Regions)
vacc.merge <- merge(region.map, dff.reg, by="region", all.x=TRUE)

## Create vacc x region plot 
figS4a <- ggplot() +
  geom_sf(data = vacc.merge, aes(fill = (vacc.cov * 100)), color = "black", size = 0.25) +
  scale_fill_gradient(low = "green", high = "blue", limits = c(40, 90), breaks = c(40, 50, 60, 70, 80,90),
                      guide = "colourbar", name = "") +
  theme(
    plot.margin = margin(0, 0.5, 0, 1, "cm"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm")
  )


## Remove any dataframes which aren't needed for additional sections for the sake of cleaning the global environment 
remove(dhs)




#### Section 3. Titers across regions in Madagascar (Fig: S4b-e, S5, S6) ####

## Calculate GMT for various different combinations 
df_titer$Regions<-df_titer$Province
df_titer$Regions[df_titer$Regions=="Amoron'i Mania"]<-"Amoron I Mania"
## GMT with only individuals older than 1 alive during an outbreak 
FullGMT_prov_outbreak<-df_titer %>% 
  dplyr::filter(YearOfBirth<=2005) %>%
  dplyr::filter(AgeCat>=1) %>%
  dplyr::group_by(Regions) %>% 
  dplyr::summarize(g_mean = Gmean(MeaslesIgGQuant), lower=lower(MeaslesIgGQuant,0.95), upper=upper(MeaslesIgGQuant,0.95),sample.size=n())

## GMT with only individuals older than 1 NEVER alive during an outbreak 
FullGMT_prov_no_outbreak<-df_titer %>% 
  dplyr::filter(YearOfBirth>2005) %>%
  dplyr::filter(AgeCat>=1) %>% 
  dplyr::group_by(Regions) %>% 
  dplyr::summarize(g_mean = Gmean(MeaslesIgGQuant), lower=lower(MeaslesIgGQuant,0.95), upper=upper(MeaslesIgGQuant,0.95),sample.size=n())

## Load in titer data for individuals alive during measles circulation
df.mapping <- FullGMT_prov_outbreak

## Merge region map
dff.reg <- df.mapping %>%
  mutate(region = Regions)
cases.merge.A <- merge(region.map, dff.reg, by="region", all.x=TRUE)

## Create plot 
figS4b <- ggplot() +
  geom_sf(data = cases.merge.A, aes(fill = g_mean), color = "black", size = 0.25) +
  scale_fill_gradient2(
    low = "yellow", high = "purple", limits = c(400, 3000),
    midpoint = 2000, mid = "red",
    guide = "colourbar", name = ""
  ) +
  theme(
    plot.margin = margin(0, 0.5, 0, 1, "cm"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm"))


## Assign information for geom text to be used in the following plot for the sample sizes
cases.merge.A$centroid <- st_centroid(cases.merge.A$geometry)
cases.merge.A <- cases.merge.A %>%
  mutate(
    centroid_x = st_coordinates(centroid)[, 1],
    centroid_y = st_coordinates(centroid)[, 2]
  )

## Create sample size plot 
figS5a <- ggplot() +
  geom_sf(data = cases.merge.A, aes(fill = sample.size), color = "black", size = 0.25) +
  scale_fill_gradient2(
    low = "yellow", high = "purple", limits = c(0, 400),
    breaks = c(0, 100, 200, 300, 400), midpoint = 200, mid = "red",
    guide = "colourbar", name = ""
  ) +
  geom_text(data = cases.merge.A, aes(x = centroid_x, y = centroid_y, label = sample.size), size = 4, color = "blue") +
  theme(
    plot.margin = margin(0, 0.5, 0, 1, "cm"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm")
  )

## Load in titer data for individuals not alive during measles circulation
df.mapping <- FullGMT_prov_no_outbreak

## Merge region map
dff.reg <- df.mapping %>%
  mutate(region = Regions)
cases.merge.C <- merge(region.map, dff.reg, by="region", all.x=TRUE)

## Create plot
figS4d <- ggplot() +
  geom_sf(data = cases.merge.C, aes(fill = g_mean), color = "black", size = 0.25) +
  scale_fill_gradient2(
    low = "yellow", high = "purple", limits = c(400, 3000),
    midpoint = 2000, mid = "red",
    guide = "colourbar", name = ""
  ) +
  theme(
    plot.margin = margin(0, 0.5, 0, 1, "cm"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm"))

## Assign information for geom text to be used in the following plot for the sample sizes
cases.merge.C$centroid <- st_centroid(cases.merge.C$geometry)
cases.merge.C <- cases.merge.C %>%
  mutate(
    centroid_x = st_coordinates(centroid)[, 1],
    centroid_y = st_coordinates(centroid)[, 2]
  )


## Create sample size plot
figS5b <-ggplot() +
  geom_sf(data = cases.merge.C, aes(fill = sample.size), color = "black", size = 0.25) +
  scale_fill_gradient2(
    low = "yellow", high = "purple", limits = c(0, 400),
    breaks = c(0, 100, 200, 300, 400), midpoint = 200, mid = "red",
    guide = "colourbar", name = ""
  ) +
  geom_text(data = cases.merge.C, aes(x = centroid_x, y = centroid_y, label = sample.size), size = 4, color = "blue") +
  theme(
    plot.margin = margin(0, 0.5, 0, 1, "cm"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks.length = unit(0, "cm")
  )


plot_grid(figS4b,figS4d,nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 10x4


plot_grid(figS5a,figS5b,nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 10x4


#### Plot region vacc coverage x GMT 

## Individuals alive during the measles circulation (FigS4C)
FullGMT_prov_outbreak_vacc<-left_join(DHS2008,FullGMT_prov_outbreak,by="Regions")
FullGMT_prov_outbreak_vacc$vacc.cov<-FullGMT_prov_outbreak_vacc$vacc.cov*100

figS4c<-ggplot(FullGMT_prov_outbreak_vacc, aes(x = vacc.cov, y = g_mean)) + 
  geom_point(aes(size=sample.size)) +
  geom_smooth(aes(weight = sample.size), method = lm, linewidth = 1, col = "red") +
  geom_hline(yintercept = 10^2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  # stat_smooth(method = "lm", col = "red") +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Vaccination coverage (%)",
       y="Measles Geometric Mean Titer \n (mIU/mL)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 15)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 15)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  #theme(legend.position="none") +
  scale_x_continuous(expand=c(0,0.7)) +
  scale_y_continuous(breaks=seq(0,3000, by=500),limits=c(0,3000))
#7x6

fit1<-lm(FullGMT_prov_outbreak_vacc$g_mean ~ FullGMT_prov_outbreak_vacc$vacc.cov, weights = FullGMT_prov_outbreak_vacc$sample.size)
summary(fit1)

## Individuals not alive during the measles circulation (FigS4E)
FullGMT_prov_no_outbreak_vacc<-left_join(DHS2008,FullGMT_prov_no_outbreak,by="Regions")
FullGMT_prov_no_outbreak_vacc$vacc.cov<-FullGMT_prov_no_outbreak_vacc$vacc.cov*100

figS4e<-ggplot(FullGMT_prov_no_outbreak_vacc, aes(x = vacc.cov, y = g_mean)) + 
  geom_point(aes(size=sample.size)) +
  geom_smooth(aes(weight = sample.size), method = lm, linewidth = 1, col = "red") +
  geom_hline(yintercept = 10^2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  # stat_smooth(method = "lm", col = "red") +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Vaccination coverage (%)",
       y="Measles Geometric Mean Titer \n (mIU/mL)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 15)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 15)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  #theme(legend.position="none") +
  scale_x_continuous(expand=c(0,0.7)) +
  scale_y_continuous(breaks=seq(0,3000, by=500),limits=c(0,3000))

fit2<-lm(FullGMT_prov_no_outbreak_vacc$g_mean ~ FullGMT_prov_no_outbreak_vacc$vacc.cov, weights = FullGMT_prov_no_outbreak_vacc$sample.size)
summary(fit2)


plot_grid(figS4c,figS4e,nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 14x6


## Plot age distribution for both groups
## GMT with only individuals older than 1 alive during an outbreak 
prov_outbreak<-df_titer %>% 
  dplyr::filter(YearOfBirth<=2005) 
## GMT with only individuals older than 1 NEVER alive during an outbreak 
prov_no_outbreak<-df_titer %>% 
  dplyr::filter(YearOfBirth>2005) 

figS6A<- ggplot(prov_outbreak, aes(x = as.factor(AgeCat))) +
  geom_bar(fill = "black", color = "black", width = 0.7) +
  labs(title = "",
       x = "Age Category",
       y = "Number of Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 0, size = 13)) +
  theme(axis.text.y = element_text(vjust = 0.6, size = 14)) +
  theme(axis.title.y = element_text(vjust = 0.6, size = 17)) +
  theme(axis.title.x = element_text(vjust = 0.6, size = 17)) +
  theme(legend.position = "none") +
  ylim(0,250)


figS6B<- ggplot(prov_no_outbreak, aes(x = as.factor(AgeCat))) +
  geom_bar(fill = "black", color = "black", width = 0.7) +
  labs(title = "",
       x = "Age Category",
       y = "Number of Samples") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 0, size = 13)) +
  theme(axis.text.y = element_text(vjust = 0.6, size = 14)) +
  theme(axis.title.y = element_text(vjust = 0.6, size = 17)) +
  theme(axis.title.x = element_text(vjust = 0.6, size = 17)) +
  theme(legend.position = "none") +
  ylim(0, 250)

plot_grid(figS6A,figS6B,nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 10x4


#### Section 4. Nested serosurvey for samples collected in 2005 and 2015 (Fig: 4A,4B, S7, S8) ####
df_2005 <- df %>%
  filter(AgeCat<10) %>%
  filter(AgeCat>=1) %>% 
  filter(YearOfBirth<2005) %>%
  filter(SampleCollectionYear==2005)
df_2005_pos <- df_2005 %>% filter(MeaslesIgGSerostatus!="negative")
df_2005_pos$serosurvey<-"2005"

df_2015 <- df %>%
  filter (AgeCat<10) %>%
  filter(YearOfBirth>2005) %>%
  filter(AgeCat>=1) %>% 
  filter(SampleCollectionYear==2015)
df_2015_pos <- df_2015 %>% filter(MeaslesIgGSerostatus!="negative")
df_2015_pos$serosurvey<-"2015"

df_nested <- rbind(df_2005,df_2015)
df_nested_pos <- rbind(df_2005_pos,df_2015_pos)


## Calculate measles IgG seroprevalence per year
am<-table(df_2005$AgeCat,df_2005$MeaslesIgGSerostatus)

## Calculate measles IgG seroprevalence per year
df_2005_prop<- data.frame(
  AgeCat=seq(1,9,1),
  seropos=am[,3],
  seroneg=am[,2],
  equiv=am[,1]
)
df_2005_prop$Seropositive=(df_2005_prop$seropos/(df_2005_prop$seroneg+df_2005_prop$seropos+df_2005_prop$equiv))*100


## Calculate measles IgG seroprevalence per year
am<-table(df_2015$AgeCat,df_2015$MeaslesIgGSerostatus)

## Calculate measles IgG seroprevalence per year
df_2015_prop<- data.frame(
  AgeCat=seq(1,9,1),
  seropos=am[,3],
  seroneg=am[,2],
  equiv=am[,1]
)
df_2015_prop$Seropositive=(df_2015_prop$seropos/(df_2015_prop$seroneg+df_2015_prop$seropos+df_2015_prop$equiv))*100

## Merge both data sets together
df_2005_prop$serosurvey<-"2005"
df_2015_prop$serosurvey<-"2015"
df_prob_both<-rbind(df_2005_prop,df_2015_prop)
df_prob_both$samplesize= df_prob_both$seropos + df_prob_both$seroneg + df_prob_both$equiv

## Create plot
fig4a<- ggplot(df_prob_both) +geom_bar(aes(x=as.factor(AgeCat),y=Seropositive, fill=serosurvey), stat="identity", position = 'dodge', width = 0.7) +
  geom_text(aes(x=as.factor(AgeCat),y=104, group=serosurvey, color=serosurvey, label = samplesize), size=4, vjust = 0.9, position = position_dodge(.9)) +
  scale_y_continuous(expand=c(.01,0),breaks=seq(0,110,by=10)) +
  labs(title="",
       x="Age at sample collection",
       y="Seroprevalence (%)")   + 
  theme(plot.title = element_text(hjust = 0.5, size=25)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  scale_color_manual(values = alpha(c("Black", "Dark Grey"), .8), breaks=c('2005', '2015')) +
  scale_fill_manual(values = alpha(c("Black", "Dark Grey"), .8), breaks=c('2005', '2015')) +
  guides(color="none")


## Chi-squared test
## Exclude equivocal results
df_nested_filtered <- df_nested[df_nested$MeaslesIgGSerostatus %in% c("positive", "negative"), ]
df_nested_filtered$serostatus<-ifelse(df_nested_filtered$MeaslesIgGSerostatus=="positive",1,0)
df_nested_filtered$serosurvey<-factor(df_nested_filtered$SampleCollectionYear)
chisq.test(table(df_nested_filtered$serostatus,df_nested_filtered$serosurvey))

## Figure S7 Measles serostatus among children 1-9 years from 2005 and 2015 serosurveys.
figS7<- ggplot(df_nested_filtered, aes(AgeCat, serostatus, color=serosurvey)) +
  geom_jitter(aes(), height=0.03) +
  geom_smooth(method=mgcv::gam,method.args=list(family=binomial)) +
  scale_x_continuous(breaks=seq(1,9, by=1)) +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Age at sample collection",
       y="Measles serostatus \n (0:seronegative; 1:seropositive)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  scale_color_manual(values = alpha(c("Black", "Dark Grey"), .9), breaks=c('2005', '2015'))
## Export as landscape 7x5

## Calculate GMT for each age category by serosurvey
FullSampleGMT<-df_nested_pos %>% 
  dplyr::group_by(AgeCat, serosurvey) %>% 
  dplyr::summarize(g_mean = log10(Gmean(MeaslesIgGQuant)), lower=log10(lower(MeaslesIgGQuant,0.95)), upper=log10(upper(MeaslesIgGQuant,0.95)), samplesize=n())

## Plot GMT for both nested serosurveys
fig4b <- ggplot(FullSampleGMT, aes(x = as.factor(AgeCat), y = g_mean, ymin = lower, ymax = upper, color = serosurvey)) +
  geom_text(aes(x = as.factor(AgeCat), y = 3.75, group = serosurvey, label = ifelse(!is.na(samplesize), as.character(samplesize), "")), 
            size = 4, vjust = 0, position = position_dodge(.9)) +
  geom_pointrange(aes(color = serosurvey), position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), labels = c(10^0, 10^1, 10^2, 10^3, 10^4), limits = c(0, max(FullSampleGMT$g_mean, na.rm = TRUE) * 1.1)) +
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.6)) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.6)) +
  labs(title = "",
       subtitle = "",
       caption = "",
       x = "Age at sample collection",
       y = "Measles Geometric \n Mean Titer (mIU/mL)") +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 0, size = 13)) +
  theme(axis.text.y = element_text(vjust = 0.6, size = 13)) +
  theme(axis.title.y = element_text(vjust = 0.6, size = 17)) +
  theme(axis.title.x = element_text(vjust = 0.6, size = 17)) +
  theme(legend.title = element_text(size = 17)) +
  theme(legend.text = element_text(size = 14)) +
  scale_color_manual(values = alpha(c("Black", "Dark Grey"), .8), breaks = c('2005', '2015'))

plot_grid(fig4a,fig4b, nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 14 x 5

## Create FigS8 Measles IgG antibody titer from both nested serosurveys in Madagascar.
## Plot arithmetic mean titer for both nested serosurveys
figS8a <- ggplot(df_nested_pos, aes(as.factor(AgeCat), MeaslesIgGLog))
figS8a <-figS8a + geom_boxplot(aes(fill=serosurvey)) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), labels = c(10^0, 10^1, 10^2, 10^3, 10^4), limits = c(0, 4)) +
  geom_hline(yintercept = 2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Age at sample collection",
       y="Measles Titer (mIU/mL)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  stat_compare_means(aes(group = serosurvey), method = "t.test", label = "p.signif") +
  scale_fill_manual(values = alpha(c("Black", "Dark Grey"), .6), breaks=c('2005', '2015'))
## Export landscape 8x6

figS8b<-ggplot(df_nested_pos, aes(FormalAge, MeaslesIgGLog, color=serosurvey)) +
  geom_jitter(aes()) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), labels = c(10^0, 10^1, 10^2, 10^3, 10^4), limits = c(0, 4)) +
  scale_x_continuous(breaks=seq(1,9, by=1)) +
  geom_hline(yintercept = 2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Age at sample collection",
       y="Measles Titer (mIU/mL)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  geom_smooth(method=mgcv::gam) +
  scale_color_manual(values = alpha(c("Black", "Dark Grey"), .9), breaks=c('2005', '2015'))

plot_grid(figS8a,figS8b, nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 14 x 5

## Remove any dataframes which aren't needed for additional sections for the sake of cleaning the global environment 
remove(df_2005_pos,df_2005_prop,df_2015_prop,df_2005,df_2015,df_prob_both)


#### Section 5. Extended 2005 nested serosurvey (Fig: S9) ####
##Create GMT plot for 1-14yo from 2015 survey
## Subset data from full dataframe 
df_2005_14 <- df %>% 
  filter(AgeCat<15) %>%
  filter(YearOfBirth<2005) %>%
  filter(SampleCollectionYear==2005)
df_2005_pos_14 <- df_2005_14 %>% filter(MeaslesIgGSerostatus!="negative")
df_2005_pos_14$serosurvey<-"2005"

## Calculate GMT
FullSampleGMT_14<-df_2005_pos_14 %>% 
  dplyr::group_by(AgeCat) %>% 
  dplyr::summarize(g_mean = log10(Gmean(MeaslesIgGQuant)), lower=log10(lower(MeaslesIgGQuant,0.95)), upper=log10(upper(MeaslesIgGQuant,0.95)), samplesize=n())


## Plot AMT for for 2005 serosurvey up to age 14
figS9a <- ggplot(df_2005_pos_14 %>% filter(AgeCat>0), aes(as.factor(AgeCat), log10(MeaslesIgGQuant)), color="black")
figS9a <- figS9a + geom_boxplot(aes()) +
  geom_hline(yintercept = 2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), labels = c(10^0, 10^1, 10^2, 10^3, 10^4), limits = c(0, 4)) +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Age at sample collection",
       y="Measles Titer (mIU/mL)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  theme(legend.position="none") +
  stat_n_text(vjust=-43.5)
## Export landscape 8x6

figS9b<-ggplot(df_2005_pos_14 %>% filter(AgeCat>0), aes(FormalAge, log10(MeaslesIgGQuant))) +
  geom_point() +
  geom_hline(yintercept = 2.521, linetype = "dashed", color = "red", na.rm = TRUE) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), labels = c(10^0, 10^1, 10^2, 10^3, 10^4), limits = c(0, 4)) +
  scale_x_continuous(breaks=seq(1,14, by=1)) +
  theme(axis.text.x = element_text(size=10, angle=0, vjust=0.6)) + 
  theme(axis.text.y = element_text(size=10, angle=0, vjust=0.6)) + 
  labs(title="", 
       subtitle="",
       caption="",
       x="Age at sample collection",
       y="Measles Titer (mIU/mL)")  +
  theme(plot.title = element_text(hjust = 0.5, size=0)) +
  theme(axis.text.x = element_text(vjust=0.6, angle= 0, size= 13)) +
  theme(axis.text.y = element_text(vjust=0.6, size= 13)) +
  theme(axis.title.y = element_text(vjust=0.6, size= 17)) +
  theme(axis.title.x = element_text(vjust=0.6, size= 17)) +
  theme(legend.title = element_text( size= 17)) +
  theme(legend.text = element_text( size= 14)) +
  geom_smooth(method=mgcv::gam) +
  scale_color_manual(values = alpha(c("Black", "Dark Grey"), .9), breaks=c('2005', '2015'))

plot_grid(figS9a,figS9b, nrow=1,ncol=2, align = "hv", scale=c(.98,.98))
## Export as landscape 14 x 6.5
