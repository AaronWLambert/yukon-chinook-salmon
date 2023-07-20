#=================================================================================
# Project Name: Yukon River INSEASON FORECAST - Integrated Bayesian Inseason Model
# Creator: Aaron Lambert, College of Fisheries and Ocean Sciences, UAF
# Date: 03.28.22
# 
# Purpose: This is a working script to generate figures for papers and presentations
library(tidyverse)
library(ggthemes)
library(wesanderson)
# Define Workflow Paths ============================================
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")

dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")

######### Import Data ###############
# Historical Canadian EOS reconstructed run
# CAN_hist <- readRDS(file.path(dir.data,"/EOS Reconstructed Canadian.RDS"))
CAN_hist <- readRDS(file.path(dir.data,"Reconstructed_CAN_Yukon_Chinook_Run_1995_2019.RDS"))
str(CAN_hist)
# This is the reconstructed data from Curry for old reconstructed modeling procedure
# CAN_hist <- readRDS(file.path(dir.data,"Can Origin Reconstructed 2Feb22.RDS"))
# PSS historical for Jun 1st = 1, to Aug 31 = 92
# Years 1995 - 2020
# PSScum_hist <- readRDS(file.path(dir.data,"/clean_counts.rds"))

# Read in PSS non-cummulative with date in date format
PSS_hist <- readRDS(file.path(dir.data,"/PSS Chinook Passage non-cumm.RDS"))
# Read in historical avg of GSI by strata (Naive estimator)
GSI_mean<-readRDS(file = file.path(dir.data,"Mean GSI by strata 2005-2020.RDS"))

# Read in historic preseason forecasts (2013 - current)
pf_hist <- readRDS(file.path(dir.data,"preseason forcast.RDS"))
# pf_hist <- readRDS(file.path(dir.data,"inv_var_weighted_forcast_v3_Jan282022.RDS"))
# Read in genetic stock identification (2005-2019) 
# (adjusted to capture early and late runs)
GSI_by_year <- readRDS(file = file.path(dir.data,"GSI by year"))

# Read in harvest data for US portion of Yukon
usYukon <- readRDS(file = file.path(dir.data,"Alaska_Harvest_of_Yukon_River_Chinook_Salmon_1961_2020.RDS"))

# Read in Canadian Harvest of Yukon kings
canHarvest <- readRDS(file = (file.path(dir.data,"Canadian Harvest 1961_2020.RDS")))


# Get average harvest 
Can_avg <- canHarvest  %>% summarise(Mean2006 = mean(total)) 

pf_hist$Year <- as.double(pf_hist$Year)
# Join can_hist and pf_hist for plotting
full_abund <-full_join(CAN_hist,pf_hist, by = "Year")
# Plot to compare previous predictions to realized runsizes
ggplot(full_abund, aes(x = Year,
                     y = can.mean,
                     ymin = can.mean - can.sd,
                     ymax = can.mean + can.sd))+
  geom_point(fill = "red")+
  geom_errorbar()+
  geom_point( aes(x = Year, y = Weighted_Forecast), color = "red")

ggplot(pf_hist, aes(x = Year, ymin = Low/1000, ymax = High/1000, y = PostSeasonEstimate/1000))+
  geom_point(col = "red")+
  geom_errorbar()+
  labs(y = "1000's of Chinook Salmon")


# PLots for harvest
#
# Create Country variable for each df
canHarvest$Country <- "Canada"
usYukon$Country <- "USA"
canHarvest$Country <- as.factor(canHarvest$Country)
usYukon$Country <- as.factor(usYukon$Country)

longHarvest_US <- usYukon %>% pivot_longer(cols = -c(Year,Country)) %>% as.data.frame()
longHarvest_Can <- canHarvest%>% pivot_longer(cols = -c(Year,Country)) %>% as.data.frame()
longHarvest_Can
fullHarvestDF <- rbind(longHarvest_Can,longHarvest_US)

totalHarvDF <- fullHarvestDF %>% 
  group_by(Year,Country) %>% 
  summarise(TotalHarvest = sum(value)) %>% 
  as.data.frame()


ggplot(totalHarvDF, aes(x = Year, y = TotalHarvest/1000, fill = Country))+
  geom_col()+
  labs(y = "1,000's of Chinook Salmon",
       fill = "")+
  scale_fill_colorblind()+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank()
  )

# Average line calc
US_avg <- longHarvest_US %>% group_by(Year) %>% summarise(Total = sum(value))
US_avg <- mean(US_avg$Total)
ggplot(longHarvest_US, aes(x = Year, y = value/1000, fill = name))+
  geom_col()+
  scale_fill_manual(values = wes_palette("Rushmore1"))+
  labs(y = "1,000's of Chinook Salmon",
       fill = "",
       title = "U.S. Harvest")+
  geom_hline(yintercept = US_avg/1000, linetype = 6, size = 2, color = "red")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle =  90)
          )+
  scale_x_continuous(breaks = seq(1961,2020,2))

# Get ave harvest for geom hline
Can_avg <- longHarvest_Can %>% group_by(Year) %>% summarise(Total = sum(value))
Can_avg <- mean(Can_avg$Total)
# Canadian harvest
ggplot(longHarvest_Can, aes(x = Year, y = value/1000, fill = name))+
  geom_col()+
  scale_fill_manual(values = wes_palette("IsleofDogs1"))+
  labs(y = "1,000's of Chinook Salmon",
       fill = "",
       title = "Can Harvest")+
  geom_hline(yintercept = Can_avg/1000, linetype = 6, size = 2, color = "red")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1, color = "grey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle =  90))+
  scale_x_continuous(breaks = seq(1961,2020,2))

ggplot(fullHarvestDF, aes(x = Year, y = value/1000, fill = factor(name)))+
  geom_col()+
  labs(x = "Year",
       y = "Number of Chinook Salmon",
       fill = "Harvest Catagory")+
  theme(legend.position = "top",
        panel.grid.major.y  =element_line(linetype = 2, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
  )+
  facet_wrap(~Country )


# Bayesian Explantation plots ###################################################################


# Normal Curve for presentation

prior <- rnorm(1e6, 60,20)
PSS <- rnorm(1e6,100,10)
posterior <- rnorm(1e6, 75, 15)
uniform <- runif(1e6,30,170)

prior.df <- data.frame("name" = "Prior", "value" = prior)
pss.df <- data.frame("name" = "PSS", "value" = PSS )
posterior.df <- data.frame("name" = "Posterior", "value" = posterior)
unif.df <- data.frame("name" = "Uniform", "value" = uniform)

df <- rbind(prior.df,pss.df,posterior.df, unif.df)

df$name <- factor(df$name, levels = c("Prior","PSS","Posterior","Uniform"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(df[df$name=="Prior",], aes(x = value, fill = name))+
  geom_density(alpha = .9, fill = cbPalette[2])+
  coord_cartesian(xlim = c(0, 200),
                  ylim = c(0,.035))+
  labs(y = "Relative Probability",
       x = "Thousands of Canadian Chinook Salmon")+
  theme_classic()+
  theme(legend.title = element_blank(),
        text = element_text(size = 20))

ggplot(df[df$name ==c("PSS","Prior") ,], aes(x = value, fill = name))+
  geom_density(alpha = .9, show.legend = F)+
  coord_cartesian(xlim = c(0, 200),
                  ylim = c(0,.035))+
  labs(y = "Relative Probability",
       x = "Thousands of Canadian Chinook Salmon")+
  scale_fill_manual(values = cbPalette[2:3])+
  theme_classic()+
  theme(legend.title = element_blank(),
        text = element_text(size = 20))


ggplot(df, aes(x = value, fill =name))+
  geom_density(alpha = .7, show.legend = F)+
  scale_fill_manual(values   = cbPalette[2:4])+
  coord_cartesian(xlim = c(0, 200),
                  ylim = c(0,.035))+
  labs(x = "Thousands of Canadian Chinook Salmon",
       y = "Relative Probability")+
  theme_classic()+
  theme(text = element_text(size = 20))

df_sub <- df$value[df$name=="PSS"]

qt<-quantile(df_sub,probs = c(.025,.975))
qt[1]
# Uniform vs prob dens
ggplot(df[df$name==c("PSS"),],
       aes(x = value, fill = name, color = name))+
  geom_density( show.legend = F)+
  geom_density(data = df[df$name==c("PSS")& 
                           df$value >= qt[1] 
                         & df$value<=qt[2],],
               fill = cbPalette[3],
               show.legend = F,
               alpha = .5)+
  # coord_cartesian(xlim = c(0, 200),
  #                 ylim = c(0,.02))+
  labs(y = "Relative Probability",
       x = "Thousands of Canadian Chinook Salmon")+
  scale_fill_manual(values = cbPalette[2])+
  scale_color_manual(values = cbPalette[2])+
  theme_classic()+
  theme(legend.title = element_blank(),
        text = element_text(size = 20))

# Credible interval Plot ##################################
set.seed(100)
PSS <- data.frame(PSS=rnorm(100000, 100, 20))

# mean1 <- mean(PSS$PSS)
rand1 <- quantile(PSS$PSS, probs = c(.025,.975))

p <- ggplot(PSS, aes(PSS)) +
  geom_density(fill="grey", color = "Black") 

d <- ggplot_build(p)$data[[1]]

p + geom_area(data = subset(d, x > rand1[1] & x<rand1[2]), 
              aes(x=x, y=y),
              color = "Black",
              fill=cbPalette[2] ,
              alpha = .5) +
  geom_segment(x=rand1[1], xend=rand1[1], 
               y=0, yend=approx(x = d$x, y = d$y, xout = rand1[1])$y,
               colour= cbPalette[3] , size=3)+
  geom_segment(x=rand1[2], xend=rand1[2], 
               y=0, yend=approx(x = d$x, y = d$y, xout = rand1[2])$y,
               colour= cbPalette[3] , size=3)+
  geom_segment(x = rand1[1], xend = rand1[2],
               y = .5*approx(x = d$x, y = d$y, xout = rand1[1])$y,
               yend= .5*approx(x = d$x, y = d$y, xout = rand1[2])$y,
               color = cbPalette[3], size = 3 )+
  
  labs(y = "Relative Probability",
       x = "Thousands of Canadian Chinook Salmon")+
  theme_classic()+
  theme(legend.title = element_blank(),
        text = element_text(size = 20))
  






#
ggplot(df[df$name==c("PSS"),], aes(x = value, fill = name))+
  geom_density(alpha = .5, show.legend = F)+
  coord_cartesian(xlim = c(0, 200),
                  ylim = c(0,.02))+
  labs(y = "Relative Probability",
       x = "Thousands of Canadian Chinook Salmon")+
  scale_fill_manual(values = cbPalette[3])+
  theme_classic()+
  theme(legend.title = element_blank(),
        text = element_text(size = 20))


# Density plots comparing pf, linear prediction, and posterior estimate
ggplot(df.prior, aes(x = value/1000, fill = par))+
  geom_density(alpha = .5)+
  # ggtitle(paste("Density plot for",myYear,", day",myDay,"\n Version", model.version))+
  # theme_classic()+
  # geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
  #                color = "End of Season Runsize"),
  #            linetype = 2,
  #            size = 1)+
  
  # xlim(c(0,300000))+
  coord_cartesian(xlim = c(0,250))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", fill = "Parameters", y = "Relative Probability")+
  # scale_fill_discrete(name = "Parameters", 
  # labels = c( "Preseason Forecast (Prior)", 
  # "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  # scale_fill_colorblind(name = "", 
  #                       labels = c( "Preseason Forecast (Prior)", 
  #                                   "PSS Prediction","Runsize"))+
  # theme_tidybayes(element_text(size = 20))+
  theme(text = element_text(size = 20, family = "serif"),
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")

# Preseason transition to inseason regression inforgraphic ######################
day <- 148:300

vect <- vector(length = length(day))
count <- 0
# vect[1:7] <- 1
for (d in 1:length(day)) {
  
  
  vect[d]<- 1-(1/(1+exp(-((d-9)-15)/12)))
  
  count <- count +.0011
}
vect[1:5] <- c(1,1,1,1,1)
df <- data.frame("day" = day[1:50], "pf" = vect[1:50])
df$PSS <- 1-df$pf
df
df_long <- pivot_longer(df, cols = -day)


ggplot(df_long, aes(x = day, y = value*100, fill = name) )+
  geom_area(show.legend = F)+
  labs(x = "Date", y = "Model Belief")+
  # scale_fill_discrete(name = "", labels = c("Preseason Forecast",
  # "PSS Sonar"))+
  scale_x_continuous(breaks = seq(152,197,7), labels = c("June 1",
                                                         "June 10",
                                                         "June 20",
                                                         "June 30",
                                                         "July 10",
                                                         "July 20",
                                                         "July 30"))+
  scale_y_continuous(breaks = c(0,25,50,75,100), labels = c("0%",
                                                            "25%",
                                                            "50%",
                                                            "75%",
                                                            "100%"))+
  scale_fill_colorblind()+
  theme(axis.text = element_text(size = 20),
        text = element_text(size = 20))

# Regression plots from model outputs ###############################3
# Compare how slope and CI change over year

#Import a retro model output
outputlist <- readRDS(file = file.path(dir.output,"OutPut_ver211_oldCan_2007_2010_2013_2021_3Nov22.RDS"))


# Put plots in one image
par(mfrow=c(1,3),mar = c(1, 2, 2, 2), oma = c(5,5,2,2))
# Day 158
# Get parameters for plots
pars <- outputlist$`2016_163`$pars
cumPSS <- outputlist$`2016_163`$cumPSS
totalEOS <- outputlist$`2016_163`$totalEOS

# Predict values to plot to axis

x <- seq(from = 0, to=3e5, by= 1000)

mat <- matrix(nrow = length(pars$alpha), ncol = length(x))

for (c in 1:length(x)) {
  
  for (r in 1:length(pars$alpha)) {
  

mat[r,c] = exp(pars$alpha[r]+ (pars$beta[r]*x[c]))

  }
}
head(mat)

quant.predPSS <- apply(X = mat, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# Prder for polygon
ord <- order(x)

# Regression plot for non-adjusted counts
plot(x = cumPSS/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = "",
     ylab = "",
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     # xlim = c(0,max(cumPSS)/1000),
     # ylim = c(min(quant.predPSS/1000), max(quant.predPSS/1000)),
     ylim = c(0,350),
     xlim=c(0,300),
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(x[ord]/1000, rev(x[ord]/1000)),
        y = c(quant.predPSS[1,ord]/1000,rev(quant.predPSS[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(x[ord]/1000, rev(x[ord]/1000)),
        y = c(quant.predPSS[2,ord]/1000,rev(quant.predPSS[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = x[ord]/1000, quant.predPSS[3,ord]/1000, col = "red", lw = 2)



# Day 178
# Get parameters for plots
pars <- outputlist$`2016_178`$pars
cumPSS <- outputlist$`2016_178`$cumPSS
totalEOS <- outputlist$`2016_178`$totalEOS

# Predict values to plot to axis

x <- seq(from = 0, to=3e5, by= 1000)

mat <- matrix(nrow = length(pars$alpha), ncol = length(x))

for (c in 1:length(x)) {
  
  for (r in 1:length(pars$alpha)) {
    
    
    mat[r,c] = exp(pars$alpha[r]+ (pars$beta[r]*x[c]))
    
  }
}
head(mat)

quant.predPSS <- apply(X = mat, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# Prder for polygon
ord <- order(x)

# Regression plot for non-adjusted counts
plot(x = cumPSS/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = "",
     ylab = "",
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     # xlim = c(0,max(cumPSS)/1000),
     # ylim = c(min(quant.predPSS/1000), max(quant.predPSS/1000)),
     ylim = c(0,350),
     xlim=c(0,300),
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(x[ord]/1000, rev(x[ord]/1000)),
        y = c(quant.predPSS[1,ord]/1000,rev(quant.predPSS[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(x[ord]/1000, rev(x[ord]/1000)),
        y = c(quant.predPSS[2,ord]/1000,rev(quant.predPSS[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = x[ord]/1000, quant.predPSS[3,ord]/1000, col = "red", lw = 2)

# Day 
# Get parameters for plots
pars <- outputlist$`2016_203`$pars
cumPSS <- outputlist$`2016_203`$cumPSS
totalEOS <- outputlist$`2016_203`$totalEOS

# Predict values to plot to axis

x <- seq(from = 0, to=3e5, by= 1000)

mat <- matrix(nrow = length(pars$alpha), ncol = length(x))

for (c in 1:length(x)) {
  
  for (r in 1:length(pars$alpha)) {
    
    
    mat[r,c] = exp(pars$alpha[r]+ (pars$beta[r]*x[c]))
    
  }
}
head(mat)

quant.predPSS <- apply(X = mat, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
# Prder for polygon
ord <- order(x)

# Regression plot for non-adjusted counts
plot(x = cumPSS/1000, y = totalEOS/1000,
     type = "p",
     pch = 21,
     bg = "red",
     xlab = "",
     ylab = "",
     # main = paste("Predicted PSS Fit","1995 to", 
     # myYear-1,"\n Day 152 to",
     # myDay,"\n Model Version",model.version),
     # xlim = c(0,max(cumPSS)/1000),
     # ylim = c(min(quant.predPSS/1000), max(quant.predPSS/1000)),
     ylim = c(0,350),
     xlim=c(0,300),
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(x[ord]/1000, rev(x[ord]/1000)),
        y = c(quant.predPSS[1,ord]/1000,rev(quant.predPSS[5,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(x = c(x[ord]/1000, rev(x[ord]/1000)),
        y = c(quant.predPSS[2,ord]/1000,rev(quant.predPSS[4,ord]/1000)),
        col=rgb(1,0,0, alpha=0.2), border=FALSE)
lines(x = x[ord]/1000, quant.predPSS[3,ord]/1000, col = "red", lw = 2)

mtext("EOS Canadian Runsize (1000's Salmon)", 
      side = 2,
      at = 0,
      adj = 0, 
      line = 2, 
      font = 3, 
      col = "Black", 
      cex = 1.5,
      outer = T)

mtext("Cummulative PSS Passage (1000's Salmon)",
      side = 1,
      at = .27, 
      adj = 0, 
      line = 2, 
      font = 3, 
      col = "Black", 
      cex = 1.5,
      outer = T)

mtext("PSS Projections for 2016", 
      side = 3,
      at = -375, 
      adj = 0, 
      line = 1,
      font = 3, col = "Black", cex = 1.5,
      outer = F)


# Bar plot comparing PF and EOS counts ############################################
#  **Need to run Yukon Inseason Forecast script first**
df <- data.frame("PF" = PF_vect, "EOS" = EOS)
df$Year <- c(2007:2010,2013:2021)

df_long <- df %>% pivot_longer( cols = -Year) %>% as.data.frame()

ggplot(df_long, aes( x=  Year, y = value, fill = name))+
  geom_col(position = "dodge")+
  scale_fill_discrete(name = "", labels = c("Runsize","Preseason Forecast"))+
  scale_x_continuous(breaks = seq(2007,2021,1))+
  labs(x = "Year", y = "Runsize")+
  theme_tidybayes()+
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 18),
        legend.position = "top")

ggplot(PSS_hist[PSS_hist$Year<=2015 &
                  PSS_hist$Day <= 187,], aes(x = Day, y = count/1000 ))+
  geom_col()+
  facet_wrap(~Year)+
  labs(y = "Thousands Chinook Salmon",
       x = "Day of Year")+
  scale_x_continuous(breaks = seq(152,188,9),
                     labels = c("June 1",
                                "June 10",
                                "June 19",
                                "June 28",
                                "July 7"))+
  theme_pubclean()+
  theme(axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 18)
  )