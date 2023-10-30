# Function for generating density plots




density.func <- function(outputList, CAN_hist ){

  
  pars <- outputList$pars
  myYear <- outputList$myYear
  myDay  <- outputList$myDay
 



# Extract estimates for prior, likelihood, and posterior 
df.run <- data.frame("par" = "Runsize",
                     "value" = pars$RunSize)

# PF (prior)
df.prior <- data.frame("par" = "Prior", 
                       "value" = pars$prior_pf)

# PSS prediction density
df.postPredPss <- data.frame("par" = "PSS Prediction", 
                             "value" = pars$post_curr_predPSS)

# Eagle prediction density
df.postPredEagle <- data.frame("par" = "Eagle Prediction", 
                               "value" = pars$post_curr_predEagle)

# Bind into one data frame with 2 columns for plotting
df.comb <- rbind(df.prior,
                 df.postPredPss,
                 df.postPredEagle,
                 df.run)
# df.comb$par <- levels(df.comb$par, levels = c("Prior","PSS_post_pred","Runsize"))
# Density plots comparing pf, linear prediction, and posterior estimate
postDense <-ggplot(df.comb, aes(x = value/1000, fill = par))+
  geom_density(alpha = .5)+
  geom_vline( aes(xintercept = median(pars$RunSize)/1000,
                  col = "Median Projection"),
              linetype = 2,
              size = 1.5)+
  geom_vline( aes(xintercept = CAN_hist$can.mean[CAN_hist$Year==myYear]/1000,
                  col = "True Runsize"),
              linetype = 2,
              size = 1.5)+
  # ylab("Relative Probability")+
  ylab("")+
  xlab("")+
  # ggtitle(paste("Version = ",Ver,"\n Day = ",myD,"\n Year = ",myY, sep = ""))+
  # xlab("1000's of Chinook Salmon")+
  coord_cartesian(xlim = c(0,200))+
  scale_fill_colorblind(name = "")+
  theme_tidybayes()+
  theme(text = element_text(size = 20, family = "serif"),plot.title = element_text(size = 12))+
  scale_color_discrete(name = "")

eyeDense <- ggplot(df.comb, aes(x = value/1000, y =fct_rev( par)))+
  stat_eye(aes(fill = after_stat(level)), .width = c(.8,.95,1), alpha = .5)+
  geom_vline(aes(xintercept = CAN_hist$can.mean[CAN_hist$Year == myYear]/1000,
                 color = "End of Season Runsize"),
             linetype = 2,
             size = 1)+
  coord_cartesian(xlim = c(20,300))+
  # ylim(c(0,500))
  labs(x = "1000's of Chinook Salmon", 
       fill = "Parameters", 
       y = "Parameter")+
  # scale_fill_discrete(name = "Parameters", 
  # labels = c( "Preseason Forecast (Prior)", 
  # "PSS Prediction","Runsize"))+
  scale_x_continuous()+
  # scale_fill_colorblind(name = "", 
  #                       labels = c( "Preseason Forecast (Prior)",
  #                                   "PSS Prediction",
  #                                   "Runsize"))+
  theme_tidybayes()+
  theme(text = element_text(size = 16, family = "serif"),
        legend.position = "top",
        panel.grid.major =element_line(),    #strip major gridlines
        panel.grid.minor = element_blank(),    #strip minor gridlines
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"))+
  scale_color_discrete(name = "")+
  guides(fill = guide_legend(override.aes = list(size = 5)))

outplot <- list("Dens.plot" = postDense,
                "Eye.plot" = eyeDense)
return(outplot)

}
