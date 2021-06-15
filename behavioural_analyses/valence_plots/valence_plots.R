    #this code to plot out acceptance rates
    #N.Garrett, sept 2018
    
    #Clear workspace
    rm(list = ls())
    
    #Get rid of exisiting plots
    dev.off()
    
    #Load packages - may need to install if not already
    library(ggplot2)
    library(data.table)
    library(cowplot)
    library(Hmisc)
    library(plyr)
    
    dat = fread("data/valence_means_fMRI.csv")
    #dat = fread("data/valence_means_behaviour.csv")

    dat$neg_combined =  (dat$desirable_neg + dat$undesirable_neg)/2
    dat$pos_combined =  (dat$desirable_pos + dat$undesirable_pos)/2
    
    #for block A, 
    #convert to long format for plotting so have all %in one column
    data_individual = melt(dat, measure.vars = c("desirable_neg", "desirable_pos", "undesirable_neg", "undesirable_pos", "neg_combined", "pos_combined"), 
                  variable.name = c("condition"), 
                  value.name = c("percent_consistent"))
    
    #code the options as factors
    data_individual$condition = as.factor(data_individual$condition)
    
    data_individual$valence = as.factor("negative")
    data_individual[condition=="desirable_pos" | condition=="undesirable_pos"| condition=="pos_combined", ]$valence = "positive"
    
    data_individual$state = as.factor("desirable")
    data_individual[condition=="undesirable_neg" | condition=="undesirable_pos", ]$state = "undesirable"
    data_individual[condition=="pos_combined" | condition=="neg_combined", ]$state = "NA"
    
    #now calculate the summary stats - means and sem - for each option for each block
    
    #note that for the initial variable (as this table hasn't been created yet), 
    #you have to do a long winded code to get the varaible name
    data_summary = data_individual[, .(mean = mean(percent_consistent)), by=c("valence", "state")]
    data_summary$sem = data_individual[, .(sd(percent_consistent)/sqrt(.N)), by=c("valence", "state")][,V1]
    
    # #order by block then by rank
    # setkey(data_summary, condition)
    
    title_text = "valence effect"
    
    #plot
    p1 <- ggplot(data = data_summary[state=="NA", ], aes(x = valence, y = mean, fill=valence)) +
      geom_bar(stat = "identity", width = .5, position="dodge") + coord_cartesian(ylim=c(0.0,1))+
      geom_errorbar(data = data_summary[state=="NA", ], aes(ymin = mean-sem, ymax = mean+sem), position = position_dodge(width = .5), width=0.25) + 
      geom_point(data = data_individual[state=="NA", ], aes(x = valence, y = percent_consistent, colour = as.factor(valence)), stroke = 2, colour = "black", alpha = 0.2, size = 1.5, pch=21, position = position_jitterdodge(dodge.width=0.5)) +
      #scale_color_manual(values = c("positive" = "blue", "negative" = "red"))+
      scale_fill_manual(values = c("positive" = "springgreen", "negative" = "salmon"))+
      ylab("consistency score") +
      xlab("outcome (previous trial)") +
      theme_cowplot()+
      #theme(legend.position="none")+
      #ggtitle(title_text)+
      #scale_fill_brewer(palette="Set1", direction=1)+
      #scale_y_continuous(limits = c(0.5, 1))+
      theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), plot.title = element_text(size = 10, face = "bold"))
    
    p1

    p2 <- ggplot(data = data_summary[state=="NA", ], aes(x = valence, y = mean, fill=valence)) +
      geom_bar(stat = "identity", width = .5, position="dodge") + coord_cartesian(ylim=c(0.25,1))+
      geom_errorbar(data = data_summary[state=="NA", ], aes(ymin = mean-sem, ymax = mean+sem), position = position_dodge(width = .5), width=0.25) + 
      geom_point(data = data_individual[state=="NA", ], aes(x = valence, y = percent_consistent, colour = as.factor(valence)), stroke = 2, colour = "black", alpha = 0.2, size = 1.5, pch=21, position = position_jitterdodge(dodge.width=0.5)) +
      scale_fill_manual(values = c("positive" = "blue", "negative" = "red"))+
      #scale_fill_manual(values = c("positive" = "green", "negative" = "orange"))+
      ylab("consistency score") +
      xlab("state (trial t-1)") +
      theme_cowplot()+
      #theme(legend.position="none")+
      ggtitle(title_text)+
      #scale_fill_brewer(palette="Set1", direction=1)+
      #scale_y_continuous(limits = c(0.5, 1))+
      theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), plot.title = element_text(size = 10, face = "bold"))
    
    p2 
    