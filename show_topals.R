
show_topals = function(fit, emphasize_sd=FALSE, hue='red') {
 
 df_grouped = data.frame(
                L = head( fit$age_group_bounds, -1),
                U = tail( fit$age_group_bounds, -1),
                N = fit$N,
                D = fit$D
               ) %>%
          mutate(logmx_obs = log(D/N))
              
 
 df_single  = data.frame(
                  age=  seq(fit$std) - .50,  # 0.5, 1.5, ...
                  std = fit$std,
                  logmx_true = ITA$logmx[1:100],
                  logmx_fit  = fit$logm
               )

  this_plot =
  ggplot(data = df_single, aes(x=age,y=logmx_true)) +
      labs(x='Age',y='Log Mortality Rate',
           title='Italy Females 1980',
           subtitle = paste(sum(fit$D),'deaths to',round(sum(fit$N)),'women')) +
      scale_x_continuous(breaks=c(0,1,seq(5,100,5)),minor_breaks = NULL) +
      scale_y_continuous(limits=c(-10,0),breaks=seq(-10,0,2),minor_breaks = NULL) +
      theme_bw()

 
 if (!emphasize_sd) {

    
 this_plot = this_plot + 
      geom_line(aes(x=age,y=std), color='black', lwd=0.5) +
      geom_line(aes(x=age,y=logmx_fit), color=hue, lwd=3, alpha=.40) +
      geom_segment(data=df_grouped,aes(x=L,xend=U,
                                       y=logmx_obs,
                                       yend=logmx_obs),
                   color=hue,lwd=1.5, alpha=.90) +
      geom_point(size=0.60, alpha=.70) 
 

 } else {

    sd = sqrt( diag( fit$B %*% fit$covar %*% t(fit$B)))
    factor = qnorm(.90)
    df_errors = data.frame(age = seq(fit$std) - 0.50,
                           L = fit$logm - factor * sd,
                           H = fit$logm + factor * sd)
    
    this_plot = this_plot +
                  geom_segment(data=df_errors,
                               aes(x=age, xend=age, y=L, yend=H), lwd=0.5, color=hue,
                               inherit.aes = FALSE)
    
  }

  print(this_plot)
} # show_topals  

