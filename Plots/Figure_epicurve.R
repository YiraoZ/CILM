#Figure 5 in Appendix - epicurves under seven scenarios
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
#set the file path
#plot epicurve
plot_epicurve<-function(){
  spatials = c("homogeneous", "heterogeneous_small_K=3", "heterogeneous_small_K=5","heterogeneous_small_K=8","heterogeneous_large_K=3","heterogeneous_large_K=5", "heterogeneous_large_K=8")
  titles = c("CSR", "Low-variance Clustered, K=3", "Low-variance Clustered, K=5", "Low-variance Clustered, K=8", "High-variance Clustered, K=3", "High-variance Clustered, K=5", "High-variance Clustered, K=8")
  filepath_prefix = "./simulation_data/"#prefix of the epidata filepath
  plot_list <-vector('list', 7)
  for(k in 1:7){
    df<-NULL
    time <- as.data.frame(0:20)
    colnames(time)<-c("time_of_infected")
    for(n in 1:10){
      filepath_data <- paste(filepath_prefix,spatials[k],"/",spatials[k],"_",toString(n),".csv",sep = "")
      data<-read.csv(filepath_data)
      daily_cases <- time%>%left_join(data %>%
                                        filter(time_of_infected != -1) %>%
                                        group_by(time_of_infected) %>%
                                        summarize(Cases = n()), by = c("time_of_infected"="time_of_infected"))%>%
        replace_na(list(Cases = 0))
      temp_df<-data.frame(x=daily_cases$time_of_infected, y=daily_cases$Cases, col=rep(n:n, each=length(daily_cases$time_of_infected)))
      df<-rbind(df, temp_df)
    }
    epicurve<-ggplot(df,aes(x=x,y=y,group=col,colour=factor(col))) + geom_line(size = 1)+scale_colour_brewer(palette = "Dark2")+
      ggtitle(titles[k])+
      theme(
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA))+labs(x = "Time", y = "Incidence")+guides(color = FALSE)
    plot_list[[k]]<-epicurve
    }
  return(plot_list)
}
myplots<-plot_epicurve()
pdf("./Plots/epicurve_plot.pdf", width = 10, height = 16)
grid.arrange(grobs = myplots, nrow = 4, ncol = 2)
dev.off()