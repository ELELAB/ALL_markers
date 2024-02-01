### This script contains functions used to do MDS on gene expression data

MDSPlot <- function(my.data, my.group, my.labels, my.cols) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  ggplot(res, aes(x=M1, y=M2)) + 
    geom_point(aes(fill = my.group, colour = my.group), shape=21, size = 5, stroke = 0.1) + 
    #geom_text(aes(x=M1,y=M2, label= my.labels)) + 
    guides(fill = guide_legend(override.aes = list(shape = 22))) + 
    scale_fill_manual(values=my.cols) + 
    scale_colour_manual(values=rep("white", length(my.cols))) + 
    theme_bw() + 
    theme(legend.title=element_blank()) + 
    theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + 
    theme(legend.position = "top") + 
    theme(axis.text=element_text(size=16, face="bold")) + 
    theme(axis.text = element_text(colour = "black"))
}