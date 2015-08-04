#!/bin/R
library(ggplot2) #the import statement. 2nd argument indicates location of libarry

coverageData=read.table("depth.txt", 
	header = F, 
	sep = "\t", 
	col.names = c("id","baseNumber","depth"))
plot=qplot(coverageData$depth, 
	geom="histogram", 
	binwidth = 5, 
	main = "Depth of Coverage Histogram", 
	xlab = "Read Depth", 
	fill=I("blue"), 
	col=I("black"),
	alpha=I(".2"),
	xlim = c(0,150)) 
ggsave(plot,file="depth.pdf")


