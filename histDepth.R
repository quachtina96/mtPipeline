#!/bin/R

#this R script is used in mtPipeline's simplePipe.py in order to generate histograms of the depth of coverage 
# for pre and post remapping of the chrM to the RCRS (Revised Cambridge Reference Sequence). 

#It reads in depth data from a file called "depth.txt".
#It creates a pdf of the histogram in "depth.pdf".

#simplepipe.py renames this depth.pdf to make it specific to the sample.


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


