#Christa Caggiano
#Fall/Winter 2016 
#Generates coverage plots for summed and averaged 
#peak data. Uses loess smoothing function with a 
#modest smoothing parameter to generate smooth 
#curves around peak center. Axes set for +/-500bp 
#and total signal per 10 million reads. 
#Code best run in R studio, line by line to check for 
#errors and in order to best customize. 

rm(list=ls()) #removes all variables in R's memory 
setwd("/Users/Christa/Documents/Rosbash Lab/Felipe Data Analysis") #sets working directory (directory coverage data is in)

average_brains = read.csv('average_clk_targets.csv', header = T, sep = ',') #loads files as dataframes from a CSV with headers
summed_brains = read.csv('summed_clk_targets.csv', header = T, sep = ',')

### average 200 genes###

#Note: graphs are generated using standard R plot function and are relatively simple
par(mai = c(1.5,  1.5, 0.5, 2.75)) #sets an invisible grid to plot on (helps to orient legend and make graph appear aesthetically pleasing)
plot(average_brains$X, average_brains$X2_200, pch=20, col='white',
     xlab='distance from clk summit (bp)', ylab='summed total signal/10 million reads', 
     main = "Top 200 CLK direct targets",
     xlim=c(-600, 600), ylim=c(0,1)) #sets X and Y axes
points(average_brains$X, average_brains$X6_200, pch=20, col='white') #handy for visualization- prints scatterplot (currently in white so it is invisible) color can be changed to see individual values 
points(average_brains$X, average_brains$X10_200, pch=20, col='white') #(not necessary for loess function)
points(average_brains$X, average_brains$X14_200, pch=20, col='white')
points(average_brains$X, average_brains$X18_200, pch=20, col='white')
points(average_brains$X, average_brains$X22_200, pch=20, col='white')
lines(loess.smooth(average_brains$X,average_brains$X2_200,span = 0.1), #uses loess smoothing function to do a local polynomial regression on data, creating a smooth curve
      col='red', lwd=1, lty=1) 
lines(loess.smooth(average_brains$X,average_brains$X6_200,span = 0.1), #span sets smoothing parameter, between 0 and 1, where 0 is no smoothing and 1 is straight lines 
      col='orange', lwd=1, lty=1) 
lines(loess.smooth(average_brains$X,average_brains$X10_200,span = 0.1),
      col='yellow', lwd=1, lty=1) 
lines(loess.smooth(average_brains$X,average_brains$X14_200,span = 0.1),
      col='green', lwd=1, lty=1) 
lines(loess.smooth(average_brains$X,average_brains$X18_200,span = 0.1),
      col='blue', lwd=1, lty=1) 
lines(loess.smooth(average_brains$X,average_brains$X22_200,span = 0.1),
      col='purple', lwd=1, lty=1) 
legend(x = 750, y = 0.7, xpd = TRUE, #prints legend (x and y sets relative to axes, cex sets relative to par grid)
       legend = c("zt2","zt6", "zt10", "zt14", "zt18", "zt22"), 
       lty = c(1,1,1,1,1,1), 
       col = c("red", "orange", "yellow", "green", "blue", "purple"), cex=0.75)

### summed 200 genes###
par(mai = c(1.5,  1.5, 0.5, 2.75)) 
plot(summed_brains$X, summed_brains$X2_200, pch=20, col='white',
     xlab='distance from clk summit (bp)', ylab='average total signal/10 million reads', 
     main = "Top 200 CLK direct targets",
     xlim=c(-600, 600), ylim=c(0,160)) 
points(summed_brains$X, summed_brains$X6_200, pch=20, col='white')
points(summed_brains$X, summed_brains$X10_200, pch=20, col='white')
points(summed_brains$X, summed_brains$X14_200, pch=20, col='white')
points(summed_brains$X, summed_brains$X18_200, pch=20, col='white')
points(summed_brains$X, summed_brains$X22_200, pch=20, col='white')
lines(loess.smooth(summed_brains$X,summed_brains$X2_200,span = 0.1),
      col='red', lwd=1, lty=1) 
lines(loess.smooth(summed_brains$X,summed_brains$X6_200,span = 0.1),
      col='orange', lwd=1, lty=1) 
lines(loess.smooth(summed_brains$X,summed_brains$X10_200,span = 0.1),
      col='yellow', lwd=1, lty=1) 
lines(loess.smooth(summed_brains$X,summed_brains$X14_200,span = 0.1),
      col='green', lwd=1, lty=1) 
lines(loess.smooth(summed_brains$X,summed_brains$X18_200,span = 0.1),
      col='blue', lwd=1, lty=1) 
lines(loess.smooth(summed_brains$X,summed_brains$X22_200,span = 0.1),
      col='purple', lwd=1, lty=1) 
legend(x = 750, y = 120, xpd = TRUE, 
       legend = c("zt2","zt6", "zt10", "zt14", "zt18", "zt22"), 
       lty = c(1,1,1,1,1,1), 
       col = c("red", "orange", "yellow", "green", "blue", "purple"), cex=0.75)

### average 200 genes###

#outputs a simple line graph, useful for showing degree of variance
par(mai = c(1.5,  1.5, 0.5, 2.75)) 
plot(average_brains$X, average_brains$X2_200, pch=20, col='white',
     xlab='distance from clk summit (bp)', ylab='average total signal/10 million reads', 
     main = "Top 200 CLK direct targets",
     xlim=c(-500, 500), ylim=c(0,1)) 
lines(average_brains$X,average_brains$X6_200, #line weight and type parameters can be changed to have fine control of lines 
      col='orange', lwd=1, lty=1)  
lines(average_brains$X,average_brains$X10_200,
      col='yellow', lwd=1, lty=1) 
lines(average_brains$X,average_brains$X14_200,
      col='green', lwd=1, lty=1) 
lines(average_brains$X,average_brains$X18_200,
      col='blue', lwd=1, lty=1) 
lines(average_brains$X,average_brains$X22_200,
      col='purple', lwd=1, lty=1) 
legend(x = 700, y = 0.7, xpd = TRUE, 
       legend = c("zt2","zt6", "zt10", "zt14", "zt18", "zt22"), 
       lty = c(1,1,1,1,1,1), 
       col = c("red", "orange", "yellow", "green", "blue", "purple"), cex=0.75)

### summed 200 genes###
par(mai = c(1.5,  1.5, 0.5, 2.75)) 
plot(summed_brains$X, summed_brains$X2_200, pch=20, col='white',
     xlab='distance from clk summit (bp)', ylab='summed total signal/10 million reads', 
     main = "Top 200 CLK direct targets",
     xlim=c(-500, 500), ylim=c(0,200)) 
lines(summed_brains$X,summed_brains$X6_200,
      col='orange', lwd=1, lty=1) 
lines(summed_brains$X,summed_brains$X10_200,
      col='yellow', lwd=1, lty=1) 
lines(summed_brains$X,summed_brains$X14_200,
      col='green', lwd=1, lty=1) 
lines(summed_brains$X,summed_brains$X18_200,
      col='blue', lwd=1, lty=1) 
lines(summed_brains$X,summed_brains$X22_200,
      col='purple', lwd=1, lty=1) 
legend(x = 700, y = 120, xpd = TRUE, 
       legend = c("zt2","zt6", "zt10", "zt14", "zt18", "zt22"), 
       lty = c(1,1,1,1,1,1), 
       col = c("red", "orange", "yellow", "green", "blue", "purple"), cex=0.75)

####################### if along with replicate averages, +/- Standard deviation data is availblae 
# the polygon feature acn be used to generate a shaded area indicating error in pooled plots 
# line weight and color may be adjusted to get clear visualization 

##with SD for averaging##

# par(mai = c(1.5,  1.5, 0.5, 2.75)) 
# plot(average_brains$x, average_brains$X2_50, pch=20, col='white',
#      xlab='distance from summit', ylab='summed signal', main = "Average Replicates zt2",
#      xlim=c(-300, 300), ylim=c(0,3500)) 
# points(average_brains$x, average_brains$X2_100, pch=20, col='white')
# points(average_brains$x, average_brains$X2_200, pch=20, col='white')
# polygon(c(average_brains$x, rev(average_brains$x)), 
#         c(average_brains$X2_50_upper, rev(average_brains$X2_50_lower)), 
#         col='salmon', border=NA)
# polygon(c(average_brains$x, rev(average_brains$x)), 
#         c(average_brains$X2_100_upper, rev(average_brains$X2_100_lower)), 
#         col='cadetblue3', border=NA)
# polygon(c(average_brains$x, rev(average_brains$x)), 
#         c(average_brains$X2_200_upper, rev(average_brains$X2_200_lower)), 
#         col='chartreuse2', border=NA)
# lines(loess.smooth(average_brains$x,average_brains$X2_50,span = 0.3),
#       col='red', lwd=3, lty=1) 
# lines(loess.smooth(average_brains$x,average_brains$X2_100,span = 0.3),
#       col='cadetblue4', lwd=3, lty=1) 
# lines(loess.smooth(average_brains$x,average_brains$X2_200,span = 0.3),
#       col='chartreuse4', lwd=3, lty=1) 
# legend(x = 300, y = 1200, xpd = TRUE, 
#        legend = c("50 genes","100 genes", "200 genes"), 
#        lty = c(1,1,1), 
#        col = c("red", "cadetblue4", "chartreuse4"),cex=0.75)


