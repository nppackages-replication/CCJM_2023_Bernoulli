library(lpcde)
library("tikzDevice")

# load data
data <- read.csv("bikeSharing.csv")

result1 <- lpcde(data$atemp2, data$count2, x=0, bw_type="imse-rot", y_grid = seq(0, 2.5, length.out=20))
plot(result1)

result2 <- lpcde(data$atemp2, data$count2, x=2.23, bw_type="imse-rot", y_grid = seq(0, 4, length.out=20))
plot(result2)

result3 <- lpcde(data$atemp2, data$count2, x=3.12878, bw_type="imse-rot", y_grid = seq(0, 3.5, length.out=20))
plot(result3)

# Put three estimates in one matrix
temp1 <- cbind(result1$Estimate[, 1], result1$Estimate[, 3], result1$Estimate[, 4], result1$Estimate[, 6])
temp2 <- cbind(result2$Estimate[, 1], result2$Estimate[, 3], result2$Estimate[, 4], result2$Estimate[, 6])
temp3 <- cbind(result3$Estimate[, 1], result3$Estimate[, 3], result3$Estimate[, 4], result3$Estimate[, 6])

# Plot
tikz(paste("fig1", ".tex", sep=""),standAlone=TRUE, height=5, width=5)
plot(temp1[, 1], temp1[, 2], ylim=c(0, 1.6), xlim=c(0, 4),
     col="#3C93C2", lwd=3,
     type="l", ylab="$\\widehat{f}(y|\\mathbf{x})$", xlab="$y_i$: number of rentals",
     xaxt="n")

lines(temp2[, 1], temp2[, 2],
      col="#F08F6E", lwd=3, lty=2
      )

lines(temp3[, 1], temp3[, 2],
      col="#AB1866", lwd=3, lty=4
)
axis(side=1, at=seq(0, 1000, 200) / sd(data$count),
     labels = seq(0, 1000, 200))

lines(c(2.5-0.6, 2.8-0.6), c(1.3, 1.3), col="#3C93C2", lwd=3)
text(3.2-0.6, 1.3, "0", pos=2, col="#3C93C2")

lines(c(2.5-0.6, 2.8-0.6), c(1.2, 1.2), col="#F08F6E", lwd=3, lty=2)
text(3.2-0.6, 1.2, "25", pos=2, col="#F08F6E")

lines(c(2.5-0.6, 2.8-0.6), c(1.1, 1.1), col="#AB1866", lwd=3, lty=4)
text(3.2-0.6, 1.1, "35", pos=2, col="#AB1866")

text(3.9-0.1, 1.45, "$\\mathbf{x}_i$: feels-like temperature ($^\\circ C$)", pos=2)

dev.off()

# Plot
tikz(paste("fig2", ".tex", sep=""),standAlone=TRUE, height=5, width=5)
plot(temp3[, 1], temp3[, 2], ylim=c(0, 1.6), xlim=c(0, 4),
     col="#AB1866", lwd=3, lty=4,
     type="l", ylab="$\\widehat{f}(y|\\mathbf{x})$", xlab="$y_i$: number of rentals",
     xaxt="n")

axis(side=1, at=seq(0, 1000, 200) / sd(data$count),
     labels = seq(0, 1000, 200))

lines(c(2.5-0.6, 2.8-0.6), c(1.1, 1.1), col="#AB1866", lwd=3, lty=4)
text(3.2-0.6, 1.1, "35", pos=2, col="#AB1866")

polygon(c(c(2.5-0.6, 2.8-0.6), rev(c(2.5-0.6, 2.8-0.6))), c(1.13, 1.13, 1.07, 1.07), col = adjustcolor("#AB1866", alpha.f=0.2),
      border=NA)
text(2.5, 1.1, "(with confidence band)", pos=4, col="#AB1866")

text(3.9-0.1, 1.45, "$\\mathbf{x}_i$: feels-like temperature ($^\\circ C$)", pos=2)

polygon(c(temp3[, 1], rev(temp3[, 1])), c(temp3[, 3] + 2.1 * temp3[, 4], rev(temp3[, 3] - 2.1 * temp3[, 4])),
       col = adjustcolor("#AB1866", alpha.f=0.2),
       border=NA)
dev.off()
