# Load deterministic results

T1_ab <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T1_alphabeta_deterministic.csv")[,2:7]
T1_Gr <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T1_Gr_deterministic.csv")[,2:7]
T1_h <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T1_lambda_deterministic.csv")[,2:6]

T2_ab <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/sens_T2_alphabeta_deterministic.csv")[,2:7]
T2_Gr <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/sens_T2_Gr_deterministic.csv")[,2:7]
T2_h <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/sens_T2_lambda_deterministic.csv")[,2:6]

T3_ab <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T3_alphabeta_deterministic.csv")[,2:7]
T3_Gr <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T3_Gr_deterministic.csv")[,2:7]
T3_h <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T3_lambda_deterministic.csv")[,2:6]

T4_ab <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T4_alphabeta_deterministic.csv")[,2:7]
T4_Gr <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T4_Gr_deterministic.csv")[,2:7]
T4_h <- read.csv("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Deterministic/Sens_T4_lambda_deterministic.csv")[,2:6]

# Save plot with 600 DPI resolution
tiff("SA_T1.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T1_ab$alpha,T1_ab$T1,col="red",pch=20,ylab=expression(T[1]^"*"),xlab=expression(alpha),main=expression("Sensitivity analysis for "~T[1]^"*"))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T1.2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T1_ab$beta,T1_ab$T1,col="red",pch=20,ylab=expression(T[1]^"*"),xlab=expression(beta),main=expression("Sensitivity analysis for "~T[1]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T1.3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T1_Gr$G,T1_Gr$T1,col="red",pch=20,ylab=expression(T[1]^"*"),xlab=expression(G),main=expression("Sensitivity analysis for "~T[1]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T1.4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T1_Gr$r,T1_Gr$T1,col="red",pch=20,ylab=expression(T[1]^"*"),xlab=expression(r),main=expression("Sensitivity analysis for "~T[1]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T1.5.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T1_h$lambda,T1_h$T1,col="red",pch=20,ylab=expression(T[1]^"*"),xlab=expression(lambda),main=expression("Sensitivity analysis for "~T[1]^"*"))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T2.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T2_ab$alpha,T2_ab$T2,col="red",pch=20,ylab=expression(T[2]^"*"),xlab=expression(alpha),main=expression("Sensitivity analysis for "~T[2]^"*"))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T2.2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T2_ab$beta,T2_ab$T2,col="red",pch=20,ylab=expression(T[2]^"*"),xlab=expression(beta),main=expression("Sensitivity analysis for "~T[2]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T2.3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T2_Gr$G,T2_Gr$T2,col="red",pch=20,ylab=expression(T[2]^"*"),xlab=expression(G),main=expression("Sensitivity analysis for "~T[2]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T2.4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T2_Gr$r,T2_Gr$T2,col="red",pch=20,ylab=expression(T[2]^"*"),xlab=expression(r),main=expression("Sensitivity analysis for "~T[2]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T2.5.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T2_h$lambda,T2_h$T2,col="red",pch=20,ylab=expression(T[2]^"*"),xlab=expression(lambda),main=expression("Sensitivity analysis for "~T[2]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T3.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T3_ab$alpha,T3_ab$T3,col="red",pch=20,ylab=expression(T[3]^"*"),xlab=expression(alpha),main=expression("Sensitivity analysis for "~T[3]^"*"))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T3.2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T3_ab$beta,T3_ab$T3,col="red",pch=20,ylab=expression(T[3]^"*"),xlab=expression(beta),main=expression("Sensitivity analysis for "~T[3]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T3.3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T3_Gr$G,T3_Gr$T3,col="red",pch=20,ylab=expression(T[3]^"*"),xlab=expression(G),main=expression("Sensitivity analysis for "~T[3]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T3.4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T3_Gr$r,T3_Gr$T3,col="red",pch=20,ylab=expression(T[3]^"*"),xlab=expression(r),main=expression("Sensitivity analysis for "~T[3]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T3.5.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T3_h$lambda,T3_h$T3,col="red",pch=20,ylab=expression(T[3]^"*"),xlab=expression(lambda),main=expression("Sensitivity analysis for "~T[3]^"*"))

# Close the device
dev.off()




# Save plot with 600 DPI resolution
tiff("SA_T4.1.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T4_ab$alpha,T4_ab$T4,col="red",pch=20,ylab=expression(T[4]^"*"),xlab=expression(alpha),main=expression("Sensitivity analysis for "~T[4]^"*"))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T4.2.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T4_ab$beta,T4_ab$T4,col="red",pch=20,ylab=expression(T[4]^"*"),xlab=expression(beta),main=expression("Sensitivity analysis for "~T[4]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T4.3.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T4_Gr$G,T4_Gr$T4,col="red",pch=20,ylab=expression(T[4]^"*"),xlab=expression(G),main=expression("Sensitivity analysis for "~T[4]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T4.4.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T4_Gr$r,T4_Gr$T4,col="red",pch=20,ylab=expression(T[4]^"*"),xlab=expression(r),main=expression("Sensitivity analysis for "~T[4]^"*"))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T4.5.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
plot(T4_h$lambda,T4_h$T4,col="red",pch=20,ylab=expression(T[4]^"*"),xlab=expression(lambda),main=expression("Sensitivity analysis for "~T[4]^"*"))

# Close the device
dev.off()





cor_T1_a <- cor(T1_ab$alpha,T1_ab$T1,method="spearman")
cor_T1_b <- cor(T1_ab$beta,T1_ab$T1,method="spearman")
cor_T1_G <- cor(T1_Gr$G,T1_Gr$T1,method="spearman")
cor_T1_r <- cor(T1_Gr$r,T1_Gr$T1,method="spearman")
cor_T1_h <- cor(T1_h$lambda,T1_h$T1,method="spearman")

cor_T2_a <- cor(T2_ab$alpha,T2_ab$T2,method="spearman")
cor_T2_b <- cor(T2_ab$beta,T2_ab$T2,method="spearman")
cor_T2_G <- cor(T2_Gr$G,T2_Gr$T2,method="spearman")
cor_T2_r <- cor(T2_Gr$r,T2_Gr$T2,method="spearman")
cor_T2_h <- cor(T2_h$lambda,T2_h$T2,method="spearman")

cor_T3_a <- cor(T3_ab$alpha,T3_ab$T3,method="spearman")
cor_T3_b <- cor(T3_ab$beta,T3_ab$T3,method="spearman")
cor_T3_G <- cor(T3_Gr$G,T3_Gr$T3,method="spearman")
cor_T3_r <- cor(T3_Gr$r,T3_Gr$T3,method="spearman")
cor_T3_h <- cor(T3_h$lambda,T3_h$T3,method="spearman")

cor_T4_a <- cor(T4_ab$alpha,T4_ab$T4,method="spearman")
cor_T4_b <- cor(T4_ab$beta,T4_ab$T4,method="spearman")
cor_T4_G <- cor(T4_Gr$G,T4_Gr$T4,method="spearman")
cor_T4_r <- cor(T4_Gr$r,T4_Gr$T4,method="spearman")
cor_T4_h <- cor(T4_h$lambda,T4_h$T4,method="spearman")


# incubation tornado plot
data <- data.frame(
  Variable = c("lambda", "G", "r", "alpha", "beta"),
  cor1 = c(cor_T1_h,cor_T1_G,cor_T1_r,cor_T1_a,cor_T1_b),
  cor2 = c(cor_T2_h,cor_T2_G,cor_T2_r,cor_T2_a,cor_T2_b),
  cor3 = c(cor_T3_h,cor_T3_G,cor_T3_r,cor_T3_a,cor_T3_b),
  cor4 = c(cor_T4_h,cor_T4_G,cor_T4_r,cor_T4_a,cor_T4_b)
)

# Order data by the absolute values of the Result column
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]
data1$Variable <- expression(lambda,G,r,alpha,beta)
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]
data2$Variable <- expression(lambda,G,r,alpha,beta)
data3 <- data[order(abs(data$cor3), decreasing = FALSE), ]
data3$Variable <- expression(lambda,G,r,alpha,eta)
data4 <- data[order(abs(data$cor4), decreasing = FALSE), ]
data4$Variable <- expression(lambda,G,r,alpha,eta)


# Save plot with 600 DPI resolution
tiff("SA_T1_det.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
barplot(
  data1$cor1, 
  names.arg = data1$Variable, 
  horiz = TRUE, 
  las = 1, # Rotate axis labels
  col = ifelse(data1$cor1 > 0, "blue", "red"), # Positive in blue, negative in red,
  xlim = c(-1, 1))

# Close the device
dev.off()



# Save plot with 600 DPI resolution
tiff("SA_T2_det.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
barplot(
  data2$cor2, 
  names.arg = data2$Variable, 
  horiz = TRUE, 
  las = 1, # Rotate axis labels
  col = ifelse(data2$cor2 > 0, "blue", "red"), # Positive in blue, negative in red,
  xlim = c(-1, 1))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T3_det.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
barplot(
  data3$cor3, 
  names.arg = data3$Variable, 
  horiz = TRUE, 
  las = 1, # Rotate axis labels
  col = ifelse(data3$cor3 > 0, "blue", "red"), # Positive in blue, negative in red,
  xlim = c(-1, 1))

# Close the device
dev.off()


# Save plot with 600 DPI resolution
tiff("SA_T4_det.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

# Increase font size globally
par(cex.lab = 2,   # Axis labels font size (1.5x default size)
    cex.main = 2,  # Title font size
    cex.axis = 1.2,  # Axis tick label font size
    mar = c(5, 5, 4, 2) + 0.1)  # Adjust margins (optional for better space)

# Produce plot
barplot(
  data4$cor4, 
  names.arg = data4$Variable, 
  horiz = TRUE, 
  las = 1, # Rotate axis labels
  col = ifelse(data4$cor4 > 0, "blue", "red"), # Positive in blue, negative in red,
  xlim = c(-1, 1))

# Close the device
dev.off()