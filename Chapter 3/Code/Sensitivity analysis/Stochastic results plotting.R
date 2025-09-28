T2_ab_s <- readRDS("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Stochastic/sens_T2_alphabeta_stochastic.rds")
T2_Gr_s <- readRDS("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Stochastic/sens_T2_Gr_stochastic.rds")
T2_h_s <- readRDS("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Stochastic/sens_T2_lambda_stochastic.rds")

T4_ab_s <- readRDS("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Stochastic/sens_T4_alphabeta_stochastic.rds")
T4_Gr_s <- readRDS("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Stochastic/sens_T4_Gr_stochastic.rds")
T4_h_s <- readRDS("~/PhD-work/Chapter 3/Data/Sensitivity analysis/Stochastic/sens_T4_lambda_stochastic.rds")

cor_T2_a <- cor(T2_ab_s$alpha,T2_ab_s$time,method="spearman")
cor_T2_b <- cor(T2_ab_s$beta,T2_ab_s$time,method="spearman")
cor_T2_G <- cor(T2_Gr_s$G,T2_Gr_s$time,method="spearman")
cor_T2_r <- cor(T2_Gr_s$r,T2_Gr_s$time,method="spearman")
cor_T2_h <- cor(T2_h_s$lambda,T2_h_s$time,method="spearman")

cor_T4_a <- cor(T4_ab_s$alpha,T4_ab_s$time,method="spearman")
cor_T4_b <- cor(T4_ab_s$beta,T4_ab_s$time,method="spearman")
cor_T4_G <- cor(T4_Gr_s$G,T4_Gr_s$time,method="spearman")
cor_T4_r <- cor(T4_Gr_s$r,T4_Gr_s$time,method="spearman")
cor_T4_h <- cor(T4_h_s$lambda,T4_h_s$time,method="spearman")

# incubation tornado plot
data <- data.frame(
  Variable = c("lambda", "G", "r", "alpha", "beta"),
  cor1 = c(cor_T2_h,cor_T2_G,cor_T2_r,cor_T2_a,cor_T2_b),
  cor2 = c(cor_T4_h,cor_T4_G,cor_T4_r,cor_T4_a,cor_T4_b)
)

# Order data by the absolute values of the Result column
data1 <- data[order(abs(data$cor1), decreasing = FALSE), ]
data1$Variable <- expression(lambda,G,r,alpha,beta)
data2 <- data[order(abs(data$cor2), decreasing = FALSE), ]
data2$Variable <- expression(lambda,G,r,alpha,beta)


# Save plot with 600 DPI resolution
tiff("SA_T2_sto.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

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
tiff("SA_T4_sto.tiff", width = 8, height = 6, units = "in", res = 600, compression = "lzw")

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
