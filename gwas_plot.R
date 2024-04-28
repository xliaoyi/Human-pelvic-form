rm(list = ls())

setwd("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/UKB_Imaging_Genetics")
library(qqman)

### Manhattan plot

data <- read.delim("PLINK_LM_OUTPUT/.acetabular_diameter.glm.linear", header = TRUE)

# pdf("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/results/gwas_hip_height.pdf", width=17, height=11)
par(mar=c(5.5, 6.5, 4.1, 2.1))
manhattan(data, chr="X.CHROM", bp="POS", p="P", snp = 'ID',
          main="Acetabular diameter", col = c("tomato", "steelblue"),
          cex.axis=1.7, cex.lab=1.7, cex.main=1.7, cex.sub=1.7)
# dev.off()


### QQ plot

# Convert P-values to chi-squared statistics
data$chisq <- qchisq(1 - data$P_BOLT_LMM_INF, df=1)

# Calculate lambda
lambda <- median(data$chisq) / qchisq(0.5, df=1)

pdf("/Users/alexxu/Library/CloudStorage/Box-Box/Narasimhan_lab/hip_shape/results/lambda_hip_height.pdf", width=11, height=8.5)
par(mar=c(5.5, 6.5, 4.1, 2.1))
qq(data$P_BOLT_LMM_INF, main="GWAS QQ Plot", pch=19, col="black", cex=1, las=1,
   cex.lab = 1.7, cex.axis = 1.7)
abline(0, 1, col="red", lwd=3)
text(x=1, y=28, labels=paste("λ =", round(lambda, 3)), cex=2)
dev.off()