# _______________________________________________
# Author: Garrett Tetrault
# Date: 2/28/2019
# _______________________________________________
library(dplyr)
library(changepoint)
library(fpop)
data(neuroblastoma, package = "neuroblastoma")

# Desired id and chromosome to examine.
id <- "4"
chr <- "2"

# Filter data for only specified id and chromosome.
profile <- filter(neuroblastoma$profiles, 
                  profile.id == id, 
                  chromosome == chr)

# _______________________________________________
#                                     changepoint

# Calculate changepoints using changepoint library.
cpt_profile <- cpt.mean(profile$logratio,
                        penalty="Manual", 
                        pen.value="log(n)")

plot(cpt_profile, 
     type = "p", pch=20, col="black", # Data point params.
     cpt.col="green", cpt.width=2,    # Changepoint line params.
     main="Changepoint intervals using cpt \n w/ penalty 'Manual'.",
     xlab="Position", ylab="LogRatio")

# _______________________________________________
#                                            fpop

# Calculate changepoints using fpop library.
fpop_profile <- Fpop(profile$logratio, lambda=1)

# Get the starts and ends of each changepoint interval.
seg_end <- fpop_profile$t.est
seg_start <- c(1, seg_end[1:length(seg_end)-1])

# Get mean over each changepoint interval.
seg_mean <- vector(length=length(seg_start))
for(i in 1:length(seg_start)) {
  seg_mean[i] <- 
    profile$logratio[seg_start[i]:seg_end[i]] %>%
    mean()
}

# We use NA to avoid jumps in plot.
fpop_segs <- c(rep(seg_mean, times=(seg_end-seg_start)), NA)
fpop_segs[seg_start] <- NA

# Plot data and changepoint segments.
plot(profile$position, profile$logratio,
     type = "p", pch=20, col="black",
     main="Changepoint intervals using fpop.",
     xlab="Position", ylab="LogRatio")

par(new=TRUE)

lines(profile$position, fpop_segs, 
      col="green", lwd=2)

