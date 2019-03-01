library(dplyr)
library(changepoint)
library(fpop)
data(neuroblastoma, package = "neuroblastoma")

# Desired id and chromosome to examine.
id <- "4"
chromo <- "2"

# Filter data for only specified id and chromosome.
profile <- filter(neuroblastoma$profiles, 
                  profile.id == id, chromosome == chromo)

# _______________________________________________
#                                     changepoint

# Calculate changepoints using changepoint library.
cpt_profile <- cpt.mean(profile$logratio, method = "PELT", 
                        penalty="Manual", pen.value="log(n)/2")

plot(cpt_profile, 
     type = "p", pch=20, col="black", # Data point params.
     cpt.col="green", cpt.width=2,    # Changepoint line params.
     main="Changepoint intervals using cpt \n w/ penalty 'Manual'.",
     xlab="Position", ylab="LogRatio")

# _______________________________________________
#                                            fpop

# Calculate changepoints using fpop library.
fpop_profile <- Fpop(profile$logratio, lambda=1)
print(fpop_profile)

# Get the starts and ends of each changepoint interval.
fpop_end <- fpop_profile$t.est
fpop_start <- c(1, fpop_end[1:length(fpop_end)-1])

# Create plottable data of the changepoint intervals.
fpop_segs <- vector(length=length(profile$position))
for(i in 1:length(fpop_start)) {
  start <- fpop_start[i]
  end <- fpop_end[i]
  
  # Set segement value to mean over the changepoint interval.
  seg_mean <- mean(profile$logratio[start:end])
  seg_length <- length(fpop_segs[start:end])
  
  fpop_segs[start:end] <- rep.int(seg_mean, times=seg_length)
  fpop_segs[start] <- NA # We use NA to avoid jumps in plot.
}

# Plot data and changepoint segments.
plot(profile$position, profile$logratio,
     type = "p", pch=20, col="black",
     main="Changepoint intervals using fpop.",
     xlab="Position", ylab="LogRatio")

par(new=TRUE)

lines(profile$position, fpop_segs, 
      col="green", lwd=2)

