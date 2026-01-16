################################################################################
# Tidy the fgl dataset

# Load
load("data-raw/fgl.rda")

# Fill in zeros
ind <- 2:9
fgl[, ind][fgl[,ind] == 0] <- min(fgl[,ind][fgl[,ind] > 0]) / 2

# Rescale to unit sum
fgl[,ind] <- fgl[,ind] / rowSums(fgl[,ind])

# Export
save(fgl, file = "data/fgl.rda")
