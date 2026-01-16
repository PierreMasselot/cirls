################################################################################
# Build dataset warming

# Read
warming <- read.csv("data-raw/warming.csv")

# Add decade
warming$decade <- factor(10 * floor(warming$year / 10))
warming <- warming[, c("year", "decade", "anomaly")]

# Export
save(warming, file = "data/warming.rda")
