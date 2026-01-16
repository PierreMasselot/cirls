# setwd("C:/Users/PierreMasselot/Downloads")
data <- read.csv("data-raw/london.csv")
london <- dplyr::mutate(data, date = as.Date(date, format = "%d/%m/%Y"))
save(london, file = "data/london.rda")
