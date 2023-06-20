## code to prepare `Nurses` dataset goes here

install.packages("haven")                             # Install haven package
library("haven")

Nurses <- read_sav("data-raw/Nurses.sav")
traum <- read.csv("data-raw/trauma.csv", sep=";")
usethis::use_data(trauma, overwrite = TRUE)
usethis::use_data(Nurses, overwrite = TRUE)
