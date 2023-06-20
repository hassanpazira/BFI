## code to prepare `Nurses` dataset goes here

install.packages("haven")    # Install haven package
library("haven")
install.packages("labelled")
library(labelled)
install.packages("foreign")
library(foreign)

# Nurses <- suppressWarnings(read.spss("data-raw/Nurses.sav",
#                                      use.value.labels = TRUE,
#                                      to.data.frame = TRUE)) # by 'foreign'
Nurses <- read_sav("data-raw/Nurses.sav") # by 'haven'
traum <- read.csv("data-raw/trauma.csv", sep=";")
usethis::use_data(trauma, overwrite = TRUE)
usethis::use_data(Nurses, overwrite = TRUE)
