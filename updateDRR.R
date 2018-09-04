#install.packages("devtools")
#install.packages("roxygen2")
#install.packages("usethis")
library("devtools")
library("roxygen2")

#setwd("H:/Projects/HPV/India")
#create("DRR")
#note that this working directory was created using "Create Project" dropdown from File (otherwise won't recognize that this is intended
#to be an R pacakge)

setwd("H:/Projects/HPV/India/DiffRelRisk")
document()

##open up github desktop
##click Repository > Show in Explorer (should pop up explorer)
##copy files into explorer from "H:/Projects/HPV/India/DRR"
##hit commit


##library("DRR")
##?DRR.CI
install_github("sampsonj74/DRR")
library("DiffRelRisk")
DRRCI(x1=20,x2=50,x3=100,n1=1000,n2=1000,n3=1000)


