########################################################
#                     Packages                         #
########################################################

# Function to test if packages that need to be loaded are already installed.
# If not installed, install using install.packages()

pkgTest <- function(x){
  if (suppressMessages(!require(x,character.only = TRUE, quietly = TRUE))){
    install.packages(x, repos = "http://cran.us.r-project.org")
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# Function to test if bioconductor packages that need to be loaded are already installed
# If not installed, install using biocLite()

bioPkgTest <- function(x){
  if (suppressMessages(!require(x,character.only = TRUE, quietly = TRUE))){
    suppressMessages(source("https://bioconductor.org/biocLite.R"))
    biocLite(x)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
