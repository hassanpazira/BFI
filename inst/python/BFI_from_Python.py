####-------------------
#### Python Instalation
####-------------------
https://blog.devgenius.io/calling-python-and-c-code-using-pybind-99ab7fefa685
### The steps to install the latest version of Python on 'MacOS' and 'Ubuntu':

## For MacOS:

# Open a terminal and install the package manager 'Homebrew' by running the following command:
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

# Install the latest version of Python (Python 3) with the following command:
brew install python@3

# Verify the installation by checking the Python version:
python3 --version

## For Ubuntu:

# Update package repositories and get latest package information
sudo apt-get update -y

# Install the latest version of Python (Python 3) using the following command in a terminal window:
sudo apt-get install python3

# Verify the installation by checking the Python version:
python3 --version

#### Installing PIP

### PIP (a Python package manager) helps us to install and use various packages/modules in Python programming.

## Install/Upgrade the latest version of PIP by running the following command:

# on MacOS
python3 -m pip install --user --upgrade pip

# on Ubuntu
sudo apt-get install python3-pip


#### Instalation of required moduls through PIP: 'pandas', 'numpy' and 'rpy2'

### Installing the 'pandas' module
pip3 install pandas
# On Ubuntu !
sudo python3 -m pip install wheel # and
sudo python3 -m pip install pandas

### Installing the 'numpy' module
pip3 install numpy

### Installing the 'rpy2' module
pip3 install rpy2
# On Ubuntu the line above may not work, if so use the following: 
sudo apt-get install python3-rpy2


### Now move on to work with 'rpy2' inside a Python script!
### First type the following command in the terminal window to go to the Python Environment:
python3

# Importing the rpy2 packages and subpackages
import pandas as pd               # importing the pandas packages
import numpy as np                # importing the numpy packages
import rpy2                       # importing the rpy2 packages
import rpy2.robjects as robjects  # importing the rpy2 subpackages

#### Installation of R packages from CRAN and GitHub:

# Importing the required functions
from rpy2.robjects.packages import importr, data   # importr() and data()
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import numpy2ri, pandas2ri
pandas2ri.activate()
numpy2ri.activate()

# Loading the packages preinstalled with R using importr()
utils = importr('utils')
base = importr('base')

## Installing required R packages from CRAN
utils.chooseCRANmirror(ind=1)  # select the first mirror in the list

package_names = ('stats', 'devtools')
utils.install_packages(StrVector(package_names))

# Loading the packages
stats = importr('stats')
devtools = importr('devtools')

## Installing the 'BFI' package from GitHub
devtools.install_github("hassanpazira/BFI", force = True)

# Loading the 'BFI' package
BFI = importr('BFI')

### Examples:
#np.random.seed(1)

## Data Simulation for center 1 from Gaussian distribution:
theta = np.array([1, 2, 3, 4, 0.75]) # model parameters: intercept is theta[0], sigma2 is theta[4]
p = 3     # number of regression parameters without intercept
n1 = 30   # sample size of center 1
X1 = np.random.normal(size=(n1, p)) # X1.shape
#XX = rpy2.robjects.DataFrame(X1) #!
X1 = pd.DataFrame(X1, columns=['X1', 'X2', 'X3']) # list(X1) or X1.columns
mu1 = theta[0] + np.dot(X1, np.delete(theta, [0, 4]))   # in gaussian \eta = \mu 
X1 = pandas2ri.py2rpy(X1)   # == base.as_data_frame(X1)
y1 = np.random.normal(loc=mu1, scale=np.sqrt(theta[4]))

## Data Simulation for center 2 from Gaussian distribution:
n2 = 50   # sample size of center 2
X2 = np.random.normal(size=(n2, p))
X2 = pd.DataFrame(X2, columns=['X1', 'X2', 'X3'])
mu2 = theta[0] + np.dot(X2, np.delete(theta, [0, 4]))   # in gaussian \eta = \mu 
X2 = pandas2ri.py2rpy(X2)
y2 = np.random.normal(loc=mu2, scale=np.sqrt(theta[4]))

# MAP estimates at center 1
Lambda1 = BFI.inv_prior_cov(X1, 0.01, 'gaussian')
fit1 = BFI.MAP_estimation(y1, X1, 'gaussian', Lambda1)
print(fit1)
theta_hat1 = fit1.rx2("theta_hat") # intercept and coefficient estimates
A_hat1 = fit1.rx2("A_hat")     # curvature matrix
summary_1 = BFI.summary_bfi(fit1, cur_mat = True)

# MAP estimates at center 2
Lambda2 = BFI.inv_prior_cov(X2, 0.01, 'gaussian')
fit2 = BFI.MAP_estimation(y2, X2, 'gaussian', Lambda2)
theta_hat2 = fit2.rx2("theta_hat")
A_hat2 = fit2.rx2("A_hat")
summary_2 = BFI.summary_bfi(fit2, cur_mat = True)

# BFI at central center
A_hats = base.list(A_hat1, A_hat2)
theta_hats = base.list(theta_hat1, theta_hat2)
Lambda1 = pd.DataFrame(Lambda1, index=fit1.rx2("names"), columns=fit1.rx2("names"))
Lambda2 = pd.DataFrame(Lambda2, index=fit2.rx2("names"), columns=fit2.rx2("names"))
Lambdas = base.list(Lambda1, Lambda2)
fit_bfi = BFI.bfi(theta_hats, A_hats, Lambdas)
print(fit_bfi)
summary_bfi = BFI.summary_bfi(fit_bfi, cur_mat = True)


numpy2ri.deactivate()
pandas2ri.deactivate()

####################

Nurses = data(BFI).fetch('Nurses')['Nurses']  # is equivalent to BFI::Nurses in R
trauma = data(BFI).fetch('trauma')['trauma']  # is equivalent to BFI::trauma in R

data0 = robjects.r("""
set.seed(1)
xx <- rnorm(5)
""")
robjects.r["xx"]
########

