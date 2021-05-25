### Using GANs to Evaluate RDD Estimators
##### Authors: Guido Imbens and Evan Munro

This repository contains code that can be used to replicate the analysis in the authors' paper.

The code is in both R and Python 3, both of which must be installed. Before replicating
any of the analysis, complete the following steps.

**1. Clone our Github repository and change working directory to the repository**
```
git clone https://github.com/evanmunro/rdd-gan/
cd rdd-gan
```

**2. Install the python requirements using pip**
```
pip install -r requirements.txt
```
This will install outside packages in Python, as well as our own Python package, `wgan`, for estimating WGANs on economic datasets, which is hosted at https://github.com/gsbDBI/ds-wgan/.

**3. Install the required R packages**

The following packages are required. They can be installed through `install.packages`, or `devtools`. 

```
install.packages(c("arrow", "rddtools", "kableExtra", "mgcv"))
devtools::install_github("swager/optrdd")
devtools::install_github("kolesarm/RDHonest")
```
While `optrdd` works with the quadprog optimizer, it runs much faster with the `mosek` optimizer ([installation instructions](https://github.com/swager/optrdd)).


**4. Download the simulated data from Google Drive**

The simulated datasets are too large to be hosted on Github. The empty folder in `data/generated/`
must be filled with the contents of the following Google Drive [folder]() before running any replication steps.

### A. Generating Data

This step re-trains the GANs and generates new simulated data. This step is computationally intensive and is best run on a GPU. If you are just interested in replicating the Monte Carlo simulation results, then skip this step.

```
python generation/estimate_gans.py
```

### B. Running Monte Carlo Simulations

```
Rscript simulation/slurm_simulation.R
```

The Monte Carlo Simulations are run on Stanford's Sherlock computing cluster.  
