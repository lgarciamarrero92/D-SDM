# D-SDM: Self-Adaptive Single-Diode Model Parameter Identification under Mismatching Conditions

## Overview
The **D-SDM Repository** provides the implementation of the methodology presented in the paper **"Self-Adaptive Single-Diode Model Parameter Identification under Mismatching Conditions"**. The objective is to robustly estimate the Single Diode Model (SDM) parameters from an experimental IV curve of a PV module operating under real conditions, potentially including mismatch effects.

## Installation Guide

### Prerequisites
We recommend using **[Anaconda](https://store.continuum.io/cshop/anaconda/)** as it bundles most of the required dependencies. If you haven't installed Anaconda yet, please download and install it from the official website before proceeding.

### Steps to Install D-SDM

#### 1. Create a New Anaconda Environment
Open a terminal (Unix) or Anaconda Prompt (Windows) and execute:
```bash
conda create --name d-sdm
```

#### 2. Activate the New Environment
```bash
conda activate d-sdm
```

#### 3. Install Required Packages

- **Pymoo** (for optimization algorithms):
  ```bash
  conda install pymoo==0.6.1
  ```

- **Pandas** (for data handling):
  ```bash
  conda install pandas
  ```

- **IPykernel** (for Jupyter integration):
  ```bash
  conda install ipykernel
  ```

#### 4. Add Jupyter Kernel for the Environment
```bash
python -m ipykernel install --user --name d-sdm --display-name "D-SDM PV"
```

#### 5. Launch Jupyter Notebook
To start a Jupyter Notebook and use the created kernel:
```bash
jupyter notebook
```
Then, select the **"D-SDM PV"** kernel from the available options.

#### 6. Deactivating the Environment
When done, deactivate the environment with:
```bash
conda deactivate
```

## Example Notebook
In the folder **notebooks**, there is an example notebook demonstrating the usage of D-SDM for parameter estimation.



