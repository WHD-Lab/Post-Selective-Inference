# Post-Selective-Inference
Lab repository for post selective inference methods and corresponding code

## Folders:

### Analysis Methods:
Stores code corresponding to the specific methodology. Currently, we only have Basic LASSO.

- #### Basic LASSO
In this folder, it stores the original copy of code for randomized WCLS with LASSO. There is a folder called "Updated Code"
within this folder, which stores updated or edited code. 

### Data Simulation & Pesudo Outcome
This folder stores function corresponding to data simulation and pesudo-outcome generation.

## Tutorial:

### Download Python code and setup virtual environment

- ### Windows system

1. Install pipx

```
python -m pip install --user pipx
```

2. In command line type: ``pipx ensurepath``

```
python -m pipx ensurepath
```

3. open new terminal (so that path changes take effect)

4. clone [Snigdha's repository](https://github.com/snigdhagit/selective-inference/tree/refactor_names), **make sure clone the refactor_names branch**. In the command, type:

``git clone -b refactor_names https://github.com/snigdhagit/selective-inference.git``

5. cd into selective-inference repo. 

6. Init and populate submodules: ``git submodule update --init``

7. Don’t use the latest version of Python it will not work. Use the older version. Create a new python3 virtual environment: ``virtualenv env3 -p python3.9`` **Make sure you have virtualenv installed, and it's installed in python3.9 environment**

8.	Active the environment: ``env3/Scripts/activate`` If it doesn’t work try ``env3\Scripts\activate``

9. In the command type: ``pip install git+https://github.com/regreg/regreg.git``. This step will tend to fail, and it will ask you to install some uninstalled packages like Cython, numpy. Make sure you use ``pip install numpy==1.22.4`` for numpy. Otherwise, you will get error message “module 'numpy' has no attribute 'asscalar'” when run the code.If you encounter error: ```python setup.py bdist_wheel did not run successfully```. Run below code

```
pip install wheel
python setup.py bdist_wheel 
```

Sometimes you need to install c++

10. In the command type: ``pip install -e . -r requirements.txt -r dev-requirements.txt``.This will gives a error “The 'sklearn' PyPI package is deprecated, use 'scikit-learn' rather than 'sklearn' for pip commands.” This can be fixed by opening the file called “requirements” in downloaded selective-inference folder. Find sklearn replacing it with scikit-learn

- ### Mac

1. Install [homebrew](https://brew.sh/)

2. Install pipx
    1. `brew install pipx`
    2. Then `pipx ensurepath`
    3. open new terminal, so the path changes take effect
  
3. Use pipx to install virtualenv: ``pipx install virtualenv``

4.  clone [Snigdha's repository](https://github.com/snigdhagit/selective-inference/tree/refactor_names), **make sure clone the refactor_names branch**. In the command, type:

``git clone -b refactor_names https://github.com/snigdhagit/selective-inference.git``

5. cd into selective-inference repo. 

6. Init and populate submodules: ``git submodule update --init``

7. Create a new python3 virtual environment: ``virtualenv env3 -p python3.9`` **Make sure you have virtualenv installed, and it's installed in python3.9 environment**

8. Active the environment: ``source env3/bin/activate``

9. In the command type: ``pip install git+https://github.com/regreg/regreg.git``. This step will tend to fail, and it will ask you to install some uninstalled packages like Cython, numpy. Make sure you use ``pip install numpy==1.22.4`` for numpy. Otherwise, you will get error message “module 'numpy' has no attribute 'asscalar'” when run the code.

10. In the command type: ``pip install -e . -r requirements.txt -r dev-requirements.txt``.This might give a error “The 'sklearn' PyPI package is deprecated, use 'scikit-learn' rather than 'sklearn' for pip commands.” This can be fixed by opening the file called “requirements” in downloaded selective-inference folder. Find sklearn replacing it with scikit-learn

**Some errors will show up if you use Mac system to run the dependent R package ```reticulate```. Here is the solution**

- If you get error that system already use default Python version instead of the version in ```env3```. 

    - unload the R package ```reticulate``` and run below code in R console

```
Sys.setenv(RETICULATE_PYTHON = "path to selective-inference/env3/bin/python")
Sys.setenv(RETICULATE_PYTHON_ENV = "path to selective-inference/env3")
library(reticulate)
use_virtualenv("path to selective-inference/env3")
```

- If you get error about conda, use below code to specify which conda environment to use

```
use_condaenv(condaenv = 'env3', conda = "path to anaconda3/bin/conda", required = T)
```

### Simulate data and test the code

**Make sure you set R working directory to ```setwd("/path/to/Post-Selective-Inference-Lab")```. Then source file ``source("./Analysis Methods/Basic LASSO/Updated Code/UI_function_LASSO.R")`` before run any of the following code**

- #### Simulate data

The data can be easily simulated using the code

```
data = generate_dataset(N = 1000, T = 50, P = 50, sigma_residual = 5, sigma_randint = 5, 
                        main_rand = 5.5, rho = 0.7,
                        beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4, 
                        beta = matrix(c(-0.2, 0.8, 0.3, 0.7, 0.3),ncol = 1),
                        theta1 = 0.8)
```
Then you should decide the control variables ($H_t$) and potential moderators ($S_t$). For example, you can write following code to randomly pick the first 50 variables are controls and the first 25 of them are potential moderators.

```
Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))
```

- #### Analysis

Then you can plug in the above simulated data to do analysis using the following code. **Make sure to replace ```virtualenv_path``` with our own path to env3**

```
# use CV LASSO to generate pseudo outcome
UI_return = DR_WCLS_LASSO(data = data, fold = 5, ID = "id", time = "decision_point", 
                          Ht = Ht, St = St, At = "action", outcome = "outcome", method_pesu = "CVLASSO", 
                          lam = NULL, noise_scale = NULL, splitrat = 0.8, 
                          virtualenv_path = "path to selective-inference folder/env3",
                          beta =  c(-0.2, 0.8, 0.3, 0.7, 0.3, rep(0, 21)), level = 0.9, core_num = 3)
                          
# use Random Forest to generate pseudo outcome                          
UI_return_RF =  DR_WCLS_LASSO(data = data, fold = 5, ID = "id", time = "decision_point", 
                              Ht = Ht, St = St, At = "action", 
                              outcome = "outcome", method_pesu = "RandomForest", 
                              lam = NULL, noise_scale = NULL, splitrat = 0.8, 
                              virtualenv_path = "path to selective-inference folder/env3",
                              beta =  c(-0.2, 0.8, 0.3, 0.7, 0.3, rep(0, 21)), level = 0.9, core_num = 3)
```

# Important Todo:
1. remove the dependence on Python code
2. Makes code can handle missing data
3. Adapt code for binary outcomes
