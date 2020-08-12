# Code that looks like this should be typed into the command line

```
type me into the command line pls
```

# Make sure Anaconda3 is installed and working, the following should list all your created envrionments
```
conda info --envs
````

# Create an enviroment for carveme, we want python 2.7.15. After running enter yes or y whenever prompted
```
conda create -n carveme_env python=2.7.15
```

# Check that new enviroment was succesfully created
```
conda info --envs
````

# Load environment
```
source activate carveme_env
````

# Now lets install carveme, available through pip. Again enter yes or y whenever prompted
```
pip install carveme
```

# Check that the packages are installed, first command shows conda packages, second command shows pip packages.
```
conda list
pip list
```

# Pip should have installed all the dependencies for carveme (framed, pandas, diamond), if not present in either of the above lists, install manually. Also see https://carveme.readthedocs.io/en/latest/installation.html , which says to run the following command after carveme installation
```
carveme_init
``` 

# Finally download and install CPLEX *inside* the carveme_env environment. To do this activate the environment.
```
source activate carveme_env
```

# Follow the instalaltion instructions of your CPLEX binary


