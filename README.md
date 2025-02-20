# cdFBA
Community Dynamic Flux Balance Analysis using process-bigraph.

Release version 0.0.1 is available from PyPI. Requires python>=3.9

To install use:

`pip install cdFBA`

### FBA Implementation

`cdFBA` uses the `cobrapy` repository to perform FBA analysis. 

### Dynamics

Static FBA optimizations are used to dynamically change the substrate concentrations in the environment at each time-step, based on substrate exchange reaction
kinetics. The subsatrate concentrations are then used to set flux boundaries for the FBA optimization in the following time-step.

### Utility Functions

cdFBA includes some basic utility functions to create the model specification dictionaries automatically. 
Just provide a list of model files and exchange reaction IDs.
The functions can be found [here](https://github.com/vivarium-collective/cdFBA/blob/main/cdFBA/utils.py).

**!!CAUTION!!**: The utility functions work only for metabolic models with the same annotation systems - i.e. same reaction name and IDs.
If you need to use models from different sources with different annotation systems, you will need to manually map the substrate/metabolite names.

### Demo Notebooks

Demo Jupyter notebooks can be found [here](https://github.com/vivarium-collective/cdFBA/tree/main/Notebooks). 
The primary demo using the process-bigraph interface can be found [here](https://github.com/vivarium-collective/cdFBA/blob/main/Notebooks/demo.ipynb).
We will be adding more notebooks to demonstrate usage reproduction of results from the literature. 



### !!This project is a work in progress!!
