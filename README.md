# CA1 Pyramidal with AMPA synapses

## Getting started

Clone this and set up required packages

``` sh
python -mvenv .venv
source .venv/bin/activate
pip install ipython jupyter pandas matplotlib seaborn numpy
# only needed when building arbor from source
pip install pybind11-stubgen scikit-build-core
# or build from source
pip install arbor==0.11 
```
From now on, you can bring the packages into scope by `source .venv/bin/activate`
and remove them from view by `deactivate`.

Then run the model
``` sh
python main.py
```
which will produce a plot of the membrane potential at the soma in `result.pdf`.
