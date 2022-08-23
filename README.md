# Bootstrap Confidence Intervals

Bootstrap confidence intervals for detecting event-related signals in time series data using python.

## Overview

* Motivated by: <https://www.frontiersin.org/articles/10.3389/fnmol.2020.00014/full>
* For the MATLAB based version check out the [repo from the authors](https://github.com/philjrdb/ERTsimulation).
* Note, the `data/ERT_testdata.mat` was downloaded from: <https://github.com/philjrdb/ERTsimulation/blob/master/ERT_testdata.mat>

## Installation

After cloning the repository you can install via `pip` by going to the repo folder in your terminal and entering:

```{bash}
pip install -e .
```

## Usage

1. Initialize the BootCI class. Input data should be of the form $m \times t$, where each of the $m$ rows is a single subject's trial-averaged data measured at $t$ different time points.
2. Next compute CI using bootstrap resampling with `n_boots` iterations (default `n_boots=1000`).
3. Finally, get the statistically significant events with `get_sig_events(num_consec)`. The `num_consec` parameter is used to ensure that the CI for the event-related transient is significant for some minimum duration.

```{python}
from bootstrap_ci import BootCI
my_boots = BootCI(my_data)
my_boots.bootstrap_CI(n_boots)
my_boots.get_sig_events(num_consec)
```

## In progress

* [ ] Implement resampling method for taking the difference between two samples.
