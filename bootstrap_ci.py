"""
Implementation of bootstrap CIs.

[See this paper](https://www.frontiersin.org/articles/10.3389/fnmol.2020.00014/full) for more info.
"""
from itertools import groupby

import numpy as np


class BootCI:
    """Bootstrap confidence intervals (CI) for event-related transients."""

    def __init__(self, data, n_boots=None):
        """
        pd.DataFrame must be orinted n x p.

        Each row is a subject's trial-averaged data consisting of p distinct time points.
        """
        self.data = data
        self.n_samples = data.shape[0]
        self.boot_samples = None
        self.ci_lower = None
        self.ci_upper = None
        self.sig_events = {}

        if n_boots:
            self._bootstrap_resample(n_boots)

    def _bootstrap_resample(self, n_boots):
        """
        Perform bootstrap resampling of input df.

        Note: Data must be (n, t), where each row is one subject's data at each of t time points.

        Args:
            df (DataFrame): DataFrame to resample.
            n_boots (int): Number of bootstrap samples.

        Returns:
            np.array: array containing bootstrap samples.
        """
        n_samples = self.data.shape[0]
        bootstraps = []
        for _ in range(n_boots):
            bootstraps.append(self.data.sample(n_samples, replace=True).mean().values)

        self.boot_samples = np.asarray(bootstraps)

    def bootstrap_CI(self, n_boots, sig=0.05):
        """
        Compute bootstrap confidence intervals (bCI) from inpute data.

        Args:
            df (DataFrame): data to compute bCIs on.
            n_boots (int): Number of bootstrap samples. Defaults to 1000.
            sig (float): Significance level used to detect events. Defaults to 0.05.

        Returns:
            dict: dictionary containing the bootstrap samples and upper/lower CI
        """
        self._bootstrap_resample(n_boots)
        boot_arr = self.boot_samples
        boot_arr.sort(axis=0)
        # get CI cutoffs
        lower_CI_idx = int(np.ceil(n_boots * (sig / 2)))
        upper_CI_idx = int(np.floor(n_boots * (1 - (sig / 2))))

        ci_lower = boot_arr[lower_CI_idx]
        ci_upper = boot_arr[upper_CI_idx - 1]
        # adjust CI for narrowness bias
        # See: https://www.tandfonline.com/doi/full/10.1080/00031305.2015.1089789
        # bootstrap CIs have a narrowness bias = sqrt(n/(n-1)), so multiply by 1/bias
        n_samples = self.n_samples
        ci_factor = np.sqrt(n_samples / (n_samples - 1))
        # get difference between the expanded CI and original CI
        ci_change = ci_factor * (ci_upper - ci_lower) - (ci_upper - ci_lower)
        self.ci_lower = ci_lower - ci_change / 2
        self.ci_upper = ci_upper + ci_change / 2

    def get_sig_events(self, num_consec, threshold=0):
        """
        Use a consecutive threshold to identify statistically significant events.

        Events must be significant (i.e, lower CI > 0, upper CI < 0) for at least `num_consec`
        consecutive data points.

        Args:
            data (dict): Dictionary with keys for "lower" and "upper" CIs from bootstrap_CI.
            num_consec (int): Number of consecutive points required above threshold.
            threshold (int, float): Array entries must be greater than threshold. Defaults to 0.

        Returns:
            dict: Dictionary with statistically significant events.
        """
        for ci in ["lower", "upper"]:
            arr = self.ci_lower
            if ci == "upper":
                arr = -self.ci_upper  # invert values to use same thresholding as lower CI
            groups = []
            uniquekeys = []
            # separate for above and below CIs
            for k, g in groupby(range(len(arr)), lambda x: arr[x] > threshold):
                groups.append(list(g))  # Store group iterator as a list
                uniquekeys.append(k)  # each group has associated boolean value

            self.sig_events[ci] = [
                g for g, g_bool in zip(groups, uniquekeys) if (g_bool and len(g) > num_consec)
            ]
