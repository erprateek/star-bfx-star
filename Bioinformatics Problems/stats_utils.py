import pandas as pd
from scipy.stats import beta

def generate_analytical_sensitivity_stats(success, total, confint=0.95):
    """[This method computes the sensitivity stats for a given P (number of successes) out of a given N (number of observations) ]

    Arguments:
        success {[int]} -- [Number of times success was obtained from the total number of cases]
        total {[int]} -- [N - The total number of observations in the study]

    Keyword Arguments:
        confint {float} -- [The % for which confidence interval is to be computed] (default: {0.95})

    Returns:
        [tuple] -- [Tuple A,B where A is the column headers and B is the values for those columns]
    """
    cols = ["N","TP","FN","Sensitivity","Min_CI","Max_CI"]
    n_missed = total-success
    if total == 0:
        return (cols, [total, success, n_missed, "-","-","-"])
    sensitivity = round(float(success)/float(total), 2)
    sensitivity = sensitivity*100
    lower_ci, upper_ci = binom_interval(success, total, confint)
    row_nums = [total, success, n_missed, sensitivity, lower_ci, upper_ci]
    cols = ["N","TP","FN","Sensitivity","Min_CI","Max_CI"]
    return (cols, row_nums)

def binom_interval(success, total, confint=0.95):
    """[This method computes the confidence interval for a given P (number of successes) out of a given N (number of observations)]

    Arguments:
        success {[int]} -- [Number of times success was obtained from the total number of cases]
        total {[int]} -- [N - The total number of observations in the study]

    Keyword Arguments:
        confint {float} -- [The % for which confidence interval is to be computed] (default: {0.95})

    Returns:
        [tuple] -- [Tuple containing lower CI and higher CI (in that order).]
    """
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total - success + 1)
    upper = beta.ppf(1 - quantile, success + 1, total - success)
    if success == total:
        upper = 1.00000
    if success == 0:
        lower = 0.00000
    return (round(lower,2)*100, round(upper,2)*100)
