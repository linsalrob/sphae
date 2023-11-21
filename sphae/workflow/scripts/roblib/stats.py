import numpy
import sys
from math import sqrt

__author__ = 'Rob Edwards'



def mean(arr):
    """
    Calculate the mean of all the values in a list

    :param arr: The list of values
    :type arr: list
    :return: The mean of all the numbers
    :rtype: float
    """
    n = len(arr)
    sm = sum(arr)
    return float(sm) / float(n)


def median(arr):
    """
    Calculate the median of all the values in a list

    :param arr: The list of values
    :type arr: list
    :return: The median of all the numbers
    :rtype: float
    """
    arr.sort()
    n = len(arr)
    mid = int(n / 2)
    return arr[mid]


def stdev(arr):
    """
    Calculate the standard deviation of all the values in a list

    :param arr: The list of values
    :type arr: list
    :return: The st dev of all the numbers
    :rtype: float
    """
    n = len(arr)
    sm = sum(arr)
    nmpyarr = numpy.array(arr)
    sum_of_squares = numpy.sum(nmpyarr * nmpyarr)
    result = 0
    try:
        result = sqrt(sum_of_squares / n - (sm / n) ** 2)
    except:
        sys.stderr.write("can't create the stdev. Only have " + str(n) + " elements in the array")

    return result

def confidence_interval_calcs(mean, stdev, n, stdevs=2):
    """
    Return a tuple of (lower, upper) confidence interavals based on the mean and the standard deviation.

    Uses [this calculation](https://stackoverflow.com/a/15034101/1165471) to get the confidence intervals:
confidence interval is `mean +/- z*sigma`, where `sigma` is the estimated standard deviation of your sample mean, given by `sigma = s / sqrt(n)`, where `s` is the standard deviation computed from your sample data and `n` is your sample size
    :param mean: the mean of the sample
    :param stdev: the standard deviation of the sample
    :param n: the number of observations in the sample
    :param stdevs: (0 < stdevs < 5) the number of standard deviations from the mean
    :return: tuple of (lower bound, upper bound)
    """

    stoz = {1: 0.34134, 2: 0.47725, 3: 0.49865, 4: 0.49997}

    if stdevs not in stoz:
        sys.stderr.write("ERROR: The number of standard deviations from the mean must be an integer between 1 and 4\n")
        return (None, None)

    sigma = stdev / sqrt(n)
    lower = mean - (stoz[stdevs] * sigma)
    upper = mean + (stoz[stdevs] * sigma)

    return (lower, upper)

def confidence_intervals(arr, stdevs=2):
    """
    Calculate the confidence interval for an array. The stdevs is the number of stdevs from the mean for the interval
    :param arr: the list of values
    :param stdevs: how many stdevs from the mean. Must be between 1 and 4
    :return: tuple of (lower bound, upper bound)
    """

    mn = mean(arr)
    st = stdev(arr)
    return confidence_interval_calcs(mn, st, len(arr), stdevs)


