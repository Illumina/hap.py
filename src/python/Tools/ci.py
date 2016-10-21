# coding=utf-8
#
# Copyright (c) 2010-2016 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt


from __future__ import division
from math import log1p

import numpy as np
import scipy.stats as stats


_VALUE_CACHE = {}




def jeffreysCI(x, n, alpha=0.05):
    '''Modified Jeffreys confidence interval for binomial proportions:
    Brown, Cai and DasGupta: Interval Estimation for a Binomial Proportion.
    2001, doi:10.1214/ss/1009213286'''

    # HAP-240 avoid division by zero
    if n == 0:
        return 0.0, 0.0, 1.0

    key = "%i_%i_%f" % (x, n, alpha)
    try:
        return _VALUE_CACHE[key]
    except:
        pass

    p = x / n
    beta = stats.distributions.beta(x+0.5, n-x+0.5)

    # lower bound
    if x == n:
        lower = (alpha/2)**(1/n)
    elif x <= 1:
        lower = 0.0
    else:
        lower = beta.ppf(alpha/2)

    # upper bound
    if x == 0:
        upper = 1-(alpha/2)**(1/n)
    elif x >= n-1:
        upper = 1.0
    else:
        upper = beta.isf(alpha/2)

    # avoid values outside the unit range due to potential numerical inaccuracy
    lower = max(lower, 0.0)
    upper = min(upper, 1.0)

    _VALUE_CACHE[key] = (p, lower, upper)

    return (p, lower, upper)


binomialCI = np.vectorize(jeffreysCI)
