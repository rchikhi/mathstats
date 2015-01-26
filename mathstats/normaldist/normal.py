'''
    Created on Jun 2, 2012
    
    @author: ksahlin
    
    This file is part of BESST.
    
    BESST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    BESST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with BESST.  If not, see <http://www.gnu.org/licenses/>.
    '''

''' Here are some functions related to the normal distribution
    '''
import math
#from math import exp, sqrt, pi, log
from decimal import Decimal, getcontext
from scipy.stats import norm
# from scipy.special import erf as errorfunction

def erf(x):
    return math.erf(x)
    # ## error function approximation with an error less than 1.5 * 10-7 for all x
    # # save the sign of x
    # sign = 1 if x >= 0 else -1
    # x = abs(x)

    # # constants
    # a1 = 0.254829592
    # a2 = -0.284496736
    # a3 = 1.421413741
    # a4 = -1.453152027
    # a5 = 1.061405429
    # p = 0.3275911

    # # A&S formula 7.1.26
    # t = 1.0 / (1.0 + p * x)
    # y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x * x)
    # return sign * y # erf(-x) = -erf(x)

def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    t = 1. / (1. + 0.5 * z)
    r = t * math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (.37409196 +
        t * (.09678418 + t * (-.18628806 + t * (.27886807 +
        t * (-1.13520398 + t * (1.48851587 + t * (-.82215223 +
        t * .17087277)))))))))
    if (x >= 0.):
        return r
    else:
        return 2. - r

def normcdf(x, mu, sigma):
    t = x - mu;
    y = 0.5 * erfcc(-t / (sigma * math.sqrt(2.0)));
    if y > 1.0:
        y = 1.0;
    return y

def normpdf(x, mu, sigma, prec=False, approx=True):
    if prec:
    #Get much better approximations with Decimal (simply more decimals)
        getcontext().prec = 100
        u = Decimal(str(x - mu)) / Decimal(str(abs(sigma)))
        y = float(str((1 / Decimal(str((math.sqrt(2 * math.pi) * abs(sigma))))) * Decimal(str(-u * u / 2)).exp()))
        return y
    elif approx:
        u = (x - mu) / abs(sigma)
        y = (1 / (math.sqrt(2 * math.pi) * abs(sigma))) * math.exp(-u * u / 2)
        return y
    else:
        return norm.pdf(x, mu, sigma)

def MaxObsDistr(nr_of_obs, prob):
    # Here, we choose the quantile that separates "normal" contigs from repeats. We want to find 
    # the k s.t. the probability of marking one or more normal contigs (from the hole set of contigs) 
    # as repeats (false positives) is less than p=0.05. This is equivalent to: Choosing p in Bin(n,p)
    # (where n = nr of contigs) s.t. P(k=0)>= 0.95 (no successes if a success is a false positive).
    # We get P(k=0)= choose(n,k)*p**n => p = 1 - (0.95*n)**n. With this value of p, if X~N(mean,sigma),
    # we want to derive k from P_x(x < mean + k*sigma) = 1-p. This gives the k that is returned from this function. 
    #from scipy.stats import norm
    def rational_approximation(t):

        # Abramowitz and Stegun formula 26.2.23.
        # The absolute value of the error should be less than 4.5 e-4.
        c = [2.515517, 0.802853, 0.010328]
        d = [1.432788, 0.189269, 0.001308]
        numerator = (c[2] * t + c[1]) * t + c[0]
        denominator = ((d[2] * t + d[1]) * t + d[0]) * t + 1.0
        return t - numerator / denominator

    def normal_CDF_inverse(p):

        assert p > 0.0 and p < 1

        # See article above for explanation of this section.
        if p < 0.5:
            # F^-1(p) = - G^-1(p)
            return -rational_approximation(math.sqrt(-2.0 * math.log(p)))
        else:
            # F^-1(p) = G^-1(1-p)
            return rational_approximation(math.sqrt(-2.0 * math.log(1.0 - p)))

    p = 1 - (prob) ** (1 / float(nr_of_obs))
    k = normal_CDF_inverse(1 - p)

    return(k)

