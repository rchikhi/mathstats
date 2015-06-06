'''
    Estimates parameters of a truncated normal distribution given observations. 
    This truncated normal distribution has a distribution on the form given in 
    http://www.ncbi.nlm.nih.gov/pubmed/22923455 . The skewness is caused by a weight function that acts on
    the normal distribution in the defined interval.
    
    Parameters of a truncated normal distribution are:
    * mean
    * variance
    * minimum value (lower bound 'a') 
    * maximum value (upper bound 'b')
    
    Application:
    Calculates gap size between two contigs assuming normal distribution of library. The 'obs' variable 
    is here mean_lib - (o1 + o2) i.e. the naive gap estimation. (change this notation at some point) Eventually change 
    function names as well to be more general.
    '''

import sys
import math
from math import exp, pi, sqrt
from mathstats.normaldist import normal


def GapEstimator(mean, sigma, read_length, mean_obs, c1_len, c2_len=None):
    '''
    Calculates the lower bound (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455). The upper bound is then 
    uniquely determined by the lower bound. 
    '''
    naive_gap = mean - mean_obs
    if c2_len is None:
        c2_len = c1_len

    d_ML = CalcMLvaluesOfdGeneral(mean, sigma, read_length, c1_len, naive_gap, c2_len)

    return d_ML

def GapEstimator_reinforced(mean, sigma, read_length, mean_obs, stddev_obs, c1_len, c2_len=None):
    '''
    Calculates the lower bound (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455). The upper bound is then 
    uniquely determined by the lower bound. 
    '''
    if c2_len is None:
        c2_len = c1_len
        
    naive_gap = mean - mean_obs
    d_ML = CalcMLvaluesOfdGeneral(mean, sigma, read_length, c1_len, naive_gap, c2_len)
    stddev_est = stddev_given_d(mean, sigma, read_length, c1_len, c2_len, d_ML)

    return d_ML,stddev_est

def PreCalcMLvaluesOfdLongContigs(mean, stdDev, readLen):
    def Nom(z, mean, stdDev):
        nom = -(1 + normal.erf((mean - d - 2 * readLen + 1) / (2 ** 0.5 * float(stdDev)))) * (pi / 2) ** 0.5
        return nom
    def Denom(d, readLen, mean, stdDev):
        first = -((pi / 2) ** 0.5) * (d + 2 * readLen - mean - 1) * (1 + normal.erf((mean - d - 2 * readLen + 1) / (2 ** 0.5 * float(stdDev))))
        second = stdDev * exp((-((mean - d - 2 * readLen + 1) ** 2) / (float(2 * stdDev ** 2))))
        denom = first + second
        return denom

    def CalcAofd(d):
        #transform to z ~N(0,1)
        z = ((d + 2 * readLen) - mean) / float(stdDev)
        nom = Nom(z, mean, stdDev)
        denom = Denom(d, readLen, mean, stdDev)
        val = nom / denom
        return val
    #we cannot span gaps outside this interval, for fair coverages with -e set to not too low value
    d_upper = int(mean + 2 * stdDev - 2 * readLen)
    d_lower = int(-4 * stdDev)
    dValuesTable = {}
    prev_obs = d_lower
    for d in xrange(d_lower, d_upper + 1):
        #Aofd=CalcAofd(d)
        func_of_d, Aofd = funcDGeneral(d, mean, stdDev, mean + 4 * stdDev, mean + 4 * stdDev, readLen)
        obs = func_of_d #d+Aofd*stdDev**2
        obs = int(round(obs, 0))
        #print 'd_ML:', d,'Observation:',obs
        dValuesTable[obs] = d
        #fill in missing values of the table here
        if abs(obs - prev_obs) > 1:
            n = abs(obs - prev_obs)
            for i in xrange(0, n):
                dValuesTable[prev_obs + i + 1] = d

        prev_obs = obs
    return dValuesTable

def funcDGeneral(d, mean, stdDev, c1Len, c2Len, readLen):
    c_min = min(c1Len, c2Len)
    c_max = max(c1Len, c2Len)

    def Nominator(d, c_min, c_max, c1Len, c2Len, readLen):

        term1 = -0.5 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) + normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
        term2 = +0.5 * (normal.erf((d + c_max + readLen - mean) / (2 ** 0.5 * float(stdDev))) + normal.erf((d + c_min + readLen - mean) / (2 ** 0.5 * float(stdDev))))
        g_prime_d = term1 + term2
        return -g_prime_d

    def Denominator(d, c1Len, c2Len, readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...

        term1 = (c_min - readLen + 1) / 2.0 * (normal.erf((c_max + d + readLen - mean) / ((2 ** 0.5) * stdDev)) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

        term2 = (c_min + c_max + d - mean + 1) / 2.0 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_max + readLen + d - mean) / (2 ** 0.5 * float(stdDev))))

        term3 = (d + 2 * readLen - mean - 1) / 2.0 * (normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

        term4 = stdDev / ((2 * pi) ** 0.5) * (exp(-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) + exp(-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2))))

        term5 = -stdDev / ((2 * pi) ** 0.5) * (exp(-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))) + exp(-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))))
        g_d = term1 + term2 + term3 + term4 + term5
        return g_d

    denominator = Denominator(d, c1Len, c2Len, readLen)
    nominator = Nominator(d, c_min, c_max, c1Len, c2Len, readLen)
    try:
        Aofd = nominator / denominator
    except ZeroDivisionError:
        Aofd = 0.0000000001

    func_of_d = d + Aofd * stdDev ** 2
    return func_of_d, Aofd

def ML_one_reference(d, mean, stdDev, reflen, readLen):
    def Nominator(d, reflen, readLen):
        term1 = exp(-((reflen + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) * (-mean * (d + 2 * readLen - 1) + (reflen + d + 1) * (d + 2 * readLen + mean - 1) + stdDev ** 2 - (reflen + d + 1) ** 2)
        term2 = exp(-((readLen + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) * (-mean * (d + 2 * readLen - 1) + (readLen + d + 1) * (d + 2 * readLen + mean - 1) + stdDev ** 2 - (readLen + d + 1) ** 2)
        #  (normal.erf((reflen + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) #- normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
        #term1 = 0.5 * (normal.erf((reflen + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
        g_prime_d = term1 - term2
        #print term1, term2
        return g_prime_d

    def Denominator(d, reflen, readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...

        #term1 = (c_min - readLen + 1) / 2.0 * (normal.erf((c_max + d + readLen - mean) / ((2 ** 0.5) * stdDev)) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

        #term2 = (c_min + c_max + d - mean + 1) / 2.0 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_max + readLen + d - mean) / (2 ** 0.5 * float(stdDev))))

        term3 = ((d + 2 * readLen - mean - 1) / 2.0) * (normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((reflen + d + 1 - mean) / ((2 ** 0.5) * stdDev)))

        term4 = stdDev / ((2 * pi) ** 0.5) * (-exp(-((reflen + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) + exp(-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2))))

        #term5 = -stdDev / ((2 * pi) ** 0.5) * (exp(-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))) + exp(-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))))
        g_d = term3 + term4
        return g_d

    denominator = Denominator(d, reflen, readLen)
    nominator = Nominator(d, reflen, readLen)
    Aofd = nominator / denominator
    func_of_d = d - Aofd * stdDev ** 2
    #print 'nom:', nominator, 'denom:', denominator
    return func_of_d


def CalcMLvaluesOfdGeneral(mean, stdDev, readLen, c1Len, naive_gap, c2Len):
    """
        returns the ML-gap as float. The ML gap is searched
        for with an accuracy of 0.1bp. This is NOT the accuracy of the ML-estimate 
        however since it depends on the number of samples etc. 
    """
    #do binary search among values
    d_upper = int(mean + 2 * stdDev - 2 * readLen)
    d_lower = int(-c1Len - c2Len) + int(max(mean - 2 * stdDev, 2 * readLen))
    #print naive_gap
    while d_upper - d_lower > 0.1:
        d_ML = (d_upper + d_lower) / 2.0
        if c2Len:
            func_of_d, Aofd = funcDGeneral(d_ML, mean, stdDev, c1Len, c2Len, readLen)
        else:
            func_of_d = ML_one_reference(d_ML, mean, stdDev, c1Len, readLen)
            #print 'current gap:', d_ML
            #print func_of_d
        if func_of_d > naive_gap:
            d_upper = d_ML
        else:
            d_lower = d_ML

    d_ML = (d_upper + d_lower) / 2.0
    return d_ML




def tr_sk_mean(d, mean, stdDev, c1Len, c2Len, readLen):
    '''
    Calculates the mean of the truncated normal distribution (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455) with lower bound d, 
     and upper bound (d+c1Len+c2Len)
    '''
    def CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen):

        term1 = -0.5 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) + normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
        term2 = +0.5 * (normal.erf((d + c_max + readLen - mean) / (2 ** 0.5 * float(stdDev))) + normal.erf((d + c_min + readLen - mean) / (2 ** 0.5 * float(stdDev))))
        g_prime_d = term1 + term2
        return g_prime_d

    def CalcGd(d, c1Len, c2Len, readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term1 = (c_min - readLen + 1) / 2.0 * (normal.erf((c_max + d + readLen - mean) / ((2 ** 0.5) * stdDev)) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

        term2 = (c_min + c_max + d - mean + 1) / 2.0 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_max + readLen + d - mean) / (2 ** 0.5 * float(stdDev))))

        term3 = (d + 2 * readLen - mean - 1) / 2.0 * (normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

        term4 = stdDev / ((2 * pi) ** 0.5) * ( exp( (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2)))) + exp( (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2)))))

        term5 = -stdDev / ((2 * pi) ** 0.5) * (exp( (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2)))) + exp( (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2)))))

        g_d = term1 + term2 + term3 + term4 + term5
        return g_d

    c_min = min(c1Len, c2Len)
    c_max = max(c1Len, c2Len)
    g_prime_d = CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen)
    g_d = CalcGd(d, c1Len, c2Len, readLen)
    a = stdDev ** 2 * g_prime_d + mean * g_d
    E_o = a / g_d - d
    return E_o



def tr_sk_std_dev(mean, stdDev, readLen, c1Len, c2Len, d):
    '''
    Calculates the standard deviation of the truncated normal distribution (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455)
     with lower bound d, upper bound (d+c1Len+c2Len).
    '''
    return stddev_given_d(mean, stdDev, readLen, c1Len, c2Len, d)

    # def E_O_square(d, mean, stdDev, c1Len, c2Len, readLen):

    #     def CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen):

    #         term1 = -0.5 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) + normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
    #         term2 = +0.5 * (normal.erf((d + c_max + readLen - mean) / (2 ** 0.5 * float(stdDev))) + normal.erf((d + c_min + readLen - mean) / (2 ** 0.5 * float(stdDev))))
    #         g_prime_d = term1 + term2
    #         return g_prime_d

    #     def CalcGd(d, c1Len, c2Len, readLen):
    #         #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
    #         term1 = (c_min - readLen + 1) / 2.0 * (normal.erf((c_max + d + readLen - mean) / ((2 ** 0.5) * stdDev)) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

    #         term2 = (c_min + c_max + d - mean + 1) / 2.0 * (normal.erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_max + readLen + d - mean) / (2 ** 0.5 * float(stdDev))))

    #         term3 = (d + 2 * readLen - mean - 1) / 2.0 * (normal.erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))) - normal.erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

    #         term4 = stdDev / ((2 * pi) ** 0.5) * ( exp( (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))))   + exp( (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2)))))

    #         term5 = -stdDev / ((2 * pi) ** 0.5) * ( exp( (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2)))) + exp( (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2)))))
    #         g_d = term1 + term2 + term3 + term4 + term5
    #         return g_d

    #     def CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d):
    #         c1 = 0
    #         c2 = normal.normpdf(d + 2 * readLen - 1, mean, stdDev) - normal.normpdf(c_min + d + readLen , mean, stdDev) - normal.normpdf(c_max + d + readLen , mean, stdDev) + normal.normpdf(c_min + c_max + d + 1 , mean, stdDev)

    #         b = stdDev ** 4 * (c1 + c2) + (mean ** 2 + stdDev ** 2) * g_d + 2 * stdDev ** 2 * mean * g_prime_d
    #         return b, c1, c2

    #     c_min = min(c1Len, c2Len)
    #     c_max = max(c1Len, c2Len)
    #     g_prime_d = CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen)
    #     g_d = CalcGd(d, c1Len, c2Len, readLen)
    #     #a = stdDev ** 2 * g_prime_d + mean * g_d
    #     b, c1, c2 = CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d)
    #     E_o_square = b / g_d - 2 * d * (stdDev ** 2 * (g_prime_d / g_d) + mean) + d ** 2
    #     return E_o_square



    # e_o = tr_sk_mean(d, mean, stdDev, c1Len, c2Len, readLen)
    # e_o_square = E_O_square(d, mean, stdDev, c1Len, c2Len, readLen)

    # if e_o_square - e_o ** 2 < 0:
    #     sys.stderr.write('Error in estimate, Std_dev(O|d =' + str(d) + ')=' + str(e_o_square) + '  ' + str(e_o ** 2) + str(c1Len) + ' ' + str(c2Len))
    #     std_dev = 0
    # else:
    #     std_dev = (e_o_square - e_o ** 2) ** 0.5
    # return(std_dev)




def erf_norm(x,mu,sigma):
    return normal.erf((x-mu)/(sqrt(2.0)*sigma))    

def e_norm_neg(x,mu,sigma):
    return exp(-(x-mu)**2/(2*sigma**2))

def e_norm_pos(x,mu,sigma):
    return exp((x-mu)**2/(2*sigma**2))

def e_sq_norm_pos(x,mu,sigma):
    try:
        return exp((mu**2+x**2)/(2*sigma**2))
    except OverflowError:
        import decimal
        decimal.getcontext().prec = 60
        return float(decimal.Decimal((mu**2+x**2)/(2*sigma**2)).exp())

def e_sq_norm_neg(x,mu,sigma):
    try:
        return exp(-(mu**2+x**2)/(2*sigma**2))
    except OverflowError:
        import decimal
        decimal.getcontext().prec = 60
        return float(decimal.Decimal(-(mu**2+x**2)/(2*sigma**2)).exp())

def e_times_pos(x,mu,sigma):
    try:
        return exp((mu*x)/sigma**2)
    except OverflowError:
        import decimal
        decimal.getcontext().prec = 60
        # print float(float(decimal.Decimal((mu*x)/sigma**2).exp())
        return float(decimal.Decimal((mu*x)/sigma**2).exp())   


 
def b1(x,mu,sigma,r,d,c_min,c_max):
    return (1/2.0)*sigma* e_sq_norm_neg(x,mu,sigma) * \
    ( sqrt(2*pi)*e_sq_norm_pos(x,mu,sigma)* ( (-mu**2 + -sigma**2)*(d +2*r -1) + mu**3 + 3*mu*sigma**2 )* erf_norm(x,mu,sigma) -\
     2*sigma*e_times_pos(x,mu,sigma) *( -mu*(d + 2*r - 1) + x*(-d + mu -2*r +1) + mu**2 +2*sigma**2 + x**2 ) )  

def b2(x,mu,sigma,r,d,c_min,c_max):
    return (1/2.0)*sigma* e_sq_norm_neg(x,mu,sigma)*(c_min -r +1) * \
    ( sqrt(2*pi)*  (mu**2 + sigma**2) * e_sq_norm_pos(x,mu,sigma) * erf_norm(x,mu,sigma) -\
     2*sigma*(mu + x)*e_times_pos(x,mu,sigma) ) 

def b3(x,mu,sigma,r,d,c_min,c_max):
    return  (1/2.0)*sigma* e_sq_norm_neg(x,mu,sigma) * \
    ( 2*sigma*e_times_pos(x,mu,sigma) *( -mu*(c_min + c_max + d +1) -x*( c_min + c_max + d -mu +1 ) + mu**2 + 2*sigma**2 + x**2 ) -\
    sqrt(2*pi)*e_sq_norm_pos(x,mu,sigma)*erf_norm(x,mu,sigma)* ( -mu**2*(c_min + c_max + d +1) -sigma**2*( c_min + c_max + d +1 ) +mu**3 + 3*mu*sigma**2 ))


def a1(x,mu,sigma,r,d,c_min,c_max):
    return sqrt(pi/2.0)*sigma*( -mu*(d +2*r -1) + mu**2 + sigma**2 )* erf_norm(x,mu,sigma) -\
    sigma**2 * e_norm_neg(x,mu,sigma) * (-d + mu  - 2*r + x +1)

def a2(x,mu,sigma,r,d,c_min,c_max):
    return (c_min -r +1)* ( sqrt(pi/2.0)*mu*sigma*erf_norm(x,mu,sigma) - sigma**2*e_norm_neg(x,mu,sigma) )

def a3(x,mu,sigma,r,d,c_min,c_max):
    return sigma**2* (-e_norm_neg(x,mu,sigma))*( c_min + c_max + d - mu -x + 1 ) -\
    sqrt(pi/2.0)*sigma*erf_norm(x,mu,sigma)*( -mu*( c_min + c_max + d +1 ) + mu**2 + sigma**2 )

def g_d1(x,mu,sigma,r,d,c_min,c_max):
    return sigma**2* (-e_norm_neg(x,mu,sigma)) - sqrt(pi/2.0)*sigma*(d - mu + 2*r -1)* erf_norm(x,mu,sigma)

def g_d2(x,mu,sigma,r,d,c_min,c_max):
    return sqrt(pi/2.0)*sigma*erf_norm(x,mu,sigma)*(c_min - r + 1)

def g_d3(x,mu,sigma,r,d,c_min,c_max):
    return sqrt(pi/2.0)*sigma * (c_min +c_max + d - mu +1 )* erf_norm(x,mu,sigma) +\
    sigma**2*e_norm_neg(x,mu,sigma)




def stddev_given_d(mean, stdDev, readLen, c1Len, c2Len, d):
    """
        b/g_d - (a/g_d)^2 according to
        http://www.biomedcentral.com/content/supplementary/1471-2105-15-281-s1.pdf
    """
    stdDev = float(stdDev)
    mean = float(mean)
    c_min = min(c1Len,c2Len)
    c_max = max(c1Len,c2Len)

    #assert c_min > readLen - 1

    # define integration points (x1, x2, x3, x4) 
    # the functions a, b and g_d each have three different 
    # intervals they need to be integrated in (because of the max(min(.,.,.)) function)
    # intervals: [x1,x2], [x2,x3], [x3,x4]. 
    x1 = d + 2 * readLen - 1
    x2 = c_min + d +readLen if c_min + d + readLen < mean + 5*stdDev else mean + 5*stdDev
    x3 = c_max + d +readLen if c_max + d + readLen < mean + 5*stdDev else mean + 5*stdDev
    x4 =  c_max +c_min + d + 1 if c_min + c_max + d + 1 < mean + 5*stdDev else mean + 5*stdDev

    args = (mean,stdDev,readLen,d,c_min,c_max)

    # print b1(x2, *args)
    # print b1(x1, *args)
    # print b1(x2, *args) - b1(x1, *args)

    # print b2(x3, *args)
    # print b2(x2, *args)
    # print b2(x3, *args) - b2(x2, *args)

    # print b3(x4, *args)
    # print b3(x3, *args)
    # print b3(x4, *args) - b3(x3, *args)

    b_term = b1(x2, *args) - b1(x1, *args) \
            + b2(x3, *args) - b2(x2, *args) \
            + b3(x4, *args) - b3(x3, *args)

    # print 'b:', b_term

    # print a1(x2, *args)
    # print a1(x1, *args)
    # print a1(x2, *args) - a1(x1, *args)

    # print a2(x3, *args)
    # print a2(x2, *args)
    # print a2(x3, *args) - a2(x2, *args)

    # print a3(x4, *args)
    # print a3(x3, *args)
    # print a3(x4, *args) - a3(x3, *args)

    a_term = a1(x2, *args) - a1(x1, *args) \
            + a2(x3, *args) - a2(x2, *args) \
            + a3(x4, *args) - a3(x3, *args)
    # print 'a:',a_term


    # print g_d1(x2, *args)
    # print g_d1(x1, *args)
    # print g_d1(x2, *args) - g_d1(x1, *args)

    # print g_d2(x3, *args)
    # print g_d2(x2, *args)
    # print g_d2(x3, *args) - g_d2(x2, *args)

    # print g_d3(x4, *args)
    # print g_d3(x3, *args)
    # print g_d3(x4, *args) - g_d3(x3, *args)

    g_d_term = g_d1(x2, *args) - g_d1(x1, *args) \
            + g_d2(x3, *args) - g_d2(x2, *args) \
            + g_d3(x4, *args) - g_d3(x3, *args)
    # print 'g_d:',g_d_term



    var = b_term / g_d_term\
    - ( a_term / g_d_term )**2

    #print sqrt(var)
    if var < 0:
        sys.stderr.write('Error in estimate, a={0}, b={1}\n'.format(a_term,b_term))
        std_dev_est = 0
    elif math.isnan(var):
        sys.stderr.write('MathOverflow - too small stddev of library ({0}bp) compared to the mean ({1}bp) - check distribution. \n'.format(stdDev,mean))
        std_dev_est = 0
    else:
        std_dev_est = sqrt(var)  
    return std_dev_est


def mean_given_d(mean, stdDev, readLen, c1Len, c2Len, d):
    """
        a/g_d - d according to
        http://www.biomedcentral.com/content/supplementary/1471-2105-15-281-s1.pdf
    """

    stdDev = float(stdDev)
    mean = float(mean)
    c_min = min(c1Len,c2Len)
    c_max = max(c1Len,c2Len)

    #assert c_min > readLen - 1

    # define integration points (x1, x2, x3, x4) 
    # the functions a, b and g_d each have three different 
    # intervals they need to be integrated in (because of the max(min(.,.,.)) function)
    # intervals: [x1,x2], [x2,x3], [x3,x4]. 
    x1 = d + 2 * readLen - 1
    x2 = c_min + d +readLen if c_min + d + readLen < mean + 5*stdDev else mean + 5*stdDev
    x3 = c_max + d +readLen if c_max + d + readLen < mean + 5*stdDev else mean + 5*stdDev
    x4 =  c_max +c_min + d + 1 if c_min + c_max + d + 1 < mean + 5*stdDev else mean + 5*stdDev

    args = (mean,stdDev,readLen,d,c_min,c_max)

    a_term = a1(x2, *args) - a1(x1, *args) \
            + a2(x3, *args) - a2(x2, *args) \
            + a3(x4, *args) - a3(x3, *args)

    g_d_term = g_d1(x2, *args) - g_d1(x1, *args) \
            + g_d2(x3, *args) - g_d2(x2, *args) \
            + g_d3(x4, *args) - g_d3(x3, *args)
        
    return a_term / g_d_term -d   
