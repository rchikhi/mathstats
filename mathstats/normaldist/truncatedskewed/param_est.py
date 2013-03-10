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
from math import exp, pi
from mathstats.normaldist import normal


def GapEstimator(mean, sigma, read_length, mean_obs, c1_len, c2_len):
    '''
    Calculates the lower bound (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455). The upper bound is then 
    uniquely determined by the lower bound. 
    '''
    naive_gap = mean - mean_obs
    d_ML = CalcMLvaluesOfdGeneral(mean, sigma, read_length, c1_len, c2_len, naive_gap)
    return d_ML

def PreCalcMLvaluesOfdLongContigs(mean, stdDev, readLen):
    def Nom(z, mean, stdDev):
        nom = -(1 + normal.erf((mean - d - 2 * readLen + 1) / (2 ** 0.5 * float(stdDev)))) * (pi / 2) ** 0.5
        return nom
    def Denom(d, readLen, mean, stdDev):
        first = -((pi / 2) ** 0.5) * (d + 2 * readLen - mean - 1) * (1 + normal.erf((mean - d - 2 * readLen + 1) / (2 ** 0.5 * float(stdDev))))
        second = stdDev * 2.718 ** (-((mean - d - 2 * readLen + 1) ** 2) / (float(2 * stdDev ** 2)))
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
    print Aofd
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
    Aofd = nominator / denominator
    func_of_d = d + Aofd * stdDev ** 2
    return func_of_d, Aofd

def CalcMLvaluesOfdGeneral(mean, stdDev, readLen, c1Len, c2Len, obs):
    #do binary search among values
    d_upper = int(mean + 2 * stdDev - 2 * readLen)
    d_lower = int(-4 * stdDev)
    while d_upper - d_lower > 1:
        d_ML = (d_upper + d_lower) / 2.0
        func_of_d, Aofd = funcDGeneral(d_ML, mean, stdDev, c1Len, c2Len, readLen)
        if func_of_d > obs:
            d_upper = d_ML
        else:
            d_lower = d_ML

    d_ML = (d_upper + d_lower) / 2.0
    return int(round(d_ML, 0))




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

        term4 = stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2))))

        term5 = -stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))))

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
    def E_O_square(d, mean, stdDev, c1Len, c2Len, readLen):

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

            term4 = stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2))))

            term5 = -stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))))
            g_d = term1 + term2 + term3 + term4 + term5
            return g_d

        def CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d):
            c1 = 0
            c2 = normal.normpdf(d + 2 * readLen - 1, mean, stdDev) - normal.normpdf(c_min + d + readLen , mean, stdDev) - normal.normpdf(c_max + d + readLen , mean, stdDev) + normal.normpdf(c_min + c_max + d + 1 , mean, stdDev)

            b = stdDev ** 4 * (c1 + c2) + mean ** 2 * g_d + stdDev ** 2 * g_d + 2 * stdDev ** 2 * mean * g_prime_d
            return b, c1, c2

        c_min = min(c1Len, c2Len)
        c_max = max(c1Len, c2Len)
        g_prime_d = CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen)
        g_d = CalcGd(d, c1Len, c2Len, readLen)
        a = stdDev ** 2 * g_prime_d + mean * g_d
        b, c1, c2 = CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d)
        E_o_square = b / g_d - 2 * d * (stdDev ** 2 * (g_prime_d / g_d) + mean) + d ** 2
        return E_o_square



    e_o = tr_sk_mean(d, mean, stdDev, c1Len, c2Len, readLen)
    e_o_square = E_O_square(d, mean, stdDev, c1Len, c2Len, readLen)
    if e_o_square - e_o ** 2 < 0:
        sys.stderr.write('Error in estimate, Std_dev(O|d =' + str(d) + ')=' + str(e_o_square) + '  ' + str(e_o ** 2) + str(c1Len) + ' ' + str(c2Len))
        std_dev = 0
    else:
        std_dev = (e_o_square - e_o ** 2) ** 0.5
    return(std_dev)




