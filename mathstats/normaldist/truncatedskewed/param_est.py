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
from math import exp, pi, sqrt,log
from mathstats.normaldist import normal


def GapEstimator(mean, sigma, read_length, mean_obs, c1_len, c2_len=None):
    '''
    Calculates the lower bound (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455). The upper bound is then 
    uniquely determined by the lower bound. 
    '''
    naive_gap = mean - mean_obs

    d_ML = CalcMLvaluesOfdGeneral(mean, sigma, read_length, c1_len, naive_gap, c2_len)

    return d_ML

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
    Aofd = nominator / denominator
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


def CalcMLvaluesOfdGeneral(mean, stdDev, readLen, c1Len, obs, c2Len):
    """
        returns the ML-gap as float. The ML gap is searched
        for with an accuracy of 0.1bp. This is NOT the accuracy of the ML-estimate 
        however since it depends on the number of samples etc. 
    """
    #do binary search among values
    d_upper = int(mean + 4 * stdDev - 2 * readLen)
    d_lower = int(-c1Len - c2Len) + int(max(mean - 4 * stdDev, 2 * readLen))
    #print obs
    while d_upper - d_lower > 0.1:
        d_ML = (d_upper + d_lower) / 2.0
        if c2Len:
            func_of_d, Aofd = funcDGeneral(d_ML, mean, stdDev, c1Len, c2Len, readLen)
        else:
            func_of_d = ML_one_reference(d_ML, mean, stdDev, c1Len, readLen)
            #print 'current gap:', d_ML
            #print func_of_d
        if func_of_d > obs:
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

            term4 = stdDev / ((2 * pi) ** 0.5) * ( exp( (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))))   + exp( (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2)))))

            term5 = -stdDev / ((2 * pi) ** 0.5) * ( exp( (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2)))) + exp( (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2)))))
            g_d = term1 + term2 + term3 + term4 + term5
            return g_d

        def CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d):
            c1 = 0
            c2 = normal.normpdf(d + 2 * readLen - 1, mean, stdDev) - normal.normpdf(c_min + d + readLen , mean, stdDev) - normal.normpdf(c_max + d + readLen , mean, stdDev) + normal.normpdf(c_min + c_max + d + 1 , mean, stdDev)

            b = stdDev ** 4 * (c1 + c2) + (mean ** 2 + stdDev ** 2) * g_d + 2 * stdDev ** 2 * mean * g_prime_d
            return b, c1, c2

        c_min = min(c1Len, c2Len)
        c_max = max(c1Len, c2Len)
        g_prime_d = CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen)
        g_d = CalcGd(d, c1Len, c2Len, readLen)
        #a = stdDev ** 2 * g_prime_d + mean * g_d
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


def a(x,d, mu, sigma, c1Len, c2Len, r):
    """ integral x min(x-d-2 r+1, a+b+d-x+1, c-r+1) e^(-(x-mu)^2/(2 sigma^2)) dx = piecewise 

    if a+b-c+d+r-x>=0 && c+d+r-x<=0:
        (c-r+1) ((sqrt(pi/2) mu sigma erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma)))/sqrt(log(e))-(sigma^2 e^(-(x-mu)^2/(2 sigma^2)))/(log(e))) 
        
    if a+b-c+d+r-x<0&& a+b+2 d+2r-2 x<=0:
        -(sigma e^(-(x-mu)^2/(2 sigma^2)) (sqrt(2 pi) e^((x-mu)^2/(2 sigma^2)) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma)) (mu log(e) (-a-b-d+mu-1)+sigma^2)+2 sigma sqrt(log(e)) (a+b+d-mu-x+1)))/(2 log^(3/2)(e)) 

    otherwise:
        (sigma e^(-(x-mu)^2/(2 sigma^2)) (sqrt(2 pi) e^((x-mu)^2/(2 sigma^2)) (mu log(e) (-d+mu-2 r+1)+sigma^2) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma))-2 sigma sqrt(log(e)) (-d+mu-2 r+x+1)))/(2 log^(3/2)(e)) 

    """
    ## specify integration points of interest
    # place1 = d - 2*readLen + 1
    # place2 = min(c1Len,c2_len) - readLen + 1 
    # place3 = max(c1Len,c2_len) + d + readLen
    # place4 = c1Len + c2Len + d + 1
    c_max = max(c1Len,c2Len)
    c_min = min(c1Len,c2Len)
    if c_max + d + r - x >=0 and c_min + d + r - x <=0:
        a_term = \
        (1.0/sqrt(2*pi*sigma**2))* (c_min-r+1) * \
        (  sqrt(pi/2.0)*mu*sigma*normal.erf((x-mu)/(sqrt(2)*sigma)) \
        - sigma**2.0 *exp(-(x-mu)**2.0/(2*sigma**2.0)) )

    elif c_max + d + r - x < 0 and c1Len + c2Len + 2*d + 2*r - 2*x <=0:
        a_term = \
        - (1.0/sqrt(2*pi*sigma**2))*(1/2.0) * sigma* exp(-(x-mu)**2.0/(2.0*sigma**2.0)) * \
        ( sqrt(2*pi) *exp( (x-mu)**2 / (2*sigma**2) ) * normal.erf( (x-mu)/(sqrt(2)*sigma)) * \
        (mu * (-c1Len - c2Len - d + mu - 1) + sigma**2.0) \
        + 2*sigma* (c1Len + c2Len + d - mu - x + 1) )

    else:
        a_term = (1.0/sqrt(2*pi*sigma**2))*(1/2.0) * (sigma * exp(-(x-mu)**2/(2*sigma**2))) *  (sqrt(2*pi)*exp((x-mu)**2/(2*sigma**2))*(mu*(-d+mu-2*r+1)+sigma**2)*normal.erf((x-mu)/(sqrt(2)*sigma))- 2*sigma*(-d+mu-2*r+x+1)) 
    return a_term



def b(x, d, mu, sigma, c1Len, c2Len, r):
    """
    integral (x^2 min(x-d-2 r+1, a+b+d-x+1, c-r+1))/e^((x-mu)^2/(2 sigma^2)) dx = piecewise
    
    if  a+b-c+d+r-x>=0 && c+d+r-x<=0 :
        (sigma (c-r+1) e^(-(mu^2+x^2)/(2 sigma^2)) (sqrt(2 pi) (mu^2 log(e)+sigma^2) e^((mu^2+x^2)/(2 sigma^2)) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma))-2 sigma sqrt(log(e)) (mu+x) e^((mu x)/sigma^2)))/(2 log^(3/2)(e)) 
    
    elif a+b-c+d+r-x<0 && a+b+2 d+2r-2 x<=0:
        (sigma e^(-(mu^2+x^2)/(2 sigma^2)) (2 sigma e^((mu x)/sigma^2) (log(e) (mu (-a-b-d+mu-1)-x (a+b+d-mu+1)+x^2)+2 sigma^2)-sqrt(2 pi) sqrt(log(e)) e^((mu^2+x^2)/(2 sigma^2)) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma)) (mu^2 log(e) (-a-b-d+mu-1)-sigma^2 (a+b+d-3 mu+1))))/(2 log^2(e)) 
    
    else:
        (sigma e^(-(mu^2+x^2)/(2 sigma^2)) (sqrt(2 pi) sqrt(log(e)) e^((mu^2+x^2)/(2 sigma^2)) (mu^2 log(e) (-d+mu-2 r+1)-sigma^2 (d-3 mu+2 r-1)) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma))-2 sigma e^((mu x)/sigma^2) (log(e) (mu (-d+mu-2 r+1)+x (-d+mu-2 r+1)+x^2)+2 sigma^2)))/(2 log^2(e)) 
    """
    c_max = max(c1Len,c2Len)
    c_min = min(c1Len,c2Len)
    if c_max + d + r - x >=0 and c_min + d + r - x <=0:
        b_term = \
        (1.0/sqrt(2*pi*sigma**2)) *(1/2.0)*sigma*(c_min-r+1)* exp(-(mu**2+x**2)/(2*sigma**2)) *\
        (sqrt(2*pi)*(mu**2 + sigma**2) * exp((mu**2+x**2)/(2*sigma**2))*normal.erf( (x-mu)/(sqrt(2)*sigma)) \
        - 2*sigma*(mu+x)*exp((mu*x)/sigma**2) )
   
    elif c_max + d + r - x < 0 and c1Len + c2Len + 2*d + 2*r - 2*x <=0:
        b_term = \
        (1.0/sqrt(2*pi*sigma**2))*(1/2.0)*sigma *  exp(-(mu**2+x**2)/(2*sigma**2)) *\
        (2*sigma*exp((mu*x)/sigma**2)*((mu *(-c1Len- c2Len -d+mu-1)- x*(c1Len+ c2Len +d-mu+1) + \
        x**2)+2*sigma**2) - sqrt(2*pi) *exp((mu**2+x**2)/(2*sigma**2))* \
        normal.erf((x-mu)/(sqrt(2)*sigma)) *(mu**2 *(-c1Len- c2Len -d+mu-1) - sigma**2 * (c1Len+ c2Len +d-3*mu+1))) 
    
    else:
        b_term = \
        (1.0/sqrt(2*pi*sigma**2))*(1/2.0)*sigma *  exp(-(mu**2+x**2)/(2*sigma**2)) *\
        (sqrt(2*pi)*exp((mu**2+x**2)/(2*sigma**2)) *  (mu**2 * (-d+mu-2*r+1) - sigma**2 *(d-3*mu+ \
        2*r-1))*normal.erf((x-mu)/(sqrt(2.0)*sigma))-2*sigma*exp((mu*x)/sigma**2)* \
        ( ( mu*(-d+mu-2*r+1)+ x*(-d+mu-2*r+1)+x**2)+2*sigma**2))
    return b_term


def g_d(x, d, mu, sigma, c1Len, c2Len, r):
    """
        integral min(x-d-2 r+1, a+b+d-x+1, c-r+1) e^(-(x-mu)^2/(2 sigma^2)) dx = piecewise 
        if  a+b-c+d+r-x>=0&&c+d+r-x<=0:
            (sqrt(pi/2) sigma (c-r+1) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma)))/sqrt(log(e)) 
        elif  a+b-c+d+r-x<0&&a+b+2 d+2 r-2 x<=0:

            (sqrt(pi/2) sigma (a+b+d-mu+1) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma)))/sqrt(log(e))+(sigma^2 e^(-(x-mu)^2/(2 sigma^2)))/(log(e))
        else:
            -(sqrt(pi/2) sigma (d-mu+2 r-1) erf((sqrt(log(e)) (x-mu))/(sqrt(2) sigma)))/sqrt(log(e))-(sigma^2 e^(-(x-mu)^2/(2 sigma^2)))/(log(e))
    """
    c_max = max(c1Len,c2Len)
    c_min = min(c1Len,c2Len)

    if c_max + d + r - x >=0 and c_min + d + r - x <=0:
        print 'g_d1'
        g_d_term = \
        sqrt(pi/2.0) * (1.0/sqrt(2*pi*sigma**2)) * sigma *(c_min-r+1)*normal.erf( (x-mu)/(sqrt(2.0)*sigma))
    
    elif c_max + d + r - x < 0 and c1Len + c2Len + 2*d + 2*r - 2*x <=0:
        print 'g_d2'
        g_d_term = \
        (1.0/sqrt(2*pi*sigma**2)) * \
        ( sqrt(pi/2.0) * sigma *(c1Len+c2Len+d-mu+1)* normal.erf((x-mu)/(sqrt(2.0)*sigma)) +\
        sigma**2* exp(-(x-mu)**2/(2*sigma**2)))
    else:
        print 'g_d3'
        g_d_term = \
        -(1.0/sqrt(2*pi*sigma**2)) * \
        ( sqrt(pi/2.0) * sigma *(d-mu+2*r-1)* normal.erf((x-mu)/(sqrt(2.0)*sigma)) -\
        sigma**2* exp(-(x-mu)**2/(2*sigma**2)))

    
    return g_d_term

def stddev_given_d(mean, stdDev, readLen, c1Len, c2Len, d):
    """
        b/g_d - (a/g_d)^2
    """
    if c1Len + c2Len + d + 1 > mean + 5*stdDev:
        upper_integration_point = mean + 5*stdDev
    else:
        upper_integration_point = c1Len + c2Len + d + 1 

    lower_integration_point =  d + 2 * readLen - 1
    
    b_ = b(upper_integration_point, d, mean, stdDev, c1Len, c2Len, readLen) - b(lower_integration_point, d, mean, stdDev, c1Len, c2Len, readLen)
    a_ = a(upper_integration_point, d, mean, stdDev, c1Len, c2Len, readLen) - a(lower_integration_point, d, mean, stdDev, c1Len, c2Len, readLen)
    g_d_lower = g_d(lower_integration_point,d, mean, stdDev, c1Len, c2Len, readLen)
    g_d_upper = g_d(upper_integration_point,d, mean, stdDev, c1Len, c2Len, readLen)
    g_d_area = g_d_upper - g_d_lower
    print g_d_lower,g_d_upper
    print a_, g_d_area, a_/ g_d_area,  ( a_ / g_d_area )**2  
    print b_, g_d_area, b_/ g_d_area, b_ / g_d_area

    var = b_ / g_d_area\
    - ( a_ / g_d_area )**2

    print var
    if var < 0:
        sys.stderr.write('Error in estimate, a={0}, b={1}'.format(a_,b_))
        std_dev_est = 0
    else:
        std_dev_est = sqrt(var)  
    return std_dev_est


