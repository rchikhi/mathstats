'''
    Estimates parameters of a truncated lognormal distribution given observations. 
    This truncated lognormal distribution has a distribution on a similar form given in 
    http://www.ncbi.nlm.nih.gov/pubmed/22923455, but with f(x) lognormal instead of normal. 
    
    Parameters of a truncated lognormal distribution are:
    * mean
    * variance
    * minimum value (lower bound 'a') 
    * maximum value (upper bound 'b')
    
    Application:
    Calculates gap size between two contigs assuming lognormal distribution of library.
    '''

import sys
import math
from math import exp, pi, sqrt
from mathstats.normaldist import normal

def calc_log_norm(observations, d, mu, sigma):
	val = sum([ ( mu - sigma**2 - math.log(o+d) ) / (sigma**2 * (o + d)) for o in observations])
	return val


def math_log_if_pos(x):
	if x > 0:
		return math.log(x)
	else:
		print 'here'
		inf = float("inf")
		return inf

def calc_g_of_d_ratio(d, mu, sigma, c1Len, c2Len, readLen):
    c_min = min(c1Len, c2Len)
    c_max = max(c1Len, c2Len)

    def Nominator(d, c_min, c_max, c1Len, c2Len, readLen):
        term1 = -0.5 * (normal.erf(( math_log_if_pos(c_min + c_max + d + 1) - mu) / (2 ** 0.5 * float(sigma))) + normal.erf((math_log_if_pos(d + 2 * readLen - 1) - mu) / (2 ** 0.5 * float(sigma))))
        term2 = +0.5 * (normal.erf(( math_log_if_pos(d + c_max + readLen) - mu) / (2 ** 0.5 * float(sigma))) + normal.erf(( math_log_if_pos(d + c_min + readLen) - mu) / (2 ** 0.5 * float(sigma))))
        g_prime_d = term1 + term2
        return -g_prime_d

    def Denominator(d, c1Len, c2Len, readLen):
    	#print 'lol', d, c1Len, c2Len, readLen
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...

        term1 = (c_min - readLen + 1) / 2.0 * (normal.erf(( math_log_if_pos(c_max + d + readLen) - mu) / ((2 ** 0.5) * sigma)) - normal.erf(( math_log_if_pos(c_min + d + readLen) - mu) / ((2 ** 0.5) * sigma)))

        term2 = (c_min + c_max + d - mu + 1) / 2.0 * (normal.erf(( math_log_if_pos(c_min + c_max + d + 1) - mu) / (2 ** 0.5 * float(sigma))) - normal.erf(( math_log_if_pos(c_max + readLen + d) - mu) / (2 ** 0.5 * float(sigma))))

        term3 = (d + 2 * readLen - mu - 1) / 2.0 * (normal.erf(( math_log_if_pos(d + 2 * readLen - 1) - mu) / (2 ** 0.5 * float(sigma))) - normal.erf(( math_log_if_pos(c_min + d + readLen) - mu) / ((2 ** 0.5) * sigma)))

        term4 = sigma / ((2 * pi) ** 0.5) * (exp(-(( math_log_if_pos(c_min + c_max + d + 1) - mu) ** 2) / (float(2 * sigma ** 2))) + exp(-((math_log_if_pos(d + 2 * readLen - 1) - mu) ** 2) / (float(2 * sigma ** 2))))

        term5 = -sigma / ((2 * pi) ** 0.5) * (exp(-((math_log_if_pos(c_max + readLen + d) - mu) ** 2) / (float(2 * sigma ** 2))) + exp(-(( math_log_if_pos(c_min + readLen + d) - mu) ** 2) / (float(2 * sigma ** 2))))
        g_d = term1 + term2 + term3 + term4 + term5
        print term1, term2,  term3,  term4,  term5
        return g_d

    denominator = Denominator(d, c1Len, c2Len, readLen)
    nominator = Nominator(d, c_min, c_max, c1Len, c2Len, readLen)
    print denominator, nominator

    try:
        g_of_d_ratio = nominator / denominator
    except ZeroDivisionError:
        g_of_d_ratio = 0.0000000001

    return g_of_d_ratio


def calc_g_prim_d(d, mu, sigma, c_min, c_max, c1Len, c2Len, r):
    term1 = -0.5 * (normal.erf(( math_log_if_pos(c_min + c_max + d + 1) - mu) / (2 ** 0.5 * float(sigma))) + normal.erf((math_log_if_pos(d + 2 * r - 1) - mu) / (2 ** 0.5 * float(sigma))))
    term2 = +0.5 * (normal.erf(( math_log_if_pos(d + c_max + r) - mu) / (2 ** 0.5 * float(sigma))) + normal.erf(( math_log_if_pos(d + c_min + r) - mu) / (2 ** 0.5 * float(sigma))))
    g_prime_d = term1 + term2
    return g_prime_d

def gd(x, d, mu, sigma, c_min, c_max, r):
	#x = o + d
	
	# in the middle part (roof of triangle house, between c_min and c_max)
	if c_max+d+r-x >= 0 and c_min+d+r-x <= 0:
		return ((c_min-r+1)/2.0 * normal.erf((math_log_if_pos(x)-mu))/(math.sqrt(2.0)*sigma))

	# x is large, on the down slope of the triangle house
	elif c_max+d+r-x < 0 and c_min + c_max + d - x > -1:
		return 1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf((mu-math_log_if_pos(x))+sigma**2)/(sqrt(2) *sigma) - (c_min+c_max+d+1)*normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma)))

	# x is small, beginning of triangle house
	elif  d+2*r-x < 1 and c_min+d+r-x > 0:
		return -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x)) / (math.sqrt(2.0)*sigma)))

	else:
		return 0


def calc_gd(d, mu, sigma, c_min, c_max, r):
	#part_1 =  from x = d+2*r to x = c_min + d + r
	part_1 = gd(c_min + d + r - 1, d, mu, sigma, c_min, c_max, r) - gd(d+2*r, d, mu, sigma, c_min, c_max, r)
	print c_min + d + r - 1, d + 2*r
	print gd(c_min + d + r - 1, d, mu, sigma, c_min, c_max, r), gd(d + 2*r, d, mu, sigma, c_min, c_max, r)

	#part_2 = from x = c_min + d + r to x = c_max + d + r
	part_2 =  gd(c_max + d + r, d, mu, sigma, c_min, c_max, r) - gd(c_min + d + r, d, mu, sigma, c_min, c_max, r)

	print gd(c_max + d + r, d, mu, sigma, c_min, c_max, r), gd(c_min + d + r, d, mu, sigma, c_min, c_max, r)

	#part_3 = from x = c_max + d + r + 1 to c_min + c_max + d 
	part_3 = gd(c_min + c_max + d, d, mu, sigma, c_min, c_max, r) - gd(c_max + d + r + 1, d, mu, sigma, c_min, c_max, r)
	print gd(c_min + c_max + d, d, mu, sigma, c_min, c_max, r), gd(c_max + d + r + 1, d, mu, sigma, c_min, c_max, r)

	print part_1
	print part_2
	print part_3
	return part_1 + part_2 + part_3

def GapEstimator(mean, sigma, read_length, observations, c1_len, c2_len=None, emperical_max=None):
    '''
    Calculates the lower bound (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455). The upper bound is then
    uniquely determined by the lower bound.
    '''

    if c2_len is None:
        c2_len = c1_len
    if emperical_max is None:
        emperical_max = 3*math.exp(mean)

    d_ML = CalcMLvaluesOfdGeneral(mean, sigma, read_length, c1_len, observations, c2_len, emperical_max)

    return d_ML


def CalcMLvaluesOfdGeneral(mu, sigma, r, c1Len, observations, c2Len, emperical_max):
    """
        returns the ML-gap as float. The ML gap is searched
        for with an accuracy of 0.1bp. This is NOT the accuracy of the ML-estimate
        however since it depends on the number of samples etc.
    """


    min_obs = min(observations)
    max_obs = max(observations)
    n = len(observations)
    d_ML = 'lol'
    c_min = min(c1Len, c2Len)
    c_max = max(c1Len, c2Len)
    min_gap = -int(min_obs) + int(2*r)
    for d_ in range(min_gap, emperical_max, 100):
		#for frac in range(10):
		d = d_ # + frac/10.0
		g_d = calc_gd(d, mu, sigma, c_min, c_max, r)
		g_prime_d = calc_g_prim_d(d, mu, sigma, c_min, c_max, c1Len, c2Len, r)

		g_d_ratio = n*g_prime_d/g_d
		other_term = calc_log_norm(observations, d, mu, sigma)
		if 0.95 < g_d_ratio /calc_log_norm(observations, d, mu, sigma) < 1.05:
			d_ML = d
			print 'MLGAP', d


		print 'gap:',d, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',other_term

    print 'MLGAP', d_ML 


    # min_obs = min(observations)
    # max_obs = max(observations)
    # n = len(observations)
    # d_ML = 'lol'
    # for d in range(-int(min_obs) + int(2*readLen) , emperical_max, 100):
    # 	g_d = n* calc_g_of_d_ratio(d, mu, sigma, c1Len, c2Len, readLen)
    # 	print 'gap:',d, 'g(d)*n', g_d , 'other:',calc_log_norm(observations, d, mu, sigma)
    # 	if 0.95 < n* calc_g_of_d_ratio(d, mu, sigma, c1Len, c2Len, readLen) / -calc_log_norm(observations, d, mu, sigma) < 1.05:
    # 		d_ML = d
    # 		print 'MLGAP', d
    # print 'MLGAP', d_ML 



   #  d_upper = emperical_max
   #  d_lower = int(-c1Len - c2Len) + int(2 * readLen)
   #  print d_upper, d_lower
   #  while d_upper - d_lower > 0.1:
   #      d_ML = (d_upper + d_lower) / 2.0
   #      print d_ML
   #      if c2Len:
			# g_of_d_ratio = calc_g_of_d_ratio(d_ML, mu, sigma, c1Len, c2Len, readLen)
			# log_norm_sum = calc_log_norm(observations, d_ML, mu, sigma)
			# print "g(d)*n",g_of_d_ratio*n
			# print "other side:", log_norm_sum
   #      # else:
   #      #     func_of_d = ML_one_reference(d_ML, mu, sigma, c1Len, readLen)

   #      # check which way it is!
   #      if n*g_of_d_ratio < log_norm_sum:
   #          d_upper = d_ML
   #      else:
   #          d_lower = d_ML

   #  d_ML = (d_upper + d_lower) / 2.0
   #  return d_ML



# evaluate this integral in 4 different points

# integral (max(0, min(x-d-2 r+1, a+b+d-x+1, c-r+1)))/((x sigma sqrt(2 pi)) e^((log(x)-mu)^2/(2 sigma^2))) dx =

# # in the middle part (roof of triangle house, between c_min and c_max)
# if c_max+d+r-x>=0 and c_min-r>-1 and c_min+d+r-x<=0:
# ((c-r+1) erf((sqrt(log(e)) (log(x)-mu))/(sqrt(2) sigma)))/(2 sqrt(log(e))) 

# # x is large, on the down slope of the triangle house
# if c_max+d+r-x<0 and a+b+2*d+2*r-2*x<=0 and a+b+d-x > -1:
# (e^(sigma^2/(2 log(e))+mu) erf((log(e) (mu-log(x))+sigma^2)/(sqrt(2) sigma sqrt(log(e))))-(a+b+d+1) erf((sqrt(log(e)) (mu-log(x)))/(sqrt(2) sigma)))/(2 sqrt(log(e)))

# # x is small, biginning of triangle house
# if  d+2*r-x<1 and a+b+2*d+2*r-2*x > 0 and c_min+d+r-x>0:
# -(e^(sigma^2/(2 log(e))+mu) erf((log(e) (mu-log(x))+sigma^2)/(sqrt(2) sigma sqrt(log(e))))-(d+2 r-1) erf((sqrt(log(e)) (mu-log(x)))/(sqrt(2) sigma)))/(2 sqrt(log(e))) 


