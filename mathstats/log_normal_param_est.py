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
	val = sum([ (mu - sigma**2 - math.log(o+d)) / (sigma**2 * (o + d)) for o in observations])
	print 'VAL:', val
	return val


def math_log_if_pos(x):
	if x > 0:
		return math.log(x)
	else:
		print 'here'
		inf = float("inf")
		return -inf

# def calc_g_prim_d_OLD(d, mu, sigma, c_min, c_max, c1Len, c2Len, r, cutoff_approx_level):

# 	x1 = d+2*r if d+2*r > 0 else 0
# 	x1 = d+2*r if d+2*r < cutoff_approx_level else cutoff_approx_level

# 	x2 = c_min + d + r - 1 if c_min + d + r - 1 > 0 else 0
# 	x2 = c_min + d + r - 1 if c_min + d + r - 1 < cutoff_approx_level else cutoff_approx_level

# 	x3 = c_min + d + r if c_min + d + r > 0 else 0
# 	x3 = c_min + d + r if c_min + d + r < cutoff_approx_level else cutoff_approx_level

# 	x4 = c_max + d + r if c_max + d + r > 0 else 0
# 	x4 = c_max + d + r if c_max + d + r < cutoff_approx_level else cutoff_approx_level

# 	x5 = c_max + d + r + 1 if c_max + d + r + 1 > 0 else 0
# 	x5 = c_max + d + r + 1 if c_max + d + r + 1 < cutoff_approx_level else cutoff_approx_level

# 	x6 = c_min + c_max + d if c_min + c_max + d < cutoff_approx_level else cutoff_approx_level

# 	term1 = -0.5 * (normal.erf(( math_log_if_pos(x6 + 1) - mu) / (2.0 ** 0.5 * float(sigma))) + normal.erf((math_log_if_pos(x1 - 1) - mu) / (2.0 ** 0.5 * float(sigma))))
# 	term2 = +0.5 * (normal.erf(( math_log_if_pos(x4) - mu) / (2.0 ** 0.5 * float(sigma))) + normal.erf(( math_log_if_pos(x3) - mu) / (2.0 ** 0.5 * float(sigma))))
# 	g_prime_d = term1 + term2
# 	return -g_prime_d

def calc_g_prim_d(d, mu, sigma, c_min, c_max, r, cutoff_approx_level):

	x1 = d+2*r if d+2*r > 0 else 0
	x1 = d+2*r if d+2*r < cutoff_approx_level else cutoff_approx_level

	x2 = c_min + d + r - 1 if c_min + d + r - 1 > 0 else 0
	x2 = c_min + d + r - 1 if c_min + d + r - 1 < cutoff_approx_level else cutoff_approx_level

	x3 = c_min + d + r if c_min + d + r > 0 else 0
	x3 = c_min + d + r if c_min + d + r < cutoff_approx_level else cutoff_approx_level

	x4 = c_max + d + r if c_max + d + r > 0 else 0
	x4 = c_max + d + r if c_max + d + r < cutoff_approx_level else cutoff_approx_level

	x5 = c_max + d + r + 1 if c_max + d + r + 1 > 0 else 0
	x5 = c_max + d + r + 1 if c_max + d + r + 1 < cutoff_approx_level else cutoff_approx_level

	x6 = c_min + c_max + d if c_min + c_max + d < cutoff_approx_level else cutoff_approx_level

	part_1 = g_prim_d(x2, d, mu, sigma, c_min, c_max, r) - g_prim_d(x1, d, mu, sigma, c_min, c_max, r)

	#part_2 = from x = c_min + d + r to x = c_max + d + r
	part_2 =  g_prim_d(x4, d, mu, sigma, c_min, c_max, r) - g_prim_d(x3, d, mu, sigma, c_min, c_max, r)

	#part_3 = from x = c_max + d + r + 1 to c_min + c_max + d 
	part_3 = g_prim_d(x6, d, mu, sigma, c_min, c_max, r) - g_prim_d(x5, d, mu, sigma, c_min, c_max, r)

	return part_1 + part_2 + part_3

def g_prim_d(x, d, mu, sigma, c_min, c_max, r):
	# in the middle part (roof of triangle house, between c_min and c_max)
	if c_max+d+r-x >= 0 and c_min+d+r-x <= 0:
		return 0

	# x is large, on the down slope of the triangle house
	elif c_max+d+r-x < 0 and c_min + c_max + d - x > -1:
		#print 'loolz', normal.erf((mu-math_log_if_pos(x))+sigma**2/(sqrt(2) * sigma)), x #normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma))
		return 0.5 * normal.erf((math_log_if_pos(x) - mu)/(sqrt(2) * sigma))

	# x is small, beginning of triangle house
	elif d+2*r-x < 1 and c_min+d+r-x > 0:
		#print 'laal', normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma))
		return -0.5 * normal.erf((math_log_if_pos(x) - mu)/(sqrt(2) * sigma))
	else:
		return 0

def gd(x, d, mu, sigma, c_min, c_max, r):
	# in the middle part (roof of triangle house, between c_min and c_max)
	if c_max+d+r-x >= 0 and c_min+d+r-x <= 0:
		return (c_min-r+1)/2.0 * normal.erf((math_log_if_pos(x)-mu)/(math.sqrt(2.0)*sigma))

	# x is large, on the down slope of the triangle house
	elif c_max+d+r-x < 0 and c_min + c_max + d - x > -1:
		#print 'loolz', normal.erf((mu-math_log_if_pos(x))+sigma**2/(sqrt(2) * sigma)), x #normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma))
		return 1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf((mu-math_log_if_pos(x))+sigma**2/(sqrt(2) * sigma)) - (c_min+c_max+d+1)*normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma)))

	# x is small, beginning of triangle house
	elif d+2*r-x < 1 and c_min+d+r-x > 0:
		#print 'laal', normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma))
		return -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x)) / (math.sqrt(2.0)*sigma)))

	else:
		return 0


def calc_gd(d, mu, sigma, c_min, c_max, r, cutoff_approx_level):

	"""
		4 different itegration points of the trapetzoid house
		from 1 possible placement.
		Values can get shaky too far out on the log normal distribution due to float aritmethic
		cutoff_approx_level is the upper point where estimations are still stable
	"""
	print 'CMON', d+2*r, c_min + d + r, c_max + d + r, c_min + c_max + d

	x1 = d+2*r if d+2*r > 0 else 0
	x1 = d+2*r if d+2*r < cutoff_approx_level else cutoff_approx_level

	x2 = c_min + d + r - 1 if c_min + d + r - 1 > 0 else 0
	x2 = c_min + d + r - 1 if c_min + d + r - 1 < cutoff_approx_level else cutoff_approx_level

	x3 = c_min + d + r if c_min + d + r > 0 else 0
	x3 = c_min + d + r if c_min + d + r < cutoff_approx_level else cutoff_approx_level

	x4 = c_max + d + r if c_max + d + r > 0 else 0
	x4 = c_max + d + r if c_max + d + r < cutoff_approx_level else cutoff_approx_level

	x5 = c_max + d + r + 1 if c_max + d + r + 1 > 0 else 0
	x5 = c_max + d + r + 1 if c_max + d + r + 1 < cutoff_approx_level else cutoff_approx_level

	x6 = c_min + c_max + d if c_min + c_max + d < cutoff_approx_level else cutoff_approx_level

	print 'JAO', x1, x2 
	print 'x1:', -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x1))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma)))
	print 'x2:', -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x2))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma)))
	print 'x1 first:',x1,(  normal.erf(((mu-math_log_if_pos(x1))+sigma**2)/(math.sqrt(2.0)*sigma)))
	print 'x2 first:',x2, (  normal.erf(((mu-math_log_if_pos(x2))+sigma**2)/(math.sqrt(2.0)*sigma)))
	print 'x1 second:',x1, -1/2.0 * ( (d+2*r-1)*normal.erf((mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma))), math_log_if_pos(x1), (mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma), math.erf((mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma))
	print 'x2 second:',x2, -1/2.0 * ( (d+2*r-1)*normal.erf((mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma))), math_log_if_pos(x2), (mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma), math.erf((mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma))
	#part_1 =  from x = d+2*r to x = c_min + d + r
	part_1 = gd(x2, d, mu, sigma, c_min, c_max, r) - gd(x1, d, mu, sigma, c_min, c_max, r)
	print gd(c_min + d + r - 1, d, mu, sigma, c_min, c_max, r), gd(d + 2*r, d, mu, sigma, c_min, c_max, r)
	print 'LOL', gd(x2, d, mu, sigma, c_min, c_max, r), gd(x1, d, mu, sigma, c_min, c_max, r)

	#part_2 = from x = c_min + d + r to x = c_max + d + r
	part_2 =  gd(x4, d, mu, sigma, c_min, c_max, r) - gd(x3, d, mu, sigma, c_min, c_max, r)

	print gd(c_max + d + r, d, mu, sigma, c_min, c_max, r), gd(c_min + d + r, d, mu, sigma, c_min, c_max, r)

	#part_3 = from x = c_max + d + r + 1 to c_min + c_max + d 
	part_3 = gd(x6, d, mu, sigma, c_min, c_max, r) - gd(x5, d, mu, sigma, c_min, c_max, r)

	print gd(c_min + c_max + d, d, mu, sigma, c_min, c_max, r), gd(c_max + d + r + 1, d, mu, sigma, c_min, c_max, r)
	print 'from to', c_min + c_max + d, c_max + d + r + 1
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
	mu, sigma, r, c1Len, c2Len = float(mu), float(sigma), float(r), float(c1Len), float(c2Len)
	# Some pruning of extreme observations would be good to have here
	min_obs = min(observations)
	max_obs = max(observations)

	n = len(observations)
	d_ML = 'lol'
	c_min = min(c1Len, c2Len)
	c_max = max(c1Len, c2Len)
	min_gap = -int(min_obs) + int(2*r)

	# set max to mean + 2*sd
	cutoff_approx_level = math.exp(mu + 0.5*sigma**2) + 5000*math.sqrt( (math.exp(sigma**2) - 1) * math.exp(2*mu + sigma**2) ) 
	upper_quantile = math.exp(mu + 0.5*sigma**2) + 3*math.sqrt( (math.exp(sigma**2) - 1) * math.exp(2*mu + sigma**2) )
	max_gap = int(upper_quantile) if upper_quantile + max_obs <  cutoff_approx_level else int(cutoff_approx_level - max_obs)  #emperical_max
	d_lower = min_gap
	d_upper = max_gap
	print d_upper, cutoff_approx_level

	# always choose the lowest root (0) to the function --> based on emperical observations
	# IF MEAN HIGHER THAN CONDITIONAL MEAN AND STDDEV LOW: Search for negative gap
	# IF MEAN HIGHER THAN CONDITIONAL MEAN AND STDDEV high (few observarions): Search for positive gap
	# IF MEAN HIGHER THAN CONDITIONAL MEAN AND STDDEV high (few observarions): Search for positive gap

	# initial values
	# use newton rhapson
	observation_term = calc_log_norm(observations, d_lower, mu, sigma)
	g_d = calc_gd(d_lower, mu, sigma, c_min, c_max, r, cutoff_approx_level)
	g_prime_d = calc_g_prim_d(d_lower, mu, sigma, c_min, c_max, r, cutoff_approx_level)
	g_d_ratio = n*g_prime_d/g_d
	print 'gap:',d_lower, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',observation_term
	# d = d_lower
	# while g_d_ratio < observation_term:
	# 	print 'CURRENT GAP:', d
	# 	d = 
	# 	other_term = calc_log_norm(observations, d_ML, mu, sigma)
	# 	g_d = calc_gd(d_ML, mu, sigma, c_min, c_max, r, cutoff_approx_level)
	# 	g_prime_d = calc_g_prim_d(d_ML, mu, sigma, c_min, c_max, r, cutoff_approx_level)
	# 	g_d_ratio = n*g_prime_d/g_d
	# 	print 'gap:',d_ML, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',other_term

	# 	if g_d_ratio > other_term:
	# 		d_upper = d_ML
	# 	else:
	# 		d_lower = d_ML

	# d_ML = (d_upper + d_lower) / 2.0

	# print d_ML
	# a = raw_input("press enter")

	##################################

	# while d_upper - d_lower > 0.1:
	# 	d_ML = (d_upper + d_lower) / 2.0
	# 	d_ML_2 = d_ML + 5
	# 	print 'D_MAX:', d_upper, "D_MIN:", d_lower
	# 	print 'CURRENT GAP:', d_ML
	# 	other_term = calc_log_norm(observations, d_ML, mu, sigma)
	# 	g_d = calc_gd(d_ML, mu, sigma, c_min, c_max, r, cutoff_approx_level)
	# 	g_prime_d = calc_g_prim_d(d_ML, mu, sigma, c_min, c_max, r, cutoff_approx_level)
	# 	g_d_ratio = n*g_prime_d/g_d
	# 	print 'gap:',d_ML, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',other_term

	# 	if g_d_ratio > other_term:
	# 		d_upper = d_ML
	# 	else:
	# 		d_lower = d_ML

	# d_ML = (d_upper + d_lower) / 2.0

	# print d_ML
	# a = raw_input("press enter")

	x = []
	y = []
	y_1 = []
	y_2 = []
	g_d_list = []
	g_prim_d_list = []

	for d_ in range(min_gap, max_gap, 100):
		#for frac in range(10):
		d = d_ # + frac/10.0
		#print "GAP", d
		g_d = calc_gd(d, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_d_list.append(g_d)
		g_prime_d = calc_g_prim_d(d, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_prim_d_list.append(g_prime_d)

		g_d_ratio = n*g_prime_d/g_d
		other_term = calc_log_norm(observations, d, mu, sigma)
		x.append(d_)
		y.append(g_d_ratio - other_term)
		y_1.append(g_d_ratio)
		y_2.append(other_term)

		# if 0.95 < g_d_ratio /calc_log_norm(observations, d, mu, sigma) < 1.05:
		# 	d_ML = d
		# 	print 'MLGAP', d


		print 'gap:',d, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',other_term
	import matplotlib.pyplot as plt
	# plt.plot(x,y,'-')
	# x1,x2,y1,y2 = plt.axis()
	# # plt.axis((x1,x2,-0.01,0.001))
	# plt.show()

	plt.plot(x,y_1,'-')
	plt.plot(x,y_2,'-r')
	x1,x2,y1,y2 = plt.axis()
	#plt.axis((x1,x2,-0.01,0.001))
	plt.show()

	# plt.plot(x,g_prim_d_list,'-')
	# x1,x2,y1,y2 = plt.axis()
	# # plt.axis((x1,x2,-0.01,0.001))
	# plt.show()

	# plt.plot(x,g_d_list,'-r')
	# x1,x2,y1,y2 = plt.axis()
	# # plt.axis((x1,x2,-0.01,0.001))
	# plt.show()

	print 'MLGAP', d_ML
	print d_upper, cutoff_approx_level


# formula from: http://www.wolframalpha.com/input/?i=integral+of+max%280%2C+min%28x-d-2r%2B1%2Ca%2Bb%2Bd-x%2B1%2Cc-r%2B1%29%29*%281%2F%28x*sigma*sqrt%282*pi%29%29%29*+e%5E-%28%28ln%28x%29+-+mu%29%5E2%2F%282*sigma%5E2%29%29+with+respect+to+x

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


