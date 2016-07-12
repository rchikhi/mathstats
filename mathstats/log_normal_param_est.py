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
	return val


def math_log_if_pos(x):
	if x > 0:
		return math.log(x)
	else:
		#print 'here'
		inf = float("inf")
		return -inf

# def calc_g_prim_d_OLD(d, mu, sigma, c_min, c_max, c_min, c2_len, r, cutoff_approx_level):

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
	# print("x1", x1)
	x2 = c_min + d + r - 1 if c_min + d + r - 1 > 0 else 0
	x2 = c_min + d + r - 1 if c_min + d + r - 1 < cutoff_approx_level else cutoff_approx_level
	# print("x2",x2)

	x3 = c_min + d + r if c_min + d + r > 0 else 0
	x3 = c_min + d + r if c_min + d + r < cutoff_approx_level else cutoff_approx_level
	# print("x3",x3)

	x4 = c_max + d + r if c_max + d + r > 0 else 0
	x4 = c_max + d + r if c_max + d + r < cutoff_approx_level else cutoff_approx_level
	# print("x4", x4)

	x5 = c_max + d + r + 1 if c_max + d + r + 1 > 0 else 0
	x5 = c_max + d + r + 1 if c_max + d + r + 1 < cutoff_approx_level else cutoff_approx_level
	# print("x5",x5)

	x6 = c_min + c_max + d if c_min + c_max + d < cutoff_approx_level else cutoff_approx_level
	# print("x6",x6)

	part_1 = g_prim_d(x2, d, mu, sigma, c_min, c_max, r) - g_prim_d(x1, d, mu, sigma, c_min, c_max, r)
	# print "g_prim_d, part_1:", part_1
	#part_2 = from x = c_min + d + r to x = c_max + d + r
	part_2 =  g_prim_d(x4, d, mu, sigma, c_min, c_max, r) - g_prim_d(x3, d, mu, sigma, c_min, c_max, r)
	# print "g_prim_d, part_2:", part_2
	#part_3 = from x = c_max + d + r + 1 to c_min + c_max + d 
	part_3 = g_prim_d(x6, d, mu, sigma, c_min, c_max, r) - g_prim_d(x5, d, mu, sigma, c_min, c_max, r)
	# print "g_prim_d, part_3:", part_3
	# print("RETURNING:", part_1 + part_2 + part_3 )
	return part_1 + part_2 + part_3

def g_prim_d(x, d, mu, sigma, c_min, c_max, r):
	# in the middle part (roof of triangle house, between c_min and c_max)
	if c_max+d+r-x >= 0 and c_min+d+r-x <= 0:
		# print("lol")
		return 0

	# x is large, on the down slope of the triangle house
	elif c_max+d+r-x < 0 and c_min + c_max + d - x > -1:
		#print 'loolz', normal.erf((mu-math_log_if_pos(x))+sigma**2/(sqrt(2) * sigma)), x #normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma))
		# print("ok")
		return 0.5 * normal.erf((math_log_if_pos(x) - mu)/(sqrt(2) * sigma))

	# x is small, beginning of triangle house
	elif d+2*r-x < 1 and c_min+d+r-x > 0:
		#print 'laal', normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma))
		return -0.5 * normal.erf((math_log_if_pos(x) - mu)/(sqrt(2) * sigma))
		# print("ok2")
	else:
		# print("lol2")
		return 0

def gd(x, d, mu, sigma, c_min, c_max, r):
	# in the middle part (roof of triangle house, between c_min and c_max)
	if c_max+d+r-x >= 0 and c_min+d+r-x <= 0:
		multiplier_height = max((c_min-r+1)/2.0, 1)
		return multiplier_height * normal.erf((math_log_if_pos(x)-mu)/(math.sqrt(2.0)*sigma))
		# return (c_min-r+1)/2.0 * normal.erf((math_log_if_pos(x)-mu)/(math.sqrt(2.0)*sigma))

	# x is large, on the down slope of the triangle house
	elif c_max+d+r-x < 0 and c_min + c_max + d - x > -1:
		#print 'loolz', normal.erf((mu-math_log_if_pos(x))+sigma**2/(sqrt(2) * sigma)), x #normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma))
		return 1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf((mu-math_log_if_pos(x)+sigma**2)/(sqrt(2) * sigma)) - (c_min+c_max+d+1)*normal.erf((mu-math_log_if_pos(x))/(math.sqrt(2.0)*sigma)))

	# x is small, beginning of triangle house
	elif d+2*r-x < 1 and c_min+d+r-x > 0:
		#print 'laal', normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma))
		return -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x)) / (math.sqrt(2.0)*sigma)))

	else:
		# print("OMG")
		return 0.0


def calc_gd(d, mu, sigma, c_min, c_max, r, cutoff_approx_level):

	"""
		4 different itegration points of the trapetzoid house
		from 1 possible placement.
		Values can get shaky too far out on the log normal distribution due to float aritmethic
		cutoff_approx_level is the upper point where estimations are still stable
	"""
	# print 'CMON', d+2*r, c_min + d + r, c_max + d + r, c_min + c_max + d

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
	# print("x1", x1)
	# print("x2", x2)
	# print("x3", x3)
	# print("x4", x4)
	# print("x5", x5)
	# print("x6", x6)

	# print x5, x6
	# print 'JAO', x1, x2 
	# print 'x1:', -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x1))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma)))
	# print 'x2:', -1/2.0 * (math.exp(sigma**2/2.0 + mu) * normal.erf(((mu-math_log_if_pos(x2))+sigma**2)/(math.sqrt(2.0)*sigma)) - (d+2*r-1)*normal.erf((mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma)))
	# print 'x1 first:',x1,(  normal.erf(((mu-math_log_if_pos(x1))+sigma**2)/(math.sqrt(2.0)*sigma)))
	# print 'x2 first:',x2, (  normal.erf(((mu-math_log_if_pos(x2))+sigma**2)/(math.sqrt(2.0)*sigma)))
	# print 'x1 second:',x1, -1/2.0 * ( (d+2*r-1)*normal.erf((mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma))), math_log_if_pos(x1), (mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma), math.erf((mu-math_log_if_pos(x1)) / (math.sqrt(2.0)*sigma))
	# print 'x2 second:',x2, -1/2.0 * ( (d+2*r-1)*normal.erf((mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma))), math_log_if_pos(x2), (mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma), math.erf((mu-math_log_if_pos(x2)) / (math.sqrt(2.0)*sigma))
	#part_1 =  from x = d+2*r to x = c_min + d + r
	part_1 = gd(x2, d, mu, sigma, c_min, c_max, r) - gd(x1, d, mu, sigma, c_min, c_max, r)
	# print gd(c_min + d + r - 1, d, mu, sigma, c_min, c_max, r), gd(d + 2*r, d, mu, sigma, c_min, c_max, r)
	# print 'LOL', gd(x2, d, mu, sigma, c_min, c_max, r), gd(x1, d, mu, sigma, c_min, c_max, r)

	#part_2 = from x = c_min + d + r to x = c_max + d + r
	part_2 =  gd(x4, d, mu, sigma, c_min, c_max, r) - gd(x3, d, mu, sigma, c_min, c_max, r)

	# print gd(c_max + d + r, d, mu, sigma, c_min, c_max, r), gd(c_min + d + r, d, mu, sigma, c_min, c_max, r)

	#part_3 = from x = c_max + d + r + 1 to c_min + c_max + d 
	part_3 = gd(x6, d, mu, sigma, c_min, c_max, r) - gd(x5, d, mu, sigma, c_min, c_max, r)

	# print gd(x6, d, mu, sigma, c_min, c_max, r), gd(x5, d, mu, sigma, c_min, c_max, r)
	# print gd(c_min + c_max + d, d, mu, sigma, c_min, c_max, r), gd(c_max + d + r + 1, d, mu, sigma, c_min, c_max, r)
	# print 'from to', c_min + c_max + d, c_max + d + r + 1
	# print "g_d, part_1:", part_1, gd(x2, d, mu, sigma, c_min, c_max, r), gd(x1, d, mu, sigma, c_min, c_max, r)
	# print "g_d, part_2:", part_2, gd(x4, d, mu, sigma, c_min, c_max, r), gd(x3, d, mu, sigma, c_min, c_max, r)
	# print "g_d, part_3:", part_3, gd(x6, d, mu, sigma, c_min, c_max, r), gd(x5, d, mu, sigma, c_min, c_max, r)
	# print "g_d, all parts sum:", part_1 + part_2 + part_3
	# print "part3:", gd(x6, d, mu, sigma, c_min, c_max, r), gd(x5, d, mu, sigma, c_min, c_max, r), x6,x5
	return part_1 + part_2 + part_3

def GapEstimator(mu, sigma, r, observations, c1_len, c2_len=None, method="NR", stepsize=None):
    '''
    Calculates the lower bound (given in http://www.ncbi.nlm.nih.gov/pubmed/22923455). The upper bound is then
    uniquely determined by the lower bound.
    r 		- read length
    Method  - Either "NR" for Newton-Rhapson or linear for going through every gap size (only for testing purposes)
    '''
    # print "LOOOOOOOOL"
    if c2_len is None:
        c2_len = c1_len

    if not stepsize:
    	stepsize = min(math.sqrt( (math.exp(sigma**2) - 1) * math.exp(2*mu + sigma**2))/float(10), 300)
    	# print stepsize
    	# sys.exit()

	mu, sigma, r, c1_len, c2_len = float(mu), float(sigma), float(r), float(c1_len), float(c2_len)
	# Some pruning of extreme observations would be good to have here
	observations.sort()
	# print observations
	min_obs = observations[0]
	# print "MIN obs:", min_obs

	max_obs_mean = sum(observations[-3:])/3.0
	# print max_obs_mean
	c_min = min(c1_len, c2_len)
	c_max = max(c1_len, c2_len)
	mode = math.exp(mu - sigma**2)
	min_gap = max(-int(min_obs/2) + 1, -c_min + int(2*r), mode - observations[-1])
	# print "min gap", min_gap
	# sys.exit()
	# set max to mu + 2*sd
	cutoff_approx_level = math.exp(mu + 0.5*sigma**2) + 5000*math.sqrt( (math.exp(sigma**2) - 1) * math.exp(2*mu + sigma**2) ) 
	# print "cutoff_approx_level", cutoff_approx_level
	upper_quantile = math.exp(mu + 0.5*sigma**2) + 5*math.sqrt( (math.exp(sigma**2) - 1) * math.exp(2*mu + sigma**2) )
	# print upper_quantile
	max_gap = int(upper_quantile - max_obs_mean) if int(upper_quantile - max_obs_mean) > 0 else int(upper_quantile - max_obs_mean)/2
	# max_gap = int(upper_quantile) if upper_quantile < cutoff_approx_level else int(cutoff_approx_level - max_obs)  #emperical_max
	d_lower = min_gap
	d_upper = max_gap
	# print "max gap", max_gap
	# print min_obs, max_obs_mean
	if min_gap > max_gap:
		# print "enter here"
		return min_gap
	# sys.exit()
	# print 'UPPER:', d_upper, cutoff_approx_level, max_gap
	# sys.exit()
    if method == "linear":
    	d_ML = get_d_ML_linear_search(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, stepsize, cutoff_approx_level)
    elif method == "NR":
    	d_ML = get_d_ML_Newton_Raphson(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, cutoff_approx_level)
    	# failed to find good root! do slower but more accurate linear search for this gap
    	if d_ML <= d_lower + 1:
			# print "HEEEREE"
			d_ML = get_d_ML_linear_search(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, stepsize, cutoff_approx_level)

    else:
    	print "Wrong method specified, specify either 'NR' (preffered) or 'linear'. "
    	return None

    return d_ML

def get_d_ML_Newton_Raphson(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, cutoff_approx_level):
	"""
	    Returns the ML-gap as float. The ML gap is searched
	    for with an accuracy of 0.5bp. This is NOT the accuracy of the ML-estimate
	    however since it depends on the number of samples etc.

	    Algorithm:

	    Step 1: Calculate the derivative of the two functions (the g(d) ratio function
	    		and the observation function. This is approximated here
	    		by calculating the function values in two close points (1bp) and using:
	    		f'(x) = (f(x_2) - f(x_1)) / (x_2 - x_1) where x_2 > x_1

	    Step 2: Find the intersection with y, i.e. the intercept, of the two functions.

	    Step 3: The derivatives help us describe the fit of a line in point x. We have
	    		y = kx+m, where k = f'(x) and m=f(x).
	    		We use the two lines to find the intersection point x' between them. The intersection
	    		point x' is used to define the new x.

	    Step 4: Restart algorithm with x'. 
	"""
	n = len(observations)

	x = d_lower
	x_prime = d_upper
	#print 'Current gap:', x #, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',observation_term
	#print 'lower,upper,nr obs:', d_lower, d_upper, n
	#print 'c_min and c_max', c_min, c_max
	#print 'observations:', observations
	import matplotlib.pyplot as plt
	dx = 1

	while True:
		# Step 1
		o_x1 = calc_log_norm(observations, x, mu, sigma)
		o_x2 = calc_log_norm(observations, x+dx, mu, sigma)

		g_d_denominator_x1 = calc_gd(x, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		# print(x, mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, cutoff_approx_level)
		g_prime_d_nominator_x1 = calc_g_prim_d(x, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		# try:
		g_d_ratio_x1 = n*g_prime_d_nominator_x1/g_d_denominator_x1
		# except ZeroDivisionError:
		# 	x += 10

		g_d_denominator_x2 = calc_gd(x+dx, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_prime_d_nominator_x2 = calc_g_prim_d(x+dx, mu, sigma, c_min, c_max, r, cutoff_approx_level)

		# try:
		g_d_ratio_x2 = n*g_prime_d_nominator_x2 /g_d_denominator_x2
		# except ZeroDivisionError:
		# 	# too small starting value
		# 	x += 10

		o_prime = (o_x2 - o_x1) / (x+dx - x)
		g_prime = (g_d_ratio_x2 - g_d_ratio_x1) / (x+dx - x)
		
		# Step 2: y = kx + m ==> m = y - kx
		intercept_o = (o_x1 + o_x2)/2.0 - o_prime*(x + (x + 1))/2.0
		intercept_g = (g_d_ratio_x1 + g_d_ratio_x2)/2.0 - g_prime*(x + (x + 1))/2.0

		#Step 3:
		x_prime = ( intercept_o - intercept_g) / (g_prime - o_prime)
		#print "x:", x, "x_prime:", x_prime, x_prime - x, "derivative","oprime:", o_prime, "g_prime", g_prime 
		# print 'Current gap:', x #, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',observation_term

		###########
		# if observations == [6539, 5449, 4185, 2700, 1294, 2449, 1010, 765, 1362, 2700, 2449, 2946, 5103]: 
		# gap_size, f_1, f_2 = get_fcns(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, 100, cutoff_approx_level)
		# fig, ax = plt.subplots()
		# plt.plot(gap_size, f_1,'-')
		# plt.plot(gap_size, f_2,'-')
		# y_0 = o_prime*-4000 + intercept_o
		# y_1 = o_prime*5000 + intercept_o
		# z_0 = g_prime*-4000 + intercept_g
		# z_1 = g_prime*5000 + intercept_g
		# ax.scatter([-4000, 5000], [y_0, y_1], marker='^', s=150, c='r')
		# ax.plot([-4000, 5000], [y_0, y_1], c='r') 
		# ax.scatter([-4000, 5000], [z_0, z_1], marker='^', s=150, c='b')
		# ax.plot([-4000, 5000], [z_0, z_1], c='b') 
		# plt.ylim((-.05, .05))
		# plt.show()
		#################

		# Step 3
		if x_prime < x:
			# catching bug here: should never be negative slope,
			# this can occur from floating point approximations in any of the 
			# used expressiond when things gets very small or large.
			x_prime = x
			# print 'x_prime moved back'
			# raw_input("press enter")
			break
		elif x_prime  > d_upper:
			# Another heuristic stopping criterion:
			# Derivatives are too sililar so the intersection
			# has high uncertainty
			x_prime = d_upper
			# print 'x_prime sprinted'
			# raw_input("press enter")
			break
		elif math.fabs(x_prime - x) < 0.5:
			break
		else:
			x = x_prime

	# print "FINAL value:", x_prime
	# sys.exit()
	d_ML = x_prime if x_prime < d_upper else d_upper
	#print d_ML
	# a = raw_input("press enter")
	return d_ML

	##################################

def get_d_ML_linear_search(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, stepsize, cutoff_approx_level):
	"""
	    Returns the ML-gap as float. The ML gap is searched
	    for with an accuracy of "stepsize" using a linear search in [d_lower, d_upper]. This is NOT the accuracy of the ML-estimate
	    however since it depends on the number of samples etc.
	"""

	n = len(observations)
	x = []
	y = []
	y_1 = []
	y_2 = []
	g_d_list = []
	g_prim_d_list = []
	d_ML = d_lower
	# this is where the two functions cross and we have our (not nexessarily unique) ML estimate
	# The fist cross is the true ML estimate based on emperical testing
	fcns_has_crossed = False
	# print int(d_lower), int(d_upper), int(stepsize)
	closest_dist_between_fcns = 2**32
	for d in range(int(d_lower), int(d_upper), int(stepsize)):
		#print "GAP", d
		if fcns_has_crossed:
			break
		g_d = calc_gd(d, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_d_list.append(g_d)
		g_prime_d = calc_g_prim_d(d, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_prim_d_list.append(g_prime_d)

		g_d_ratio = n*g_prime_d/g_d
		other_term = calc_log_norm(observations, d, mu, sigma)
		if other_term < g_d_ratio and not fcns_has_crossed:
			# the ML estimate occurs at the first crossing
			d_ML = (d + (d-stepsize))/2.0
			fcns_has_crossed = True
		else:
			if other_term > g_d_ratio and other_term - g_d_ratio < closest_dist_between_fcns:
				closest_dist_between_fcns = other_term - g_d_ratio
				d_ML = (d + (d-stepsize))/2.0

		x.append(d)
		y.append(g_d_ratio - other_term)
		y_1.append(g_d_ratio)
		y_2.append(other_term)
		#print
		#print 'gap:',d, 'g(d)*n', g_d_ratio, "g_d", g_d, "g_prime", g_prime_d, 'other:',other_term


	# import matplotlib.pyplot as plt
	# # plt.plot(x,y,'-')
	# # x1,x2,y1,y2 = plt.axis()
	# # # plt.axis((x1,x2,-0.01,0.001))
	# # plt.show()

	# plt.plot(x,y_1,'-')
	# plt.plot(x,y_2,'-r')
	# x1,x2,y1,y2 = plt.axis()
	# #plt.axis((x1,x2,-0.01,0.001))
	# plt.show()

	# # plt.plot(x,g_prim_d_list,'-')
	# # x1,x2,y1,y2 = plt.axis()
	# # # plt.axis((x1,x2,-0.01,0.001))
	# # plt.show()

	# # plt.plot(x,g_d_list,'-r')
	# # x1,x2,y1,y2 = plt.axis()
	# # # plt.axis((x1,x2,-0.01,0.001))
	# # plt.show()

	# print 'MLGAP', d_ML
	return d_ML



def get_fcns(mu, sigma, r, c_min, observations, c_max, d_lower, d_upper, stepsize, cutoff_approx_level):
	"""
	    Returns the ML-gap as float. The ML gap is searched
	    for with an accuracy of "stepsize" using a linear search in [d_lower, d_upper]. This is NOT the accuracy of the ML-estimate
	    however since it depends on the number of samples etc.
	"""

	n = len(observations)
	x = []
	y = []
	y_1 = []
	y_2 = []
	g_d_list = []
	g_prim_d_list = []
	# this is where the two functions cross and we have our (not nexessarily unique) ML estimate
	# The fist cross is the true ML estimate based on emperical testing
	fcns_has_crossed = False

	for d in range(int(d_lower), int(d_upper), int(stepsize)):
		#print "GAP", d
		g_d = calc_gd(d, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_d_list.append(g_d)
		g_prime_d = calc_g_prim_d(d, mu, sigma, c_min, c_max, r, cutoff_approx_level)
		g_prim_d_list.append(g_prime_d)

		g_d_ratio = n*g_prime_d/g_d
		other_term = calc_log_norm(observations, d, mu, sigma)
		if other_term < g_d_ratio and not fcns_has_crossed:
			# the ML estimate occurs at the first crossing
			d_ML = (d + (d-stepsize))/2.0
			fcns_has_crossed = True

		x.append(d)
		y.append(g_d_ratio - other_term)
		y_1.append(g_d_ratio)
		y_2.append(other_term)


	return x, y_1, y_2











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


