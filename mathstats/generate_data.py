'''
    Generate data from gapest model with normal or log-normal fragment distribution. 
    '''

import sys
import math
import random
import numpy as np
from math import exp, pi, sqrt
from mathstats.normaldist import normal
from scipy.stats import norm,lognorm

def generate(n, mu, sigma, gap, distribution="normal"):

	if distribution == 'normal':
		samples = norm.rvs(loc=mu, scale=sigma, size=2*n)

	elif distribution == 'lognormal':
		logsample = norm.rvs(loc=mu, scale=sigma, size=max(1,int(gap)/100)*n)
		samples = np.exp(logsample)
	else:
		print("Specify normal, lognormal or do not set this argument.")
		return None

	min_sample = min(samples)
	max_sample = float(max(samples))


	mean_samples =  sum(samples)/len(samples)
	print 'Mean all observations:', mean_samples
	std_dev_samples = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_samples + mean_samples ** 2), samples))) / (len(samples) - 1)) ** 0.5
	print 'STDDEV all samples:', std_dev_samples

	#observations_over_gap = [ int(round(max(s-gap,0),0)) for s in samples]
	#observations_over_gap = filter(lambda x: x>0, observations_over_gap)
	#print sum(observations_over_gap)/len(observations_over_gap)

	samples_kept = []
	for s in samples:
		p = random.uniform(0,1)
		if p < (s-gap) / max(0,(max_sample-gap)):
			samples_kept.append(s)

	print len(samples_kept)
	print "gap:", gap

	mean_samples =  sum(samples_kept)/len(samples_kept)
	print 'Mean conditional fragment size:', mean_samples
	std_dev_samples = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_samples + mean_samples ** 2), samples_kept))) / (len(samples_kept) - 1)) ** 0.5
	print 'STDDEV conditional fragment size:', std_dev_samples

	#print 'Mean conditional fragment length:', sum(samples_kept)/len(samples_kept)

	observations_over_gap = [ int(round(max(s-gap,0),0)) for s in samples_kept]
	observations_kept = filter(lambda x: x>0, observations_over_gap)
	mean_obs =  sum(observations_kept)/len(observations_kept)
	print 'Mean conditional observed size:', mean_obs
	std_dev_obs = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_obs + mean_obs ** 2), observations_kept))) / (len(observations_kept) - 1)) ** 0.5
	print 'STDDEV conditional observed size:', std_dev_obs
	print
	print
	return observations_kept


def plot_conditional_means(n, mu, sigma, min_gap, max_gap, distr="normal"):
	import matplotlib.pyplot as plt
	conditional_means = []
	x = [] 
	sample_size = []
	for gap in range(min_gap, max_gap, 500):
		s = generate(n, mu, sigma, gap, distribution=distr)
		conditional_means.append(sum(s)/len(s))
		#sample_size.append(len(s))
		x.append(gap)

	plt.plot(x, conditional_means)
	#plt.plot(x, sample_size)
	plt.show()


