'''
    Generate data from gapest model with normal or log-normal fragment distribution. 
    '''

import argparse
import sys
import os
import math
import random
import numpy as np
from math import exp, pi, sqrt
sys.path.insert(0,'/Users/kxs624/Documents/workspace/mathstats/')
from mathstats.normaldist import normal
from scipy.stats import norm,lognorm

from mathstats import log_normal_param_est as lnpe
from mathstats.normaldist.truncatedskewed import param_est as pe
import matplotlib.pyplot as plt

def generate(n, mu, sigma, gap, c_min, c_max, distribution="lognormal"):
	r=100 #readlength
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
	# print 'Mean all observations:', mean_samples
	std_dev_samples = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_samples + mean_samples ** 2), samples))) / (len(samples) - 1)) ** 0.5
	# print 'STDDEV all samples:', std_dev_samples

	#observations_over_gap = [ int(round(max(s-gap,0),0)) for s in samples]
	#observations_over_gap = filter(lambda x: x>0, observations_over_gap)
	#print sum(observations_over_gap)/len(observations_over_gap)

	samples_kept = []
	for s in samples:
		if s > c_min + c_max + gap or s < gap or s < -gap or c_min <= -gap:
			continue
		p = random.uniform(0,1)
		# print 'lol',max_sample, gap, (s-gap-2*r) / max(0,(max_sample-gap-2*r))
		if p < (s-gap-2*r) / max(0,(max_sample-gap-2*r)):
			samples_kept.append(s)

	# print len(samples_kept)
	if len(samples_kept) <= 1:
		return []
	# print "gap:", gap

	mean_samples =  sum(samples_kept)/len(samples_kept)
	# print 'Mean conditional fragment size:', mean_samples
	std_dev_samples = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_samples + mean_samples ** 2), samples_kept))) / (len(samples_kept) - 1)) ** 0.5
	# print 'STDDEV conditional fragment size:', std_dev_samples


	observations_over_gap = [ int(round(max(s-gap,0),0)) if gap > 0 else int(round(max(s - gap,0),0)) for s in samples_kept]
	observations_kept = filter(lambda x: x>0, observations_over_gap)
	mean_obs =  sum(observations_kept)/len(observations_kept)
	# print 'Mean conditional observed size:', mean_obs
	std_dev_obs = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_obs + mean_obs ** 2), observations_kept))) / (len(observations_kept) - 1)) ** 0.5
	# print 'STDDEV conditional observed size:', std_dev_obs
	# print
	# print
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

def plot_estimate_uncertainty(n, mu, sigma, c_min, c_max, distribution="lognormal"):
	if distribution == 'lognormal':
		gaps = []
		estimates_log = []
		estimates_normal = []
		for gap in range(max(-2000,-c_min+200),7001, 500):
			# print "GAP:", gap
			for i in range(50):
				sample = generate(n, mu, sigma, gap, c_min, c_max)
				if not sample:
					continue
				gap_est = lnpe.GapEstimator(mu, sigma, 100, sample, c_min, c2_len=c_max)

				mu_norm = math.exp(mu + sigma**2/2.0)
				sigma_norm = math.sqrt( (math.exp(sigma**2) - 1)* math.exp(2*mu+sigma**2))
				mean_sample = sum(sample)/float(len(sample))
				gap_est_normal = pe.GapEstimator(mu_norm, sigma_norm, 100, mean_sample, c_min, c2_len=c_max)

				gaps.append(gap)
				estimates_log.append(gap_est)
				estimates_normal.append(gap_est_normal)

		plt.plot(gaps, estimates_log,'o', alpha=.5)
		plt.plot(gaps, estimates_normal,'^', alpha=.5)
		plt.plot(gaps, gaps,'-')
		plt.show()
	else:
		print("Specify lognormal or do not set this argument.")
		return None

def main(args):
	# observations_kept = generate(args.n, args.mu, args.sigma, args.gap, args.c_min, args.c_max, args.distribution)
	plot_estimate_uncertainty(args.n, args.mu, args.sigma, args.c_min, args.c_max, args.distribution)
	# plot_conditional_means(args.n, args.mu, args.sigma, args.min_gap, args.max_gap, args.distr)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate data from gapest model with normal or log-normal fragment distribution.")
	#n, mu, sigma, gap, c_min, c_max,
	parser.add_argument('n', type=int, help='Number of samples')
	parser.add_argument('mu', type=float, help='mean')
	parser.add_argument('sigma', type=float, help='std deviation')
	parser.add_argument('gap', type=float, help='gap')
	parser.add_argument('c_min', type=int, help='c_min')
	parser.add_argument('c_max', type=int, help='c_max')
	parser.add_argument('--distribution', type=str, default="lognormal", help='Either normal or lognormal')

	# parser.add_argument('transcripts', type=str, nargs="+", help='Path to the transcript fasta file')
	# parser.add_argument('--blast_xml', dest="blast_xml", type=str, nargs="+", help='Path to the blast records')
	# parser.add_argument('--pickled_parse', dest="pickled_parse", action="store_true", help='Path to the blast records')

	parser.add_argument('outfolder', type=str, help='Output path of results')

	args = parser.parse_args()

	outfolder = args.outfolder
	if not os.path.exists(outfolder):
	    os.makedirs(outfolder)
	main(args)
