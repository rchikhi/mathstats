'''
    Generate data from gapest model with normal or log-normal fragment distribution. 
    '''

import sys
import math
from math import exp, pi, sqrt
from mathstats.normaldist import normal
from scipy.stats import norm,lognorm
def generate(n, mu, sigma, gap, distribution="normal"):

	if
	samples = norm.rvg()

	observations_over_gap = [max(s-gap,0) for sample in samples]
	observations_over_gap = filter(lambda x: x>0, samples_over_gap)

	min_observation = min(observations_over_gap)
	max_observation = max(observations_over_gap)

	samples_kept = []
	for o in observations_over_gap:
		p = math.randv(0,1)
		if p < o / max_observation:
			samples_kept.append(o)

	return samples_kept



