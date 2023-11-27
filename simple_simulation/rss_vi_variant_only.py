import numpy as np 
import os
import sys
import pdb








class RSS_VI_VARIANT_ONLY(object):
	def __init__(self, LD=None, z=None, max_iter=100):
		self.LD = LD
		self.z = z
		self.max_iter = max_iter
	def fit(self):
		# Quick error check
		if self.LD is None or self.z is None:
			raise ValueError('RSS_VI_VARIANT_ONLY requires LD and z-scores')

