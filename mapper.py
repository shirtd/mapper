#!/usr/bin/python
import sys
from util import *
from math import *
from numpy import *

class Filter:
	def __init__(self, q, k):
		self.q = float(q)
		self.k = float(k)
	def apply(self, x):
		x.f = filter_fun(x, self.q, self.k)
	def image(self, domain):
		for x in domain.samples:
			self.apply(x)

class Cover:
	def __init__(self, l, p):
		self.l = float(l)
		self.p = float(p)
		self.n = 0
		self.min = None
		self.max = None
		self.coverset = []
	def build(self, a, b):
		self.min = a
		print "		l = %f, p = %f" % (self.l, self.p)
		r = ((b-a)/self.l - 1)/(1-self.p)
		s = int(ceil(r))
		print "		r = %f => s = %d coversets" % (r, s)
		self.min = a
		self.max = a + (s + 1)*self.l - s*self.l*self.p
		print "	adjusted range: [%f, %f]" % (self.min, self.max)
		for i in range (0, int(s)+1):
			sys.stdout.write("	")
			aa = self.min + i*self.l - i*self.l*self.p
			bb = self.min + (i+1)*self.l - i*self.l*self.p
			self.coverset.append(Levelset(aa, bb))
			self.coverset[self.n].i = self.n
			print "added levelset %d" % self.coverset[self.n].i
			self.n = self.n + 1

	def levelsets(self, domain):
		print "computing levelsets..."
		for x in domain.samples:
			for c in self.coverset:
				c.add(x)
			x.clusters = [ None for c in range(len(x.levelsets)) ];
	def distance_mat(self):
		for c in self.coverset:
			c.distance_mat()
	def cluster(self, ll):
		for c in self.cluster:
			c.cluster(ll)

class Levelset:
	def __init__(self, a, b):
		self.a = a
		self.b = b
		self.i = None
		self.mid = (a + b)/2
		self.samples = []
		self.clustering = None
		self.D = []
	def add(self, x):
		if (self.a <= x.f) & (x.f <= self.b):
			self.samples.append(x)
			x.levelsets.append(self)
			return True
		else:
			return False
	def distance_mat(self):
		nn = len(self.samples)
		if (nn != 0):
			self.D = [[0. for x in range(nn)] for x in range(nn)] 
			for i in range(0,nn):
				for j in range(i+1,nn):
					if i == j:
						self.D[i][j] = 0.
					else:
						self.D[i][j] = correlation_dist(self.samples[i], self.samples[j])
						self.D[j][i] = self.D[i][j]
	def cluster(self):
		print
		print "CLUSTERING LEVELSET WITH %d SAMPLES" % len(self.samples)		
		# if len(self.samples) > 1:
		print "computing distance matrix..."
		self.distance_mat()
		print matrix(self.D)
		print
		print "clustering... "
		self.clustering = Clustering(self)
		return self.clustering

class Clustering:
	def __init__(self, levelset):
		self.L = 0
		self.m = 0
		self.levelset = levelset
		self.D = self.levelset.D
		self.clusters = []
		self.best = None
		for i in range(0, len(self.levelset.samples)):
			self.clusters.append(Cluster(self, i))
		print "	%d clusters" % len(self.clusters)
		while len(self.clusters) > 2:
			self.process()
			print "updating matrix..."
			if self.update():
				print matrix(self.D)
				print
			else: break
		print
	def process(self):
		best_d = float("inf")
		for i in range(0, len(self.clusters)):
			cur = self.clusters[i]
			cur.getClosest(self)
			if cur.min_d < best_d:
				self.best = cur
				best_d = cur.min_d
			# 	print "	       > NEW BEST: clusters %d and %d with %f" % (cur.i, cur.closest.i, best_d)
			# print "		best: cluster[%d].min_d = %f" % (self.best.i, self.best.min_d)
			# print "			cluster[%d].min_d = %f" % (i, cur.min_d) 
		if self.best.i != min([ self.best.i, self.best.closest.i]):
			self.best = self.best.closest
		i = self.best.i
		j = self.best.closest.i
		# print "	best pair: clusters %d and %d a distance %f (%f) apart" % (i, j, self.best.min_d, self.best.closest.min_d)
		# print "		sanity check... D[%d][%d] = %f, D[%d][%d] = %f" % (i, j, self.D[i][j], j, i, self.D[j][i])
		# print "			cluster %d is closest to cluster %d; cluster %d is closest to cluster %d" % (self.best.i, self.best.closest.i, self.best.closest.i, self.best.closest.closest.i)
		# print
	def update(self):
		min_i = self.best.i
		max_i = self.best.closest.i
		self.best.merge()
		self.best = None
		print "	removing cluster %d" % max_i		
		del self.clusters[max_i]
		nn = len(self.clusters)
		for k in range(0, nn):
			print "	cluster %d -> %d" % (self.clusters[k].i, k)
			self.clusters[k].i = k
		if nn > 2:
			D = [[0. for x in range(nn)] for x in range(nn)]
			j_mod = 0
			for i in range(0, nn):
				if i >= max_i: i_mod = 1
				else: i_mod = 0
				for j in range(i, nn):
					if i == j:
						D[i][j] = 0.
					else:
						if j >= max_i: j_mod = 1
						else: j_mod = 0
						if (i != min_i) & (j != min_i):
							D[i][j] = self.D[i+i_mod][j+j_mod]
							D[j][i] = D[i][j]
						elif i == min_i:
							D[i][j] = min([self.D[min_i][j+j_mod], self.D[max_i][j+j_mod]])
							D[j][i] = D[i][j]
						elif j == min_i:
							D[i][j] = min([self.D[i+i_mod][min_i], self.D[i+i_mod][min_i]])
							D[j][i] = D[i][j]
			self.D = D
			return True
		else: return False

class Cluster:
	def __init__(self, clustering, i):
		self.clustering = clustering
		self.i = i
		self.n = 0
		self.samples = []
		self.addSample(self.clustering.levelset.samples[self.i])
		self.min_d = float("inf")
		self.closest = None
		self.vertex = None
	def addSample(self, s):
		self.samples.append(s)
		s.setCluster(self)
		self.n = self.n + 1
	def getClosest(self, clustering):
		self.min_d = float("inf")
		self.clustering = clustering
		for j in range(0, len(self.clustering.clusters)):
			if (self.clustering.D[self.i][j] < self.min_d) & ( self.i != j ):
				self.min_d = self.clustering.D[self.i][j]
				self.closest = self.clustering.clusters[j]
	def merge(self):
		print "	merging clusters %d and %d to position %d" % (self.i, self.closest.i, self.i)
		for s in self.closest.samples:
			self.addSample(s)
		self.min_d = float("inf")
		self.closest = None
		return self
