#!/usr/bin/python
import sys
from mapper import *
from util import *
from math import *
from numpy import *

class Mapper:
	def __init__(self, file_name, q, k, l, p, ll):
		self.f = Filter(q, k)
		self.data = Data(file_name)
		domain = self.data.domain
		self.f.image(domain)
		domain.bounds()
		print "	range: [%f, %f]" % (domain.min, domain.max)

		self.cover = Cover(l, p)
		self.cover.build(domain.min, domain.max)

		print
		sys.stdout.write("cover:")
		for t in range (0, len(self.cover.coverset)):
			sys.stdout.write("	")
			print "[ %f, %f ]" % (self.cover.coverset[t].a, self.cover.coverset[t].b)
		print

		self.cover.levelsets(domain)
		for t in range(0, len(self.cover.coverset)):
			print "	coverset %d has %d samples" % (t, len(self.cover.coverset[t].samples))
		print

		self.complex = Complex(self.cover)
		self.complex.cluster()
		print "  ------------------------------------------------------ "
		for c in self.complex .cover.coverset:
			samples = c.samples
			clusters = c.clustering.clusters
			print
			print "coverset [%f, %f] has %d samples and %d clusters" % (c.a, c.b, len(samples), len(clusters))
			for i in range(0,len(clusters)):
				print "	cluster %d has %d samples" % (i, len(clusters[i].samples))
		print
		print "  ------------------------------------------------------ "
		print

		print "%d samples total" % len(domain.samples)
		print

		for s in domain.samples:
			if len(s.clusters) > 1:
				for c in s.clusters:
					for t in s.clusters:
						if c != t:
							self.complex.newEdge(c.vertex, t.vertex)

		self.complex.spring_embedding()


# --------------------------------------- #

print # data import
file_name = "data/TCGA-BRCA-L3-S35.txt"
# file_name = "data/test.txt"
# print "importing data file " + file_name +"..."
# data = Data(file_name)

q = 2		# filter function inner-power, root denominator
k = 1		# filter function root numerator
l = 0.5		# cover-set interval length
p = 1/7.	# cover-set interval percent overlap

# o = mapper(data, q, k, l, p) # [ [ filter_values, levelset_edges ], [ cover, levelsets ] ]

# sample = o[0] # filter values and levelset edges of the domain
# image = o[1] # the covering of the filter function over the data matrix, and the levelsets of the covering
# F = sample[0] 	# [ f(V1), ... , f(Vi), ... , f(Vn) ]
# C = image[0] 	# [ [a, a + lp], ... , [a + il(1 - p), a + l(i(1 - p) + i], ... , [a + sl(1 - p), a + l(s(1 - p) + 1]]
# L = image[1]	# [ { i | Vi in f-1(C0) }, ... , { i | Vi in f-1(Ct) }, ... , { i | Vi in f-1(Cs) } ]
# T = sample[1] 	# [ { t | f(V1) in Ct }, ..., { t | f(Vi) in Ct }, ... , { t | f(Vn) in Ct } ]

# ---------------------------------------------- #

mapper = Mapper(file_name, q, k, l, p, l*p)