#!/usr/bin/python
import sys
from mapper import *
from util import *
from math import *
from numpy import *

class Mapper:
	# def __init__(self, file_name, q, k, l, p, ll):
	def __init__(self, data, q, k, l, p, ll, file_out):
		self.f = Filter(q, k)
		self.data = data
		# self.data = Data(file_name)
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

		self.complex.spring_embedding(file_out)
		print domain.min
		print domain.max
# --------------------------------------- #

def load_data(file_name):
	file = open(file_name)
	data = []
	with file as f:
		i = 0
		for line in f:
			data.append([])
			string = line.split()
			j = 0
			for x in string:
				data[i].append(float(x))
			i = i + 1
	return data

def normalize_healthy(normal, tumor):
	normal_mat = matrix(normal)
	tumor_mat = matrix(tumor)
	normal_T = normal_mat.T.tolist()
	tumor_mat_T = tumor_mat.T
	i = 0
	for col in normal_T:
		summ = 0
		for x in col:
			summ = summ + x
		tumor_mat_T[i] = tumor_mat_T[i] - (summ/len(col))
		i = i + 1
	return tumor_mat_T.T.tolist()

def scale_col(data, scale):
	data_mat_T = matrix(data).T
	i = 0
	for col in data_mat_T.tolist():
		maxx = float("-inf")
		for x in col:
			if x > maxx:
				maxx = x
		data_mat_T[i] = ((data_mat_T[i]/maxx)*2 - 1) * scale
		i = i + 1
	return data_mat_T.T.tolist()

# --------------------------------------- #

print # data import
file_name = "data_new/protein/tumor/all_protein_tumor.txt"
tumor_data = "data_new/protein/all_tumor_protein.txt"
normal_data = "data_new/protein/all_normal_protein.txt"
# file_name = "data_new/protein/normal/all_protein_normal.txt"
# file_name = "data/TCGA-BRCA-L3-S35.txt"
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

# mapper = Mapper(file_name, q, k, l, p, l*p)

tumor_norm_output = "data_new/protein/tumor_protein_normalized.txt"

tumor_mat = load_data(tumor_data)
normal_mat = load_data(normal_data)
tumor_mat_norm = normalize_healthy(normal_mat, tumor_mat)
# tumor_mat_norm_1 = scale_col(tumor_mat_norm, 1)
savetxt(tumor_norm_output, tumor_mat_norm, delimiter='\t')

data = Data(None)
data.n = len(tumor_mat_norm)
data.m = len(tumor_mat_norm[0])
data.matrix = []
data.data = tumor_mat_norm
data.rows = []
data.cols = []
data.samples = []
data.domain = Domain(data)

q = input('q: ')
k = input('k: ')
p = input('p: ')
l_min = input('l_min: ')
l_max = input('l_max: ')
l_d = input('l_d: ')
count = 0;

while l_min + l_d*count <= l_max:
	file_out = "images/image"
	if count < 10000:
		file_out = file_out + "0"
		if count < 1000:
			file_out = file_out + "0"
			if count < 100:
				file_out = file_out + "0"
				if count < 10:
					file_out = file_out + "0"
	file_out = file_out + str(count+1) + ".png"
	mapper = Mapper(data, q, k, l, p, l*p, file_out)
	count = count + 1
