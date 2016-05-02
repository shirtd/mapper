#!/usr/bin/python
import sys
from mapper import *
from math import *
from numpy import *
import igraph

# ------------- filter function ------------- #

def filter_fun(x, q, k):
	f = 0;
	for y in x.v:
		f = f + pow(abs(float(y)), float(q))
	return f ** (float(k)/float(q))


# ------------ utility functions ------------ #

def bounds(samples):
	minimum = float("inf")
	maximum = float("-inf")
	for x in samples:
		if x.f < minimum:
			minimum = x.f
		if x.f > maximum:
			maximum = x.f
	return [ minimum, maximum ]
def norm(u):
	prod = 0
	for i in range(0, len(u)):
		prod = prod + (u[i] ** 2.)
	prod = prod ** (1./2)
	return prod
def dot(u, v):
	prod = 0
	for i in range(0, len(u)):
		prod = prod + u[i]*v[i]
	return prod
def correlation_dist(x, y):
	u = (matrix(x.v) - x.mean).tolist()[0]
	v = (matrix(y.v) - y.mean).tolist()[0]
	return (1 - (dot(u,v)/(norm(u)*norm(v))))

  
# ------------- data structures ------------- #

class Data:
	def __init__(self, file_name):
		if file_name != None:
			self.file_name = file_name
			self.file = open(file_name)
			self.n = None
			self.m = None
			self.matrix = []
			self.data = []
			self.rows = []
			self.cols = []
			self.samples = []
			self.domain = None
			self.process(self.file)
	def process(self, file):
		i = 0
		with file as f:
			for line in f:
				self.matrix.append([])
				self.data.append([])
				string = line.split()
				j = 0
				for x in string:
					if (i == 0) & (j != 0):
						self.cols.append(x)
						self.matrix.append(x)
					elif (i != 0) & (j == 0):
						self.rows.append(x)
						self.matrix.append(x)
					elif (i != 0) & (j != 0 ):
						self.data[i].append(float(x))
						self.matrix.append(float(x))
					else:
						self.matrix.append(x)
					j = j + 1
				i = i + 1
		del self.data[0]
		self.m = len(self.cols)
		self.n = len(self.rows)
		self.domain = Domain(self)

class Domain:
	def __init__(self, data):
		self.samples = []
		print len(data.data)
		print len(data.rows)
		i = 0
		for v in data.data:
			s = Sample(v, data.rows[i])
			self.samples.append(s)
			i = i + 1
		self.min = None
		self.max = None
	def bounds(self):
		bdy = bounds(self.samples)
		self.min = bdy[0]
		self.max = bdy[1]
class Sample:
	def __init__(self, v, label):
		self.v = v
		self.label = label
		self.levelsets = []
		self.clusters = []
		self.f = None
		self.norm = norm(self.v)
		self.mean()
	def f(self, f):
		self.f = f.apply(self.v)
	def mean(self):
		avg = 0
		for x in self.v:
			avg = avg + x
		self.mean = avg/len(self.v)
	def addCluster(self, c):
		for d in self.clusters:
			if d == c:
				return False
		self.clusters.append(c)
		return True
	def setCluster(self, c):
		levelset = c.clustering.levelset
		i = 0
		for l in self.levelsets:
			if l == levelset:
				self.clusters[i] = c
			i = i + 1


# ---------------- geometry ---------------- #

class Vertex:
	def __init__(self, cluster, i):
		self.x = None
		self.y = None
		self.cluster = cluster
		self.i = i
	def d(self,vertex):
		mn = float("inf")
		for s in self.cluster.samples:
			for t in vertex.cluster.samples:
				d = correlation_dist(s, t)
				if d < mn:
					mn = d
class Edge:
	def __init__(self, u, v):
		self.u = u
		self.v = v
class Complex:
	def __init__(self, cover):
		self.cover = cover
		self.l = self.cover.n
		self.min_size = 10
		self.max_size = 50
		self.min_samples = float("inf")
		self.max_samples = float("-inf")
		self.min_f = self.cover.min
		self.max_f = self.cover.max
		self.pal = igraph.RainbowPalette(self.l)
		self.vertices = []
		self.edges = []
		self.edge_mat = []
		self.weights = []
		self.dim = 2
		self.n = 0
		self.m = 0
		self.graph = igraph.Graph()
		self.width = 500
		self.height = 500
	def cluster(self):
		cover_samples = 0
		for s in self.cover.coverset:
			self.addClusters(s.cluster())
			cover_samples = cover_samples + len(s.samples)
		for v in self.vertices:
			size = self.size(v.cluster.n)
			self.graph.vs[v.i]['size'] = size
			# print "	color of vertex %d: " % v.i
			# sys.stdout.write("	")
			# print self.graph.vs[v.i]['color']
	def color(self, cluster):
		levelset = cluster.clustering.levelset
		color = self.pal.get(levelset.i)
		# print "	color of vertex %d: " % cluster.vertex.i
		# sys.stdout.write("	")
		# print color
		return color
	def size(self, n):
		size = ((float(n) - self.min_samples)/self.max_samples)*(self.max_size - self.min_size) + self.min_size
		# print "size of vertex %d: %f pixels" % (n, size) 
		return size
	def addClusters(self, clustering):
		for c in clustering.clusters:
			self.newVertex(c)
	def newVertex(self, cluster):
		# print " adding vertex %d" % self.n
		v = Vertex(cluster, self.n)
		cluster.vertex = v
		if cluster.n < self.min_samples:
			self.min_samples = cluster.n
		if cluster.n > self.max_samples:
			self.max_samples = cluster.n
		self.vertices.append(v)
		self.graph.add_vertex()		
		self.graph.vs[self.n]['color'] = self.color(cluster)
		self.graph.vs[self.n]['label'] = "class"+str(self.n) # +"\nlevelset = ["+str(cluster.clustering.levelset.a)+", "+str(cluster.clustering.levelset.a)+"]"
		self.n = self.n + 1
		new_row = [None]
		for r in self.edge_mat:
			new_row.append(None)
			r.append(None)
		self.edge_mat.append(new_row)
		return v
	def newEdge(self, u, v):
		if u != v:
			e = self.edge_mat[u.i][v.i]
			if e == None:
				# print " linking vertices %d and %d" % (u.i, v.i)
				e = Edge(u,v)
				self.edges.append(e)
				self.edge_mat[u.i][v.i] = e
				self.edge_mat[v.i][u.i] = e
				self.graph.add_edge(u.i, v.i)
				self.weights.append(u.d(v))
				self.m = self.m + 1
			return e
		else: return None
	def spring_embedding(self, file_name):
		minx = 0
		maxx = self.width
		miny = 0
		maxy = self.height
		area = self.width*self.height	# default is the number of vertices
		maxiter = 500	# number of iterations to perform (default 500)
		maxdelta = self.n	# maximum distance to move a vertex in an iteration (default is the number of vertices)
		coolexp = 1.5	# cooling component of the simulated annealing (default 1.5)
		repulserad = self.n**3.	# radius at which vertex-vertex repulsion cancels out attraction of adjacent vertices (default len(vertices)^3)
		#seed = None 	# if None, uses a random starting layout for the algorithm. If a matrix (list of lists) uses the given matrix as the starting position
		# print self.weights
		# self.graph.layout_fruchterman_reingold(self.weights, maxiter, maxdelta, area, coolexp, repulserad,minx, maxx, miny, maxy, 0, 0, seed, self.dim)
		# self.complex.embed(self.graph)
		s = []
		height = 600
		width = 600
		d_x = (width - 100)/self.l
		pos_x = 50		
		pos_y = 50
		for l in self.cover.coverset:
			if (len(l.clustering.clusters) == 0):
				pos_y = height/2
				s.append([pos_x, height/2])
			else:
				d_y = (height - 100)/len(l.clustering.clusters)
				for c in l.clustering.clusters:
					s.append([pos_x, pos_y])
					pos_y = pos_y + d_y
			pos_x = pos_x + d_x
			pos_y = 50
		layout = self.graph.layout_fruchterman_reingold(seed=s)
		igraph.plot(self.graph, file_name, layout = layout)
		print self.pal.get(0)
		print self.pal.get(self.l-1)
		# self.graph.__plot__(context=None, bbox=(0,0,600,600), palette=self.pal, layout=layout, target=file_name, vertex_order_by='asc')

# class Drawing:
# 	def __init__(self, width, height, K):
# 		self.width = width	# can be a sequence or iterable, or even an edge attribute name
# 		self.height = height
# 		self.complex = K
# 		self.n = len(self.complex.vertices)
# 		self.graph = igraph.Graph()
# 		self.weights = self.complex.weights
# 		self.dim = self.complex.dim
# 	def addVertex(self, vertex):
# 		self.graph.add_vertex
# 		g.vs[0]['color'] = cl_red
# 	def spring_embedding(self):
# 		self.graph.add_vertices(len(self.complex.vertices))
# 		for e in self.complex.edges:
# 			print "	adding edge from vertex %d to %d" % (e.u.i, e.v.i)
# 			self.graph.add_edge(e.u.i, e.v.i)
# 		minx = 0
# 		maxx = self.width
# 		miny = 0
# 		maxy = self.height
# 		area = self.width*self.height	# default is the number of vertices
# 		maxiter = 500	# number of iterations to perform (default 500)
# 		maxdelta = self.n	# maximum distance to move a vertex in an iteration (default is the number of vertices)
# 		coolexp = 1.5	# cooling component of the simulated annealing (default 1.5)
# 		repulserad = self.n**3.	# radius at which vertex-vertex repulsion cancels out attraction of adjacent vertices (default len(vertices)^3)
# 		seed = None 	# if None, uses a random starting layout for the algorithm. If a matrix (list of lists) uses the given matrix as the starting position
# 		# print self.weights
# 		# self.graph.layout_fruchterman_reingold(self.weights, maxiter, maxdelta, area, coolexp, repulserad,minx, maxx, miny, maxy, 0, 0, seed, self.dim)
# 		# self.complex.embed(self.graph)
# 		layout = self.graph.layout_fruchterman_reingold()
# 		igraph.plot(self.graph, layout = layout)
