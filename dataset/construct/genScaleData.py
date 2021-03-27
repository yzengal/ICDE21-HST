import os
import sys
from random import randint, sample, shuffle
import numpy as np

class constForScale:
	caseN = 3
	Range = 10**7
	dim = 2
	nList = map(int, [1e5,3e5,5e5,7e5,1e6,3e6,5e6])
	alphas = [2]
	
class CFS(constForScale):
	pass

class baseGenerator:
	def gen(self, n):
		pass	
	
class randomGenerator(baseGenerator):

	def __init__(self, mx):
		self.mx = mx

	def setMx(self, mx):
		self.mx = mx

	def gen(self, n):
		ret = [0] * n
		for i in xrange(n):
			x = randint(-self.mx, self.mx)
			ret[i] = x
		return ret

def genPoints(n):
	gtor = randomGenerator(CFS.Range)
	points = [[0]*CFS.dim for i in xrange(n)]
	for j in xrange(CFS.dim):
		vals = gtor.gen(n)
		vals = map(int, vals)
		for i in xrange(n):
			points[i][j] = vals[i]
	return points
		

def genSynData(n, points, desPath, caseN):
	perms = range(0, n)
	prefix = "sca_%d" % (n)
	tmpFilePath = os.path.join(desPath, prefix)
	if not os.path.exists(tmpFilePath):
		os.mkdir(tmpFilePath)	
	alpha = CFS.alphas[0]
	for i in xrange(0, caseN):
		desFileName = "data_%d.dat" % (i)
		desFileName = os.path.join(tmpFilePath, desFileName)
		with open(desFileName, "w") as fout:
			shuffle(perms)
			beta = randint(10000,alpha*10000)*1.0 / (alpha*10000)
			fout.write("%d %.4f %d\n" % (n, beta, alpha))
			line = " ".join(map(str, perms)) + "\n"
			fout.write(line)
			for j in xrange(n):
				fout.write(" ".join(map(str, points[j]))+"\n")
	
	
def genSynDataSet(desPath, caseN):
	if not os.path.exists(desPath):
		os.mkdir(desPath)
	nmax = max(CFS.nList)
	points = genPoints(nmax)
	for n in CFS.nList:
		genSynData(n, points, desPath, caseN)
	
	
def exp0():
	desPath = "scaleData"
	genSynDataSet(desPath, CFS.caseN)
	
	
if __name__ == "__main__":
	exp0()
	
