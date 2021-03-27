import os
import sys
from random import randint, sample, shuffle
import numpy as np

class constForSyn():
	caseN = 5
	Range = 10**7
	tot_num = 10**5
	dim = 10
	muList = [int(p*Range) for p in [0.1, 0.3, 0.5, 0.7, 0.9]]
	sigmaList = [int(p*Range) for p in [0.10, 0.15, 0.20, 0.25, 0.3]]
	pointFiles = ["uni", "nor", "exp"]
	alphas = [2]
	
class CFR(constForSyn):
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

class normalGenerator(baseGenerator):

	def __init__(self, mu, sigma):
		self.mu = mu
		self.sigma = sigma

	def gen(self, n, lb = None, rb = None):
		# print lb, rb
		ret = np.random.normal(self.mu, self.sigma, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		# print ret
		return ret

	def setMu(self, mu):
		self.mu = mu

	def setSigma(self, sigma):
		self.sigma = sigma


class uniformGenerator(baseGenerator):

	def __init__(self, low, high):
		self.low = low
		self.high = high

	def gen(self, n, lb = None, rb = None):
		ret = np.random.uniform(self.low, self.high, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setLow(self, low):
		self.low = low

	def setHigh(self, high):
		self.high = high

class expGenerator(baseGenerator):

	def __init__(self, mu):
		self.mu = mu

	def gen(self, n, lb = None, rb = None):
		ret = np.random.exponential(self.mu, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setMu(self, mu):
		self.mu = mu

def genPoints(gtor, n):
	points = [[0]*CFR.dim for i in xrange(n)]
	for j in xrange(CFR.dim):
		vals = gtor.gen(n)
                vals = map(int, vals)
		for i in xrange(n):
			points[i][j] = vals[i]
	return points
		

def genSynData(points, desPath, prefix, caseN):
	n = len(points)
	perms = range(0,n)
	for k in CFR.alphas:
		tmpFilePath = "%s_%d" % (prefix, k)
		tmpFilePath = os.path.join(desPath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)		
		for i in xrange(2,caseN):
			desFileName = "data_%d.dat" % (i)
			desFileName = os.path.join(tmpFilePath, desFileName)
			with open(desFileName, "w") as fout:
				shuffle(perms)
				beta = randint(10000,k*10000)*1.0 / (k*10000)
				fout.write("%d %.4f %d\n" % (n, beta, k))
				line = " ".join(map(str, perms)) + "\n"
				fout.write(line)
				for j in xrange(n):
					fout.write(" ".join(map(str, points[j]))+"\n")
	
def genSynDataSet(desPath, caseN):
	if not os.path.exists(desPath):
		os.mkdir(desPath)
	for prefix in CFR.pointFiles:
		if prefix=="uni":
			for mu in CFR.muList:
				gtor = randomGenerator(mu)
				points = genPoints(gtor, CFR.tot_num)
				genSynData(points, desPath, prefix+"_"+str(mu), caseN)
		elif prefix=="nor":
			mu, sig = CFR.muList[2], CFR.sigmaList[2]
			for mu in CFR.muList:
				gtor = normalGenerator(mu, sig)
				points = genPoints(gtor, CFR.tot_num)
				genSynData(points, desPath, prefix+"_mu"+"_"+str(mu), caseN)			
			mu, sig = CFR.muList[2], CFR.sigmaList[2]
			for sig in CFR.sigmaList:
				gtor = normalGenerator(mu, sig)
				points = genPoints(gtor, CFR.tot_num)
				genSynData(points, desPath, prefix+"_sig"+"_"+str(sig), caseN)			
		elif prefix=="exp":
			for mu in CFR.muList:
				gtor = expGenerator(mu)
				points = genPoints(gtor, CFR.tot_num)
				genSynData(points, desPath, prefix+"_"+str(mu), caseN)		
	
	
def exp0():
	desPath = "synData"
	genSynDataSet(desPath, CFR.caseN)
	
	
if __name__ == "__main__":
	exp0()
	
