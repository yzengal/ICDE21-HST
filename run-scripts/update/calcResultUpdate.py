#!/usr/bin/env python
import os
import sys
import datetime
from genUpdateData import constForReal
from copy import deepcopy

class constForCalcResult(constForReal):
	caseN = 4
	execNames = ["reHSTO", "HSF"]
	infos = ["distortion", "runtime", "memusage"]
	infoN = len(infos)
	begBid = 1
	endBid = 5 * infoN
	
class CFCR(constForCalcResult):
	pass


def parseResult(resFileName, line):
	tmpList = line.strip().split(" ")
	if len(tmpList) < 5:
		print "result is wrong: %s: %s" % (resFileName, line)
		return [0.0] * CFCR.infoN
	maxDistor, avgDistor, runtime, memused = tmpList[1:]
	maxDistor, avgDistor, runtime, memused = float(maxDistor), float(avgDistor), float(runtime), float(memused)/(1024.0*1024)
	return maxDistor, runtime, memused


def fetchResult(fileName):
	ret = []
	with open(fileName, "r") as fin:
		lines = fin.readlines()
		for line in lines:
			if not line:
				continue
			tmpList = parseResult(fileName, line)
			ret.append(tmpList)
	return ret[1:-1]


def genBlock(execNames, infoAll, cntAll, bid):
	lines = []
	
	for iid,infoName in enumerate(CFCR.infos):
		tmpLines = []
		tmpLines.append( "%%%%%%%% %s\n" % (infoName) )
		for eid,execName in enumerate(execNames):
			tmpList = []
			for uid in xrange(CFCR.insertN):
				tot = infoAll[eid][uid][iid]
				cnt = cntAll[eid]
				if cnt == 0:
					avg = 0.0
				else:
					avg = tot*1.0/cnt
				tmpList.append(avg)
			tmpLine = ", ".join(map(lambda s:"%.2f" % (s), tmpList))
			line = "%6s%d = [%s];\n" % (execName, bid, tmpLine)
			tmpLines.append(line)
		lines += tmpLines
		bid += 1
	return lines


def calcResultByVariable(srcFilePath, resFilePath, bid):
	execNames = CFCR.execNames
	caseN = CFCR.caseN
	infoN = CFCR.infoN
	lines = []
	
	infoTmpList = [[0.0]*infoN for i in xrange(CFCR.insertN)]
	infoList = [deepcopy(infoTmpList) for i in xrange(len(execNames))]
	cntList = [0] * len(execNames)
	for eid,execName in enumerate(execNames):
		for caseId in range(caseN):
			tmpFileName = "data_%d.txt" % (caseId)
			resFileName = os.path.join(srcFilePath, execName, resFilePath, tmpFileName)
			if not os.path.exists(resFileName):
				print "%s has no log for %s/%s" % (execName, resFilePath, tmpFileName)
				continue
			tmpList = fetchResult(resFileName)
			for i in xrange(CFCR.insertN):
				for k in xrange(infoN):
					infoList[eid][i][k] += tmpList[i][k]
			cntList[eid] += 1
	tmpLines = genBlock(execNames, infoList, cntList, bid)
	lines += tmpLines
	
	return lines


def calcResult(srcFilePath, desFileName):
	alpha = 2
	
	## initialization
	lines = []
	bid, step, vspace = CFCR.begBid, CFCR.infoN, ["\n\n"]

	lines += ["%% the synthetic dataset of Uniform\n"]
	resFilePath =  "uni_%d" % (alpha)
	tmpLines = calcResultByVariable(srcFilePath, resFilePath, bid)
	lines += tmpLines + vspace
	bid += step

	lines += ["%% the synthetic dataset of Normal\n"]
	resFilePath =  "nor_%d" % (alpha)
	tmpLines = calcResultByVariable(srcFilePath, resFilePath, bid)
	lines += tmpLines + vspace
	bid += step

	lines += ["%% the synthetic dataset of Exponential\n"]
	resFilePath =  "exp_%d" % (alpha)
	tmpLines = calcResultByVariable(srcFilePath, resFilePath, bid)
	lines += tmpLines + vspace
	bid += step
	
	lines += ["%% the real dataset of TKY\n"]
	resFilePath = "checkinTKY_%d" % (alpha)
	tmpLines = calcResultByVariable(srcFilePath, resFilePath, bid)
	lines += tmpLines + vspace
	bid += step

	lines += ["%% the real dataset of Haikou\n"]
	resFilePath = "gaiaHaikou_%d" % (alpha)
	tmpLines = calcResultByVariable(srcFilePath, resFilePath, bid)
	lines += tmpLines + vspace
	bid += step
	
	if desFileName:
		with open(desFileName, "w") as fout:
			fout.writelines(lines)

	return lines
	

def exp0():
	srcFilePath = "./updateResult"
	desFileName = "./updateResult.res"
	calcResult(srcFilePath, desFileName)


if __name__ == "__main__":
	exp0()
