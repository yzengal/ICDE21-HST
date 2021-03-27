#!/usr/bin/env python
import os
import sys
import datetime


class constForCalcResult:
	caseN = 5
	alphas = [2,3,4,5,6]
	execNames = ["BF", "ODP", "OPT"]
	infoN = 2
	infos = ["runtime", "memusage"]
	begBid = 1
	endBid = 8
	
class CFCR(constForCalcResult):
	pass


def parseResult(resFileName, line):
	tmpList = line.strip().split(" ")
	if len(tmpList) != 3:
		print "result is wrong: %s: %s" % (resFileName, line)
		return None
	runtime, memused = tmpList[-2:]
	runtime, memused = float(runtime), float(memused)/(1024*1024)
	return runtime, memused


def fetchResult(fileName):
	ret = ""
	with open(fileName, "r") as fin:
		for line in fin:
			line = line.strip()
			if line:
				ret = line
	return ret


def genBlock(execNames, infoAll, cntAll, bid):
	lines = []
	
	for iid,infoName in enumerate(CFCR.infos):
		tmpLines = []
		tmpLines.append( "%%%%%%%% %s\n" % (infoName) )
		for eid,execName in enumerate(execNames):
			tmpList = []
			for vid in xrange(len(infoAll)):
				infoList = infoAll[vid][eid]
				cntList = cntAll[vid]
				tot = infoList[iid]
				cnt = cntList[eid]
				if cnt == 0:
					avg = 0.0
				else:
					avg = tot*1.0/cnt
				tmpList.append(avg)
			tmpLine = ", ".join(map(lambda s:"%.2f" % (s), tmpList))
			line = "%s%d = [%s];\n" % (execName, bid, tmpLine)
			tmpLines.append(line)
		lines += tmpLines
		bid += 1
	return lines


def calcResultByVariable(srcFilePath, resFilePaths, bid):
	execNames = CFCR.execNames
	caseN = CFCR.caseN
	infoN = CFCR.infoN
	lines = []
	
	tmpList1, tmpList2 = [], []
	for resFilePath in resFilePaths:
		infoList = [[0.0]*infoN for i in xrange(len(execNames))]
		cntList = [0] * len(execNames)
		for eid,execName in enumerate(execNames):
			for caseId in range(caseN):
				tmpFileName = "data_%d.dat" % (caseId)
				resFileName = os.path.join(srcFilePath, execName, resFilePath, tmpFileName)
				if not os.path.exists(resFileName):
					print "%s has no log for %s/%s" % (execName, resFilePath, tmpFileName)
					continue
				line = fetchResult(resFileName)
				tmpList = parseResult(resFileName, line)
				if not tmpList:
					continue
				# print tmpList
				for k in xrange(infoN):
					infoList[eid][k] += tmpList[k]
				cntList[eid] += 1
		tmpList1.append(infoList)
		tmpList2.append(cntList)
	tmpLines = genBlock(execNames, tmpList1, tmpList2, bid)
	lines += tmpLines
	
	return lines


def calcResult(srcFilePath, desFileName):
	## initialization
	lines = []
	bid, step, vspace = CFCR.begBid, CFCR.infoN, ["\n\n"]

	# the NYC real datset
	lines += ["%% the real dataset of NYC\n"]
	resFilePaths = []
	for alpha in CFCR.alphas:
		tmpFilePath = "checkinNYC_%d" % (alpha)
		resFilePaths.append(tmpFilePath)
	tmpLines = calcResultByVariable(srcFilePath, resFilePaths, bid)
	lines += tmpLines + vspace
	bid += step
	
	# the TYK real datset
	lines += ["%% the real dataset of TKY\n"]
	resFilePaths = []
	for alpha in CFCR.alphas:
		tmpFilePath = "checkinTKY_%d" % (alpha)
		resFilePaths.append(tmpFilePath)
	tmpLines = calcResultByVariable(srcFilePath, resFilePaths, bid)
	lines += tmpLines + vspace
	bid += step

	# the Haikou real datset
	lines += ["%% the real dataset of Haikou\n"]
	resFilePaths = []
	for alpha in CFCR.alphas:
		tmpFilePath = "gaiaHaikou_%d" % (alpha)
		resFilePaths.append(tmpFilePath)
	tmpLines = calcResultByVariable(srcFilePath, resFilePaths, bid)
	lines += tmpLines + vspace
	bid += step
	
	# the Haikou real datset
	lines += ["%% the real dataset of Chengdu\n"]
	resFilePaths = []
	for alpha in CFCR.alphas:
		tmpFilePath = "gaiaChengdu_%d" % (alpha)
		resFilePaths.append(tmpFilePath)
	tmpLines = calcResultByVariable(srcFilePath, resFilePaths, bid)
	lines += tmpLines + vspace
	bid += step
	
	if desFileName:
		with open(desFileName, "w") as fout:
			fout.writelines(lines)

	return lines
	

def exp0():
	srcFilePath = "./realResult"
	desFileName = "./realResult.res"
	calcResult(srcFilePath, desFileName)


if __name__ == "__main__":
	exp0()
