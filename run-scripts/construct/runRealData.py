import os
import sys
import random
import commands
import multiprocessing

execNames = ["ODP", "OPT", "BF"]

def run(execName, dataFile, logFile):
	cmdLine = "./%s %s" % (execName, dataFile)
	print(cmdLine)
	line = commands.getoutput(cmdLine)
	print line
	print "[Finish]", cmdLine
	with open(logFile, "w") as fout:
		fout.write(line)


def batchRun(srcFilePath, desFilePath, dataSetN, nprocess):
	pool = multiprocessing.Pool(processes = nprocess)
	
	if not os.path.exists(desFilePath):
		os.mkdir(desFilePath)
	dirNames = os.listdir(srcFilePath)
	for dataSetId in range(dataSetN):
		for execName in execNames:
			tmpFilePath = os.path.join(desFilePath, execName)
			if not os.path.exists(tmpFilePath):
				os.mkdir(tmpFilePath)
			for dirName in dirNames:
				tmpFilePath = os.path.join(desFilePath, execName, dirName)
				if not os.path.exists(tmpFilePath):
					os.mkdir(tmpFilePath)
				srcFileName = "data_%d.dat" % (dataSetId)
				desFileName = os.path.join(desFilePath, execName, dirName, srcFileName)
				srcFileName = os.path.join(srcFilePath, dirName, srcFileName)
				if os.path.exists(desFileName):
					continue
				pool.apply_async(run, (execName, srcFileName, desFileName, ))

	pool.close()
	pool.join()

def test():
	dataSetN = 1
	nprocess = 1
	srcFilePath = "./realData"
	desFilePath = "./realResult"
	batchRun(srcFilePath, desFilePath, dataSetN, nprocess)

if __name__ == '__main__':
	test()