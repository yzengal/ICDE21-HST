import os
import sys
import math
from random import randint, sample, shuffle


class constForReal():
	caseN = 5
	# pointFiles = ["checkinTKY.txt", "gaiaHaikou.txt"]
	pointFiles = ["uni.txt", "nor.txt", "exp.txt"]
	alphas = [2]
	insertN = 10

class CFR(constForReal):
	pass


def genRealData(desPath, srcFile, caseN):
	n = 0
	points = []
	with open(srcFile, "r") as fin:
		for line in fin:
			line = line.strip()
			if not line:
				continue
			if n == 0:
				n = int(line)
			else:
				tmpList = map(int, line.split())
				points.append(tmpList)


	n = len(points)
	shuffle(points)
	perms = range(0, n)
	prefix = srcFile[:-4]
	for k in CFR.alphas:
		tmpFilePath = "%s_%d" % (prefix, k)
		tmpFilePath = os.path.join(desPath, tmpFilePath)
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)
		for i in xrange(caseN):
			desFileName = "data_%d.txt" % (i)
			desFileName = os.path.join(tmpFilePath, desFileName)
			print("start make %s" % desFileName)
			opt_cnt = 0
			with open(desFileName, "w") as fout:
			# shuffle(perms)
			# beta = randint(10000, k * 10000) * 1.0 / (k * 10000)
			# fout.write("%d %.4f %d\n" % (n, beta, k))
			# line = " ".join(map(str, perms)) + "\n"
			# fout.write(line)
			# for i in xrange(n):
			#	fout.write("%d %d\n" % (points[i][0], points[i][1]))
				# shuffle(perms)
				lines = []
				
				static_coef = 0.8
				id_prob = 50
				ins_max_coef = 0.03
				ins_min_coef = 0.01
				dec_min_coef = 0.01
				dec_max_coef = 0.5
			# generate config

				ins_max = int(math.ceil(n * ins_max_coef))
				ins_min = int(math.ceil(n * ins_min_coef))

				static_n = int(math.ceil(n * static_coef))

			# get the buffer
				static_set = perms[0:static_n]
				buffer = perms[static_n : n]
				buffer_size = len(buffer)
				inserted = []
				
			# static insert
				opt_cnt += 1
				lines.append("I %d\n" % static_n)
				line = ""
				for i in static_set:
					line += "%d %s\n" % (i, " ".join(map(str, points[i])))
				lines.append(line)
				
			# pre-calculate the number of inserted points in each ten insertion
				insertNums = [1] * CFR.insertN
				left = buffer_size - CFR.insertN
				for iii in xrange(CFR.insertN):
					more = randint(ins_min, ins_max)
					more = min(more, left)
					insertNums[iii] += more
					left -= more
				insertNums[CFR.insertN-1] += left
				insertId = 0
					
			# dynamic insert and decrease
				while buffer_size != 0:
					op = randint(1,100)
					if op <= id_prob: # insert
						opt_cnt += 1
						dy_num = randint(ins_min, ins_max)
						dy_num = min(dy_num, buffer_size)
						dy_num = insertNums[insertId]
						insertId += 1
						# shuffle(buffer)
						lines.append("I %d\n" % (dy_num))
						line = ""
						for i in buffer[0:dy_num]:
							line += "%d %s\n" % (i, " ".join(map(str, points[i])))
						lines.append(line)
						inserted += buffer[0:dy_num]
						buffer = buffer[dy_num:]
						buffer_size = len(buffer)

					elif len(inserted) > 0: # delete
						opt_cnt += 1
						l = len(inserted)
						dec_min = math.floor(l * dec_min_coef)
						dec_max = math.ceil(l * dec_max_coef)
						dy_num = randint(dec_min, dec_max)
						dy_num = min(dy_num, len(inserted))
						shuffle(inserted)
						lines.append("D %d\n" % dy_num)
						line = ""
						for i in inserted[0:dy_num]:
							line += "%d " % (i)
						lines.append(line + "\n")
						# buffer += inserted[0:dy_num]
						inserted = inserted[dy_num:]

			# global parameter
				fout.write("%d %d %d\n" % (n, opt_cnt, k))
				fout.writelines(lines)


def genRealDataSet(desPath, caseN):
	if not os.path.exists(desPath):
		os.mkdir(desPath)
	for srcFile in CFR.pointFiles:
		genRealData(desPath, srcFile, caseN)


def exp0():
	desPath = "updateData_zyx"
	genRealDataSet(desPath, CFR.caseN)


if __name__ == "__main__":
	exp0()
