#!/usr/bin/env python
import commands
import os
import sys
from random import randint
from calcResultUpdate import constForCalcResult, calcResult

class constForPlotResult(constForCalcResult):
	fontsz = 11
	marksz = 11
	legendsz = 10
	linewd = 2
	vspace = "\n\n"
	figWidth = 5
	figHeight = 3
	labelsize = 14
	colorN = 10
	lineStyles = {
		"reHSTO":'s-', 
		"HSF":'o-',
	}
	lineColors = {
		"reHSTO":3, 
		"HSF":1,
	}
	lineNames = {
		"reHSTO":"reHST", 
		"HSF":"HST+HSF",
	}
	xLabels = ["Uniform (No. of insertion)", "Normal (No. of insertion)", "Exponential (No. of insertion)", "Fsq-TKY (No. of insertion)", "Didi-Haikou (No. of insertion)"]
	yLabels = ["distortion", "running time (Secs)", "memory usage (MB)"]
	
class CFPR(constForPlotResult):
	pass
	

def genHead():
	lines = []
	lines.append("fontsize = %s;\n" % (CFPR.fontsz))
	lines.append("markersize = %s;\n" % (CFPR.marksz))
	lines.append("legendsize = %s;\n" % (CFPR.legendsz))
	lines.append("linewidth = %s;\n" % (CFPR.linewd))
	lines.append("figwidth = %s;\n" % (CFPR.figWidth))
	lines.append("figheight = %s;\n" % (CFPR.figHeight))
	lines.append("labelsize = %s;\n" % (CFPR.labelsize))
	lines.append(CFPR.vspace)
	return lines

	
def genFigureBlock(bid, xList, xlabel, ylabel):
	lines = []
	# lines.append("%%%%%02d---------------------------------------------------------\n" % (bid))
	lines.append("figure;\n")
	lines.append("hold on;\n")
	lines.append("figuresize(%s, %s, 'inches');\n" % ("figwidth", "figheight"))
	lines.append("box on;\n")
	# append the X&Y
	line = "x%d = [%s];\n" % (bid, ", ".join(xList))
	lines.append(line)
	line = "mpdc_ = distinguishable_colors(%d);\n" % (CFPR.colorN)
	lines.append(line)
	for i in xrange(len(CFPR.execNames)):
		execName = CFPR.execNames[i]
		line = "plot(x%d, %s%d, '%s', 'Color', mpdc_(%s,:), 'LineWidth', linewidth, 'MarkerSize', markersize);\n" % (bid, execName, bid, CFPR.lineStyles[execName], CFPR.lineColors[execName])
		lines.append(line)
	# add xlabel & ylabel
	lines.append("set(gca, 'XTickLabel', x%d);\n" % (bid))
	lines.append("rotateXLabels(gca(), 60);\n")
	lines.append("set(gca, 'FontSize', fontsize, 'Xtick', x%d);\n" % (bid))
	if "distortion" in ylabel:
		lines.append("set(gca, 'YLim', [0.0, max(reHSTO%d)+2000]);\n" % (bid))
	elif "memory" in ylabel:
		if xlabel.startswith("Didi"):
			lines.append("set(gca, 'YLim', [10.0, 30.0]);\n")
		else:
			lines.append("set(gca, 'YLim', [0.0, 10.0]);\n")
	lines.append("set(gca, 'XLim', [min(x%d), max(x%d)]);\n" % (bid, bid))
	if "time" in ylabel:
		lines.append("set(gca, 'YScale', 'log');\n")
	lines.append("xlabel('%s', 'FontSize', labelsize);\n" % (xlabel))
	lines.append("ylabel('%s', 'FontSize', labelsize);\n" % (ylabel))

	# add legend
	lines.append("set(gca, 'FontName', 'Arial');")
	lineNames = map(lambda s: CFPR.lineNames[s], CFPR.execNames)
	tmpList = map(lambda s:"'%s'" % (s), lineNames)
	line = "h_legend = legend(%s, 'Location', 'SouthEast');\n" % (", ".join(tmpList))
	lines.append(line)
	lines.append("set(h_legend,'color', 'none');\n")
	lines.append("set(h_legend, 'FontSize', legendsize);\n")
	lines.append("hold off;\n\n")
	# if bid==4:
		# lines.append("return;\n\n")
	
	return lines


def genFigureBlocks(bid, eid, xLists):
	ret = []
	for i in xrange(bid, eid+1):
		# print i, i-bid, CFPR.infoN, len(xLists)
		xList = xLists[(i-bid)/CFPR.infoN]
		xlabel = CFPR.xLabels[(i-bid)/CFPR.infoN]
		ylabel = CFPR.yLabels[(i-bid)%CFPR.infoN]
		lines = genFigureBlock(i, xList, xlabel, ylabel)
		ret += lines
	return ret
	
def genFigure(srcFilePath, desFileName):
	xLists = []
	xLists.append(map(str, range(1,CFPR.insertN+1)))
	xLists.append(map(str, range(1,CFPR.insertN+1)))
	xLists.append(map(str, range(1,CFPR.insertN+1)))
	xLists.append(map(str, range(1,CFPR.insertN+1)))
	xLists.append(map(str, range(1,CFPR.insertN+1)))
	bid, eid = CFPR.begBid, CFPR.endBid
	with open(desFileName, "w") as fout:
		### head
		lines = genHead()
		fout.writelines(lines)
		lines = calcResult(srcFilePath, "")
		fout.writelines(lines)
		fout.write(CFPR.vspace)
		lines = genFigureBlocks(bid, eid, xLists)
		fout.writelines(lines)


def exp0():
	srcFilePath = "./updateResult"
	desFileName = "./updateResult.m"
	genFigure(srcFilePath, desFileName)
	
	
if __name__ == "__main__":
	exp0()
	