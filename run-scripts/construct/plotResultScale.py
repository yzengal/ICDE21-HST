#!/usr/bin/env python
import commands
import os
import sys
from random import randint
from calcResultScale import constForCalcResult, calcResult


class constForPlotResult(constForCalcResult):
	fontsz = 15
	marksz = 13
	legendsz = 15
	linewd = 3
	vspace = "\n\n"
	figWidth = 5
	figHeight = 3
	labelsize = 15
	colorN = 10
	lineStyles = {
		"BF":'s-', 
		"ODP":'^-', 
		"OPT":'o-',
		"SPAA":'+-',
	}
	lineColors = {
		"BF":3, 
		"ODP":2, 
		"OPT":1,
		"SPAA":5,
	}
	lineNames = {
		"BF":"BASE", 
		"ODP":"HST+DP", 
		"OPT":"HST+DPO",
		"SPAA":"SPAA",
	}
	"runtime", "memusage"
	xLabels = ["n"]
	yLabels = ["Running time (Secs)", "Memory usage (MB)"]
	
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
	lines.append("set(gca, 'XLim', [min(x%d), max(x%d)]);\n" % (bid, bid))
	if not ylabel.startswith("Memory"):
		lines.append("set(gca, 'YScale', 'log');\n")
	lines.append("set(gca, 'FontName', 'Arial');\n")
	lines.append("xlabel('%s', 'FontSize', labelsize);\n" % (xlabel))
	lines.append("ylabel('%s', 'FontSize', labelsize);\n" % (ylabel))

	# add legend
	lineNames = map(lambda s: CFPR.lineNames[s], CFPR.execNames)
	tmpList = map(lambda s:"'%s'" % (s), lineNames)
	line = "h_legend = legend(%s, 'Location', 'SouthEast');\n" % (", ".join(tmpList))
	lines.append(line)
	lines.append("set(h_legend,'color', 'none');\n")
	lines.append("set(h_legend, 'FontSize', legendsize);\n")
	lines.append("set(h_legend,'Box','off');\n")
	lines.append("hold off;\n\n")
	
	return lines


def genFigureBlocks(bid, eid, xLists):
	ret = []
	for i in xrange(bid, eid+1):
		xList = xLists[(i-bid)/CFPR.infoN]
		xlabel = CFPR.xLabels[(i-bid)/CFPR.infoN]
		ylabel = CFPR.yLabels[(i-bid)%CFPR.infoN]
		lines = genFigureBlock(i, xList, xlabel, ylabel)
		ret += lines
	return ret
	
def genFigure(srcFilePath, desFileName):
	xLists = []
	xLists.append(map(str, map(int, CFPR.nList)))
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
	srcFilePath = "./scaleResult"
	desFileName = "./scaleResult.m"
	genFigure(srcFilePath, desFileName)
	
	
if __name__ == "__main__":
	exp0()
	
