##
#	Makefile for the manuscript submitted to ICDE'21,
#	Author: Trasier
##
CC = gcc
CXX = g++
LIBS = 
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium

all: BF ODP OPT
monitor.o: monitor.h monitor.cpp
	$(CXX) $(CFLAGS) -c monitor.cpp $(LIBS) -D WATCH_MEM
	
global.o: global.h global.cpp
	$(CXX) $(CFLAGS) -c global.cpp $(LIBS) -D WATCH_MEM
	
HST.o: HST.h HST.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST.cpp $(LIBS) -D WATCH_MEM

HST_opt.o: HST_opt.h HST_opt.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_opt.cpp $(LIBS) -D WATCH_MEM

BF: BF.cpp global.h global.cpp global.o HST.h HST.cpp HST.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o BF monitor.o global.o HST.o BF.cpp $(LIBS) -D WATCH_MEM

OPT: OPT.cpp global.h global.cpp global.o HST_opt.h HST_opt.cpp HST_opt.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o OPT monitor.o global.o HST_opt.o OPT.cpp $(LIBS) -D WATCH_MEM

ODP: ODP.cpp global.h global.cpp global.o HST_opt.h HST_opt.cpp HST_opt.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o ODP monitor.o global.o HST_opt.o ODP.cpp $(LIBS) -D WATCH_MEM
	
.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm BF ODP OPT 2> /dev/null || true
		
