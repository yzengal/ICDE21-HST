##
#	Makefile for the manuscript submitted to ICDE'21,
#	Author: Trasier
##
CC = gcc
CXX = g++
LIBS = 
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium

all: reHSTO HSF
monitor.o: monitor.h monitor.cpp
	$(CXX) $(CFLAGS) -c monitor.cpp $(LIBS) -D WATCH_MEM
	
global.o: global.h global.cpp
	$(CXX) $(CFLAGS) -c global.cpp $(LIBS) -D WATCH_MEM

HST_dpo.o: HST_dpo.h HST_dpo.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_dpo.cpp $(LIBS) -D WATCH_MEM	
	
reHSTO: reHSTO.cpp global.h global.cpp global.o HST_dpo.h HST_dpo.cpp HST_dpo.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o reHSTO monitor.o global.o HST_dpo.o reHSTO.cpp $(LIBS) -D WATCH_MEM
	
HST_opt.o: HST_opt.h HST_opt.cpp global.h global.cpp global.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -c global.o monitor.o HST_opt.cpp $(LIBS) -D WATCH_MEM
		
HSF: HSF.cpp global.h global.cpp global.o HST_opt.h HST_opt.cpp HST_opt.o monitor.h monitor.cpp monitor.o
	$(CXX) $(CFLAGS) -o HSF monitor.o global.o HST_opt.o HSF.cpp $(LIBS) -D WATCH_MEM
	
	
.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm reHSTO HSF 2> /dev/null || true
		
