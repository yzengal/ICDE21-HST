/**
	\author:	Trasier
	\date:		2019.6.18
*/
#include "global.h"

int nV = 0;
location_t* V = NULL;
const double speed = 1.0;
const double EPS = 1e-5;
const double INF = 1e20;
// const int MAX_DIM = 20;
double usedTime = 0.0;
int usedMemory = 0;
double timeLimit = 10 * 24 * 60 * 60; // 10 days

int dcmp(double x) {
	if (fabs(x) < EPS)
		return 0;
	return x>0 ? 1:-1;
}

double dist(location_t *V, int x, int y) { 
	if (x == y) return 0;
	
	location_t& a = V[x];
	location_t& b = V[y];
	double ret = 0.0; 
	
	for (int i=0; i<MAX_DIM; ++i) {
		ret += (a.x[i]-b.x[i])*(a.x[i]-b.x[i]);
	}
	
	return sqrt(ret);
} 

double dist(location_t &a, location_t &b) {
	double ret = 0.0; 
	
	for (int i=0; i<MAX_DIM; ++i) {
		ret += (a.x[i]-b.x[i])*(a.x[i]-b.x[i]);
	}
	
	return sqrt(ret);
}

void freeLocation() {
	delete[] V;
}
