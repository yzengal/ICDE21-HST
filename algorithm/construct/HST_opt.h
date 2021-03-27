/**
	\Author: Trasier
	\Date:	2019/06/18
**/

#ifndef HST_OPT_H
#define HST_OPT_H

#include <bits/stdc++.h>
#include "global.h"
using namespace std;
#ifdef WATCH_MEM
#include "monitor.h"
#endif

//#define LOCAL_DEBUG
#define USE_TIMER

extern const int MAX_SAMPLE;
extern const int MAX_HEIGHT;
extern int H;
extern int alpha;
extern double* expks;
extern double* sumks;
extern int* pi;
extern int* reverse_pi;
extern double dmax;
extern double beta;


struct Node_t{
	int nid;		// the idx of the node
	int nomIdx; 	// the id of norminator point
	int lev;		// the level of the node before compression
	Node_t *far;	// the parent node
	int cid;		// the position of this node in its parent's child list.
	float wei;		// the weight of the edge from its parent 
	vector<Node_t*> child;
	
	Node_t(int nid_ = 0, int nomIdx_=0, int lev_=0, Node_t* far_=NULL, float wei_=0) {
		nid = nid_;
		nomIdx = nomIdx_;
		lev = lev_;
		far = far_;
		wei = wei_;
		cid = 0;
	}
} ;

extern Node_t* rt;
extern Node_t** leaves;

void mergeChild(Node_t* far, Node_t* u);
void addChild(Node_t* far, Node_t* u);
void initLocation(string &fileName);
void initMemory_HST(int n);
void freeMemory_HST();
void freeHST(Node_t* &rt);
void constructHST_fast_dp(bool load, clock_t startClock);
void constructHST_fast_opt(bool load, clock_t startClock);
void randomization();
void calcDmaxPrune();
double distAtLevel(int level);
double distOnHST(int u, int v);
double distOnHST(Node_t* u, Node_t* v);
pair<Node_t*, int> getLCA(int u, int v);
pair<Node_t*, int> getLCA(Node_t* u, Node_t* v);
int levelOfLCA(int u, int v);
int levelOfLCA(Node_t* u, Node_t* v);


#endif
