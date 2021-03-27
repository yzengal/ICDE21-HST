/**
	\Author: Trasier
	\Date:	2019/10/9
**/

#ifndef HST_H
#define HST_H

#include <bits/stdc++.h>
using namespace std;
#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

//#define LOCAL_DEBUG
#define USE_TIMER

typedef struct Node_t {
	int nid;		// the idx of the node
	int nomIdx; 	// the id of norminator point
	int lev;		// the level of the node before compression
	Node_t *far;	// the parent node
	float wei;		// the weight of the edge from its parent
	int cid;		// the position of this node in its parent's child list.
	vector<Node_t*> child;
	
	Node_t() {}
	Node_t(int nid_, int nomIdx_, int lev_, Node_t* far_, float wei_) {
		nid = nid_;
		nomIdx = nomIdx_;
		lev = lev_;
		far = far_;
		wei = wei_;
		cid = 0;
	}
} ;

extern const int MAX_HEIGHT;
extern int H;
extern double* expks;
extern double* sumks;
extern double logAlpha;
extern int* pi;
extern int* reverse_pi;
extern double dmax;
extern double gamma_or;
extern Node_t* rt;
extern Node_t** leaves;


void initLocation(string &fileName);
void initMemory_HST(int n);
void freeMemory_HST();
void freeHST(Node_t*& rt);
void addChild(Node_t* far, Node_t* u);
void constructHST_fast2(bool load, clock_t startClock);
void constructHST_fast(bool load, clock_t startClock);
void constructHST_slow(bool load, clock_t startClock);
void constructHST_fast_am(bool load, clock_t startClock);
void randomization();
void calcDmax();
double distAtLevel(int level);
double distOnHST(int u, int v);
double distOnHST(Node_t* u, Node_t* v);
pair<Node_t*, int> getLCA(int u, int v);
pair<Node_t*, int> getLCA(Node_t* u, Node_t* v);
int levelOfLCA(int u, int v);
int levelOfLCA(Node_t* u, Node_t* v);

#endif
