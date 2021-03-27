/**
	\Author: Trasier
	\Date:	2019/10/9
**/
#ifndef HST_DPO_H
#define HST_DPO_H


#include <bits/stdc++.h>
using namespace std;

#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

struct Node_t {
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
};

struct Tree_t {
	int nV;
	int H;
	int node_n, merge_n;
	double logAlpha;
	double dmax;
	double beta;
	double alpha;
	double distor, distor_, cdistor;
	double* expks;
	double* sumks;
	Node_t* rt;
	vector<int> pi;
	vector<int> reverse_pi;
	vector<location_t> V;
	vector<Node_t*> leaves;
	vector<bool> deleted; // mark whether the points are deleted
	
public:
	void initLocation(int nV, vector<location_t>& V, double alpha);  // init the parameter of the tree 
	void constructHST(clock_t startClock); // startClock is used to control whether the algorithm is TLE.
	void freeMemory(); // free the memory usage of the structure
	pair<double,double> getDistortion(); // get the distortion of the tree						
	double distOnHST(int u, int v); // query the distance on the tree
	pair<Node_t*, int> getLCA(int u, int v); // query the LCA of two boths
	void deleteOne(int idx);
	long long  memoryCost();
	
private:
	void initMemory(int n);
	void randomization();
	void freeHST(Node_t*& rt);
	double distAtLevel(int level);
	double distOnHST(Node_t* u, Node_t* v);
	pair<Node_t*, int> getLCA(Node_t* u, Node_t* v);
	int levelOfLCA(int u, int v);
	int levelOfLCA(Node_t* u, Node_t* v);
	
	inline double log2(double x) {
		return log(x) / logAlpha;
	}

	inline double pow2(int i) {
		return expks[i];
	}
	
	inline int getLevel(double H, double dist) {
		if (dist < 1.0) return H+1;
		
		int k = ceil(log2(dist/beta));
		
		if (expks[k]*beta == dist)
			++k;
		
		return H+1-k;
	}
	
	inline void addChild(Node_t* far, Node_t* u) {
		u->far = far;
		u->cid = far->child.size();
		far->child.push_back(u);
	}

	inline void mergeChild(Node_t* far, Node_t* u) {
		if (far == rt || far->child.size()!=1) 
			return ;
		
		Node_t* gfar = far->far;
		int gcid = far->cid;
		
		u->far = far->far;
		u->cid = gcid;
		u->wei += far->wei;
		
		gfar->child[gcid] = u;
		delete far;
		merge_n++;
	}
};

#endif
