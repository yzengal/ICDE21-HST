/**
	\Author: Trasier
	\Date:	2019/06/18
**/

#include "HST_opt.h"
#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

const int MAX_HEIGHT = 40;
int H = 0;
int alpha = 2;
double logAlpha = 1.0;
int* pi = NULL;
int* reverse_pi = NULL;
double dmax = -1.0;
double* expks = NULL;
double* sumks = NULL;
Node_t* rt = NULL;
Node_t** leaves = NULL;
double beta = 1.0;

inline double log2(double x) {
	return log(x) / logAlpha;
}

inline double pow2(int i) {
	return (i<0) ? 1.0 : expks[i];
}

inline int getLevel(double H, double dist) {
	if (dist < 1.0) return H+1;
	
	int k = ceil(log2(dist/beta));
	
	if (expks[k]*beta == dist)
		++k;
	
	return H+1-k;
}

void initMemory_HST(int n) {
	nV = n;
	V = new location_t[nV];
	pi = new int[nV];
	reverse_pi = new int[nV];
	expks = new double[MAX_HEIGHT];
	expks[0] = 1.0;
	for (int i=1; i<MAX_HEIGHT; ++i)
		expks[i] = expks[i-1] * alpha;
	
	sumks = new double[MAX_HEIGHT+1];
	sumks[0] = expks[0];
	for (int i=1; i<=MAX_HEIGHT; ++i)
		sumks[i] = sumks[i-1] + expks[i];
	
	leaves = new Node_t*[nV];
	
	logAlpha = log(alpha);
}

void freeMemory_HST() {
	freeLocation();
	delete[] pi;
	delete[] reverse_pi;
	delete[] expks;
	delete[] sumks;
    delete[] leaves;
	freeHST(rt);
}

void freeHST(Node_t*& rt) {
	if (rt == NULL) return ;
	
	Node_t* chi;
	
	for (int i=0; i<rt->child.size(); ++i) {
		chi = rt->child[i];
		freeHST(chi);
	}
	
	delete rt;
	rt = NULL;
}

void initLocation(string &fileName) {
	ifstream fin(fileName.c_str(), ios::in);   
    
	if (!fin.is_open()) {
		fprintf(stderr, "FILE %s IS INVALID.", fileName.c_str());
		exit(1);
	}
	
	fin >> nV;
    fin >> beta;
    fin >> alpha;
    string a;
	initMemory_HST(nV);
    for(int i=0;i<nV;i++){
        int tmp = i;
#ifndef LOCAL_DEBUG
        fin >> tmp;
#endif
        pi[i] = tmp;
        reverse_pi[tmp] = i;
    }

    for (int i=0; i<nV; ++i) {
		for (int j=0; j<MAX_DIM; ++j)
			fin >> V[i].x[j];
    }
    fin.close();
}

void randomization() {
	// generate the permutation pi
	for (int i=0; i<nV; ++i){
		pi[i] = i;
	}
	random_shuffle(pi, pi+nV);
	for (int i=0; i<nV; ++i){
		reverse_pi[pi[i]] = i;
	}
	// generate the parameter gamma
	beta = (rand()%alpha + 1) * 1.0	/ alpha;
}

inline void addChild(Node_t* far, Node_t* u) {
	u->far = far;
	u->cid = far->child.size();
	far->child.push_back(u);
}

static int merge_n;

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

void constructHST_fast_opt(bool load, clock_t startClock) {	
	// label the far[][] with node ID
	int nid = 1;
	int *cnt = NULL, *cen = NULL, *next_cen = NULL;
	Node_t** nodes = NULL;
	Node_t** nodes_ = NULL; 
	double *cenl = NULL;
	int *p = NULL, *q = NULL, *far = NULL, *far_ = NULL;
	
	// cnt[i]: count
	// p: x
	// q: 
	cnt = new int[nV];
	cen = new int[nV];
	cenl = new double[nV];
	p = new int[nV];
	q = new int[nV];
	far = new int[nV];
	far_ = new int[nV];
	nodes = new Node_t*[nV];
	nodes_ = new Node_t*[nV];
	
	if (rt != NULL)
		freeHST(rt);
	if (!load)
		randomization();

	// Prunning 1: calcDmax()
	dmax = 0;
	for (int i=0; i<nV; ++i) {
		cen[i] = 0;
		cenl[i] = dist(V, i, pi[cen[i]]);
		dmax = max(dmax, cenl[i]);
	}
	
	// initialization
	H = ceil(log2(dmax+EPS)) + 1;	
	double radius = pow2(H)*beta;
	
	// construct the root
	rt = new Node_t(0, pi[0], 1, NULL, radius);
	for (int i=0; i<nV; ++i) {
        nodes_[i] = rt;
		far[i] = far_[i] = 0;
		p[i] = i;
	}
	
	merge_n = 0;
	for (int k=2; k<=H+1; ++k) {
		radius /= alpha;
		for (int i=0; i<nV; ++i) {
			if (cenl[i] < radius)
				continue;
			
			int pid;
			while (cenl[i] >= radius) {
				pid = pi[++cen[i]];
				if (cen[pid] <= reverse_pi[i]) {
					cenl[i] = dist(V, i, pi[cen[i]]);
				}
			}
		}
		
		int bid = nid, i = 0, bi;
		memset(q, -1, sizeof(int)*nV);
		while (i < nV) {
			q[cen[p[i]]] = far[p[i]] = nid++;
			bi = i++;
			while (i<nV && far_[p[i]]==far_[p[i-1]]) {
				int pid = cen[p[i]];
				if (q[pid] == -1) {
					q[pid] = nid++;
				}
				far[p[i]] = q[pid];
				++i;
			}
			while (bi < i) {
				int pid = cen[p[bi]];
				q[pid] = -1;
				++bi;
			}
		}
		memset(cnt, 0, sizeof(int)*nV);
		for (int i=0; i<nV; ++i) {
			++cnt[far[i]-bid];
		}
		for (int j=1; j<nid-bid; ++j) {
			cnt[j] += cnt[j-1];
		}
		for (int i=nV-1; i>=0; --i) {
			int j = far[p[i]] - bid;
			q[--cnt[j]] = p[i];
		}
		// create the new node at the $k$-th level
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			nodes[q[i]] = new Node_t(far[q[i]], pi[cen[q[i]]], k, nodes_[q[i]], radius);
			addChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				nodes[q[j]] = nodes[q[i]];
				++j;
			}
		}
		// merge the new node with its parent 
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			mergeChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				++j;
			}
		}
		// fill in the leaves
		if (k == H+1) {
			for (int i=0; i<nV; ++i) {
				leaves[q[i]] = nodes[q[i]];
			}
		}
		swap(far, far_);
		swap(p, q);
		swap(nodes, nodes_);

#ifdef USE_TIMER
		if (k%2 == 1) {
			double usedTime = (clock() - startClock) / CLOCKS_PER_SEC;
			if (usedTime > timeLimit) {
				printf("Finished the %d-th level, %d left.\n", k, H+1-k);
				printf("Time Limitation Exceed\n");
				break; 
			}
		}
#endif
	}
	
	usedMemory = 0;
	usedMemory += (nid-merge_n) * sizeof(Node_t);
	printf("nV=%d, nid=%d, H=%d,  merge_n=%d\n", nV, nid, H, merge_n);
	usedMemory += nV*2*sizeof(int);
	
	delete[] cnt;
	delete[] cen;
	delete[] cenl;
	delete[] p;
	delete[] q;
	delete[] far;
	delete[] far_;
	delete[] nodes;
	delete[] nodes_;
}

double distAtLevel(int lev) {
	if (lev >= H+1) return 0.0;
	return sumks[H-lev] * beta * 2.0;

}

double distOnHST(int u, int v) {
	int level = levelOfLCA(u, v);
	return distAtLevel(level);
}

double distOnHST(Node_t* u, Node_t* v) {
	int level = levelOfLCA(u, v);
	return distAtLevel(level);
}

int levelOfLCA(int u, int v) {
	return levelOfLCA(leaves[u], leaves[v]);
}

int levelOfLCA(Node_t* u, Node_t* v) {
	if (u==NULL || v==NULL)
		return -1;
	
	while (u!=NULL && v!=NULL && u->lev!=v->lev) {
		if (u->lev > v->lev) {
			u = u->far;
		} else {
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		return -1;
	
	return u->lev;
}

pair<Node_t*,int> getLCA(int u, int v) {
	return getLCA(leaves[u], leaves[v]);
}

pair<Node_t*,int> getLCA(Node_t* u, Node_t* v) {
	if (u==NULL || v==NULL)
            return make_pair(rt, -1);
	
	while (u!=NULL && v!=NULL && u!=v) {
		if (u->lev > v->lev) {
			u = u->far;
		} else {
			v = v->far;
		}
	}
	
	if (u==NULL || v==NULL)
		return make_pair(rt, -1);
	
	return make_pair(u, u->lev);
}

void calcDmaxPrune() {
	dmax = 0;
	for (int j=1; j<nV; ++j) {
		dmax = max(dmax, dist(V, pi[0], j));
	}
}

void constructHST_fast_dp(bool load, clock_t startClock) {
	// perm[i][j]: center of node i at level j
	// reverse LEVEL in the tree as 

	if (!load)
		randomization();
	calcDmaxPrune();
	
	// initialization
	H = ceil(log2(dmax+EPS)) + 1;
	double radius = pow2(H)*beta;
	int** perm = new int*[nV];
	for (int i=0; i<nV; ++i) {
		perm[i] = new int[H+2];
	}
	
	// using DP to construct the HST
	for (int i=0; i<nV; i++) {
		perm[i][1] = 0;
		for (int j=2; j<=H+1; ++j){	
			perm[i][j] = reverse_pi[i];		
		}
		for (int j=0; j<reverse_pi[i]; j++) {	
			double curd = dist(V, i, pi[j]);
			int k = getLevel(H, curd);
			perm[i][k] = min(perm[i][k], j);
		}
		
		perm[i][H+1] = reverse_pi[i];
		for (int k=H; k>=1; --k){ 
			perm[i][k] = min(perm[i][k], perm[i][k+1]);
		}
	}

	// label the far[][] with node ID
	int nid = 1;
	int *cnt = NULL, *p = NULL, *q = NULL;
	int *far = NULL, *far_ = NULL;
	Node_t** nodes = NULL;
	Node_t** nodes_ = NULL; 
	
	// cnt[i]: count
	// p: x
	// q: result
	cnt = new int[nV];
	far = new int[nV];
	far_ = new int[nV];
	p = new int[nV];
	q = new int[nV];
	nodes = new Node_t*[nV];
	nodes_ = new Node_t*[nV];
	
	// construct the root
	rt = new Node_t(0, pi[0], 1, NULL, 0);
	for (int i=0; i<nV; ++i) {
        nodes_[i] = rt;
		p[i] = i;
		far_[i] = 0;
	}
	
	for (int k=2; k<=H+1; ++k) {
		radius /= alpha;
		int bid = nid, i = 0, bi;
		memset(q, -1, sizeof(int)*nV);
		while (i < nV) {
			q[perm[p[i]][k]] = far[p[i]] = nid++;
			bi = i++;
			while (i<nV && far_[p[i]]==far_[p[i-1]]) {
				int pid = perm[p[i]][k];
				if (q[pid] == -1) {
					q[pid] = nid++;
				}
				far[p[i]] = q[pid];
				++i;
			}
			while (bi < i) {
				int pid = perm[p[bi]][k];
				q[pid] = -1;
				++bi;
			}
		}
		memset(cnt, 0, sizeof(int)*nV);
		for (int i=0; i<nV; ++i) {
			++cnt[far[i]-bid];
		}
		for (int j=1; j<nid-bid; ++j) {
			cnt[j] += cnt[j-1];
		}
		for (int i=nV-1; i>=0; --i) {
			int j = far[p[i]] - bid;
			q[--cnt[j]] = p[i];
		}
		// create the new node at the $k$-th level
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			nodes[q[i]] = new Node_t(far[q[i]], pi[perm[q[i]][k]], k, nodes_[q[i]], radius);
			addChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				nodes[q[j]] = nodes[q[i]];
				++j;
			}
		}
		// merge the new node with its parent 
		for (int i=0,j=0; i<nV; i=j) {
			j = i;
			mergeChild(nodes_[q[i]], nodes[q[i]]);
			while (j<nV && far[q[j]]==far[q[i]]) {
				++j;
			}
		}		
		if (k == H+1) {
			for (int j=0; j<nV; ++j) {
				leaves[j] = nodes[j];
			}
		}		
		swap(q, p);
		swap(far, far_);
		swap(nodes, nodes_);

#ifdef USE_TIMER
		if (k%2 == 1) {
			double usedTime = (clock() - startClock) / CLOCKS_PER_SEC;
			if (usedTime > timeLimit) {
				printf("Finished the %d-th level, %d left.\n", i, H+1-i);								
				printf("Time Limitation Exceed\n");
				break; 
			}
		}
#endif
	}
 
	usedMemory = 0;
	usedMemory += nid * sizeof(Node_t);
	printf("nV=%d, nid=%d, H=%d\n", nV, nid, H);
	usedMemory += nV*(H+2)*sizeof(int);
	
	delete[] cnt;
	delete[] far;
	delete[] far_;
	delete[] p;
	delete[] q;
	delete[] nodes;
	delete[] nodes_;
	for (int i=0; i<nV; ++i)
		delete[] perm[i];
	delete[] perm;
}
