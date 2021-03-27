/**
	\Author: Trasier
	\Date:	2019/10/9
**/

#include "HST.h"
#include "global.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

const int MAX_HEIGHT = 40;
int H = 0;
int *pi = NULL;
int *reverse_pi = NULL;
double dmax = -1.0;
double gamma_or = -1.0;
double* expks = NULL;
double* sumks = NULL;
double alpha = 2.0;
double logAlpha = 1.0;
Node_t* rt = NULL;
Node_t** leaves = NULL;

inline double log2(double x) {
	return log(x) / logAlpha;
}

inline double pow2(int i) {
	return expks[i];
}

inline int getLevel(double H, double dist) {
	if (dist < 1.0) return H+1;
	
	int k = ceil(log2(dist/gamma_or));
	
	if (expks[k]*gamma_or == dist)
		++k;
	
	return H+1-k;
}

void addChild(Node_t* far, Node_t* u) {
	u->far = far;
	u->cid = far->child.size();
	far->child.push_back(u);
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
	sumks[0] = 0.0;
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
    fin >> gamma_or;
    fin >> alpha;
    string a;
    // getline(fin,a);
	initMemory_HST(nV);
    for(int i=0;i<nV;i++){
        int tmp = i;
        fin >> tmp;
        pi[i] = tmp;
        reverse_pi[tmp] = i;
    }

    for (int i=0; i<nV; ++i) {
		for (int j=0; j<MAX_DIM; ++j)
			fin >> V[i].x[j];
    }
    fin.close();
}

void calcDmax() {
	if (dmax >= 0.0) return ;
	
	// initialize the parameters for decomposition
	dmax = 0.0;
	for (int i=0; i<nV; i++) {
		for (int j=i+1; j<nV; j++) {
			dmax = max(dmax, dist(V, i, j));
		}
	}
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
	gamma_or = (rand()%((int)alpha) + 1) * 1.0	/ alpha;
}

// donot forget to revise freeMemory_HST() to delete the whole tree
void constructHST_slow(bool load, clock_t startClock) {
	if (!load)
		randomization();
	calcDmax();
	
	// initialization
	H = ceil(log2(dmax+EPS)) + 1;
	double radius = pow2(H)*gamma_or;
	
	/*** revise the following lines***/	
	
	int* visit = NULL;
	Node_t** nodes = NULL;
	Node_t** nodes_ = NULL; 
	visit = new int[nV];
	nodes = new Node_t*[nV];
	nodes_ = new Node_t*[nV];
	rt = new Node_t(0, pi[0], 1, NULL, 0);
	for (int i=0; i<nV; ++i) {
        nodes_[i] = rt;
	}
	/*** revise the above lines***/
	
	vector<vector<int> > preC, curC;
	vector<int> vtmp;
	for (int i=0; i<nV; ++i)
		vtmp.push_back(i);
	preC.push_back(vtmp);
	vtmp.clear();
	
	// construct the HST by brute-force
	int nid = 1;
	
	for (int i=2; i<=H+1; ++i) {
		radius /= alpha;

		memset(visit, 0, sizeof(int)*nV);
		for (int cid=0; cid<preC.size(); ++cid) {
			vector<int> cluster = preC[cid];
			for (int j=0; j<nV; ++j) {
				vector<int> newCluster;
				for (int uid=cluster.size()-1; uid>=0; --uid) {
					int u = cluster[uid];
					if (!visit[u] && dist(V, u, pi[j])<radius) {
						newCluster.push_back(u);
						visit[u] = true;
					}
				}
				if (!newCluster.empty()) {
					curC.push_back(newCluster);
					/*** revise the following lines***/
					Node_t* par = nodes_[newCluster[0]];
					Node_t* child = new Node_t(nid, pi[j], i, par, radius);
					///// add child
					addChild(par, child);
					for (int vid=newCluster.size()-1; vid>=0; --vid) {
						nodes[newCluster[vid]] = child;
					}
					/*** revise the above lines***/
					++nid;
				}
			}
		}
		/*** revise the following lines***/
		if (i == H+1) {
			for (int j=0; j<nV; ++j) {
				leaves[j] = nodes[j];
			}
		}
		swap(nodes, nodes_);

#ifdef WATCH_MEM
	watchSolutionOnce(getpid(), usedMemory);
#endif
		/*** revise the above lines***/

		preC = curC;
		curC.clear();
#ifdef USE_TIMER
		if (i%2 == 1) {
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
	
	delete[] visit;
	delete[] nodes;
	delete[] nodes_;
}

double distAtLevel(int lev) {
	if (lev >= H+1) return 0.0;
	return sumks[H-lev] * gamma_or * 2.0;

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
