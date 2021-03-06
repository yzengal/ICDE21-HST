/**
	\Author: Trasier
	\Date:	2019/05/27
**/
#include "HST_dpo.h"

// #define LOCAL_DEBUG

Tree_t bs_tree;
int* mp;
int* rmp;
bool* deleted;
int loc_n;
int tot_n;
int insert_n =0;
double alpha;
vector<location_t> locs;
vector<location_t> clocs;
vector<int> cids;
pair<double,double> obj(1.0, 0.0);
clock_t start_time, end_time, last_time;
double totUsedTime = 0.0;

pair<double,double> getPartDistortion() {
	double d, dt, rat;
	double ret = 1.0+EPS;
	double ret_ = 0.0;
	double cnt = 0.0;
	bool checkFlag = true;
	
	for (int i=0; i<clocs.size(); ++i) {
		int u = cids[i];
		for (int j=0; j<locs.size()-clocs.size()+i; ++j) {
			d = dist(clocs[i], locs[j]);
			dt = bs_tree.distOnHST(u, j);
			if (d > 0) {
				rat = max(dt/d, 1.0);
				ret = max(ret, rat);
				ret_ += rat;
				cnt += 1.0;
			}
			
			if (dt < d) {
				checkFlag = false;
			}
		}
	}
	if (cnt <= 0) {
		ret_ = 1.0;
	} else {
		ret_ /= cnt;
	}
	
	// if (!checkFlag) {
		// puts("false");
	// }
	
	return make_pair(ret, ret_);	
}

void rebuild() {
	last_time = clock();
	bs_tree.freeMemory();
	bs_tree.initLocation(locs.size(), locs, alpha);
	bs_tree.constructHST(start_time);
	end_time = clock();
	
	pair<double,double> tmp = getPartDistortion();
	obj.first = max(obj.first, tmp.first);
	obj.second += tmp.second;
	
	double timeCost = (double)(end_time-last_time)/CLOCKS_PER_SEC;
	long long usedMemory_ = bs_tree.memoryCost();
	printf("reHSTO %.4lf %.4lf %.4lf %lld\n", tmp.first, tmp.second, timeCost, usedMemory_);
	fflush(stdout);
	totUsedTime += timeCost;
	insert_n += 1;
}

inline void insert_loc(int id, location_t& loc) {
	mp[id] = locs.size();
	rmp[locs.size()] = id;
	cids.push_back(locs.size());
	clocs.push_back( loc );
	locs.push_back( loc );
}

void delete_loc(int id) {
	int end = locs.size()-1, now = mp[id];
	int eid = rmp[end];
	
	// we can safely ignore ``now''
	locs[now] = locs[end];
	mp[eid] = now;
	rmp[now] = eid;
	locs.pop_back();
	deleted[id] = true;
}

void static_input() {
	char op[10];
	int num;
	int id;
	location_t loc;
	
	clocs.clear();
	scanf("%s %d", op,&num);
	for (int i = 0; i < num; ++i) {
		scanf("%d", &id);
		for (int j=0; j<MAX_DIM; ++j) {
			scanf("%lf", &loc.x[j]);
		}
		insert_loc(id, loc);
	}
	
	last_time = clock();
	bs_tree.initLocation(locs.size(), locs, alpha);
	bs_tree.constructHST(start_time);

	end_time = clock();
	pair<double,double> tmp = make_pair(bs_tree.distor, bs_tree.distor_);//bs_tree.getDistortion();
	obj.first = max(obj.first, tmp.first);
	obj.second += tmp.second;
	
	double timeCost = (double)(end_time-last_time)/CLOCKS_PER_SEC;
	long long usedMemory_ = bs_tree.memoryCost();
	printf("reHSTO %.4lf %.4lf %.4lf %lld\n", tmp.first, tmp.second, timeCost, usedMemory_);
	fflush(stdout);
	totUsedTime += timeCost;
	insert_n += 1;
}

void dec_or_inc() {
	char op[10];
	int num;
	location_t loc;
	int id;
	
	scanf("%s %d",op,&num);
	if (op[0] == 'I') {
		clocs.clear();
		cids.clear();
		for (int i = 0; i < num; ++i) {
			scanf("%d", &id);
			for (int j=0; j<MAX_DIM; ++j) {
				scanf("%lf", &loc.x[j]);
			}
			insert_loc(id, loc);
		}
		rebuild();
	}
	
	if (op[0] == 'D') {
		for (int i = 0; i < num; ++i) {
			scanf("%d", &id);
			delete_loc(id);
		}
	}
}

int main(int argc, char **argv) {
	string srcFileName;
	string desFileName;
	int nop;
	
	if (argc > 1)
		srcFileName = string(argv[1]);
	else
		srcFileName = "data_0.txt";
	if (argc > 2)
		desFileName = string(argv[2]);
	freopen(srcFileName.c_str(), "r", stdin);
	
	start_time = clock();
	scanf("%d %d %lf", &tot_n, &nop, &alpha);
	mp = new int[tot_n+5];
	rmp = new int[tot_n+5];
	deleted = new bool[tot_n+5];
	memset(deleted, false, sizeof(bool)*(tot_n)+5);
	
	static_input();
	
	/*
	* rebuild can only happen in the insert part
	* So I remove it to inside the dec_or_inc
	*/

	// rebuild(alpha,time_limit);
	for (int i = 2; i <= nop; ++i) { // the first op is static input...
		dec_or_inc();
		// rebuild(alpha,time_limit);
	}
	
#ifdef WATCH_MEM
    watchSolutionOnce(getpid(), usedMemory);
#endif
	printf("reHSTO %.4lf %.4lf %.4lf %lld\n", obj.first, obj.second/insert_n, totUsedTime, usedMemory);
	fclose(stdout);
	
	bs_tree.freeMemory();
	delete[] mp;
	delete[] rmp;
	delete[] deleted;
	
	return 0;
}
