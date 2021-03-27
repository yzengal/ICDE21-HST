/**
	\Author: Trasier
	\Date:	2019/10/9
**/
#include <bits/stdc++.h>
using namespace std;

#include "global.h"
#include "HST_opt.h"
#ifdef WATCH_MEM
#include "monitor.h"
#endif

int main(int argc, char **argv) {
	string srcFileName;
	string desFileName;
	
	if (argc > 1)
		srcFileName = string(argv[1]);
	if (argc > 2)
		desFileName = string(argv[2]);
	
    clock_t start_time, end_time;
	
    initLocation(srcFileName);

    start_time = clock();
	constructHST_fast_opt(true, start_time);
    end_time = clock();

    clock_t lastTime = end_time-start_time;
	printf("OPT %.4lf %lld\n", (double)lastTime/CLOCKS_PER_SEC, usedMemory);
    
	
	freeMemory_HST();
	fclose(stdout);
	
	return 0;
}
