#include <iostream>
#include <sstream>
#include <chrono>
#include <thread>
#include "Duct.h"

using namespace std;

int main() 
{
	int longitud = 50;
	int radio    = 10;
	int n_gens   = 50;
	int n_hilos  = 4;

	chrono::time_point<chrono::steady_clock> tic, toc;
	
	tic = chrono::steady_clock::now();
	Duct duct1(longitud, radio, true, n_hilos);
	duct1.execute(n_gens);
	toc = chrono::steady_clock::now();
	double elapsed = chrono::duration_cast<chrono::nanoseconds>(toc-tic).count();
	cout << elapsed * 1e-9 << endl;
	return 0;

	/*
	vector<thread> v;
	for(int i=0 ; i<6 ; i++)
		v.push_back(thread (foo));
	for(int i=0 ; i<6 ; i++)
		v[i].join();

	cout << "sacabo\n";
	return 0;
	*/
}