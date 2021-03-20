#pragma once

#include <iostream>
#include <cstdio>
#include <bitset>
#include <cmath>
#include <vector>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <iomanip>
const double pi = acos(-1.0);


const int maxbit = 50;

using std::vector;
typedef std::bitset<maxbit> data;
typedef std::pair<data, data> Pair;
typedef double (*Func)(double, double);
typedef long long ll;

class BitConvertTools
{
protected:
	double l, r;
	int bit, block;
	double range, precision;
	double get_2pwr(int x);
public:
	BitConvertTools(){}
	BitConvertTools(double _l, double _r, double _precision);
	void setData(double x, data&);
	double restoreData(data&);
	int getBit();
	static void concat(data& srcA, data& srcB, BitConvertTools bctA, BitConvertTools bctB, data& dest);
	static Pair split(data&, BitConvertTools&, BitConvertTools&);
	static void CrossOver(data& srcA, data& srcB, int bit, int pos);
};

class GeneticAlgo
{
protected:
	BitConvertTools a, b;
	int bit;
	int N, Gmax;
	double pc, pm;
	Func* F;
	vector<data> *Colony;
public:
	GeneticAlgo(Func*, BitConvertTools& _a, BitConvertTools& _b,
		int _N, double _pc, double _pm, int _Gmax);
	double _eval(data&);
	void generateColony();
	double _step();
	double Solve();
	~GeneticAlgo();
};
