#include "GeneticAlgorithmExample.h"

double BitConvertTools::get_2pwr(int x)
{
	double ret = 1;
	for (int i = 0; i < x; i++)
	{
		ret *= 2;
	}
	return ret;
}

BitConvertTools::BitConvertTools(double _l, double _r, double _precision)
{
	l = _l;
	r = _r;
	precision = _precision;
	double delta = r - l;
	block = delta / precision;
	bit = ceil(log2(block));
	range = get_2pwr(bit) - 1;
}

void BitConvertTools::setData(double x, data& data)
{
	x -= l;
	int interval = (x / (r-l)) * range;
	int pos = 0;
	while (interval)
	{
		data[pos++] = interval % 2;
		interval /= 2;
	}
}

double BitConvertTools::restoreData(data& data)
{
	double ret = 0;
	int interval = 0;
	int base = 1;
	for (int i = 0; i < bit; i++)
	{
		interval += data[i] * base;
		base <<= 1;
	}
	ret = l + 1.0 * interval / range * (r - l);
	return ret;
}

int BitConvertTools::getBit()
{
	return bit;
}

void BitConvertTools::concat(data& srcA, data& srcB, BitConvertTools bctA, BitConvertTools bctB, data& dest)
{
	for (int i = 0; i < bctB.bit; i++)
	{
		dest[i] = srcB[i];
	}
	for (int i = 0; i < bctA.bit; i++)
	{
		dest[i + bctB.bit] = srcA[i];
	}
}




Pair BitConvertTools::split(data& src, BitConvertTools& bctA, BitConvertTools& bctB)
{
	data divider;
	data retA, retB;
	for (int i = 0; i < bctB.bit; i++)
	{
		divider[i] = 1;
	}
	retB = src & divider;
	divider.reset();
	for (int i = 0; i < bctA.bit; i++)
	{
		divider[i] = 1;
	}
	retA = (src >> bctB.bit) & divider;
	return Pair(retA, retB);
}

void BitConvertTools::CrossOver(data& srcA, data& srcB, int bit, int pos)
{
	data retA, retB;
	data divider;
	for (int i = 0; i < pos; i++)
	{
		divider[i] = 1;
	}
	retA |= srcB & divider;
	retB |= srcA & divider;
	divider.reset();
	for (int i = pos; i < bit; i++)
	{
		divider[i] = 1;
	}
	retA |= srcA & divider;
	retB |= srcB & divider;
}






GeneticAlgo::GeneticAlgo(Func* _F, BitConvertTools& _a, BitConvertTools& _b, int _N, double _pc, double _pm, int _Gmax)
{
	a = _a;
	b = _b;
	N = _N;
	pc = _pc;
	pm = _pm;
	Gmax = _Gmax;
	F = _F;
	bit = a.getBit() + b.getBit();
	srand(time(0));
	Colony = new vector<data>;
}

double GeneticAlgo::_eval(data& data)
{
	Pair x = BitConvertTools::split(data, a, b);
	double x1 = a.restoreData(x.first);
	double x2 = b.restoreData(x.second);
	return (*F)(x1, x2);
}

void GeneticAlgo::generateColony()
{
	for (int i = 0; i < N; i++)
	{
		data to_ins;
		for (int j = 0; j < bit; j++)
		{
			to_ins[j] = rand() % 2;
		}
		Colony->push_back(to_ins);
	}
}


double x1, x2;
double GeneticAlgo::_step()
{
	double best = DBL_MIN;
	double tot = 0;
	double* choice_prob = new double[N];
	double* ovl_prob = new double[N];
	double* val = new double[N];
	vector <data> *newColony = new vector <data>;
	for (int i = 0; i < N; i++)
	{
		val[i] = _eval((*Colony)[i]);
		tot += val[i];
		best = std::max(best, val[i]);
	}
	for (int i = 0; i < N; i++)
	{
		choice_prob[i] = val[i] / tot;
		ovl_prob[i] = (i == 0 ? choice_prob[i] : ovl_prob[i - 1] + choice_prob[i]);
	}
	//Spin
	for (int i = 0; i < N; i++)
	{
		double prob = 0.9999 * rand() / RAND_MAX;
		int pos = std::upper_bound(ovl_prob, ovl_prob + N, prob) - ovl_prob;
		newColony->push_back((*Colony)[pos]);
	}
	//Crossover
	vector <int> to_crossover;
	for (int i = 0; i < N; i++)
	{
		double prob = 1.0 * rand() / RAND_MAX;
		if (prob < pc)
		{
			to_crossover.push_back(i);
		}
	}
	if (to_crossover.size() & 1)
	{
		to_crossover.pop_back();
	}
	int size = to_crossover.size();
	for (int i = 0; i < size; i++)
	{
		int pos = (1.0 * rand() / RAND_MAX) * bit;
		data& dA = (*newColony)[to_crossover[i]];
		data& dB = (*newColony)[to_crossover[size - i - 1]];
		BitConvertTools::CrossOver(dA, dB, bit, pos);
	}
	//Metamorphosis
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < bit; j++)
		{
			double prob = 1.0 * rand() / RAND_MAX;
			if (prob < pm)
			{
				(*newColony)[i][j] = !(*newColony)[i][j];
			}
		}
	}
	delete Colony;
	Colony = newColony;
	delete[] choice_prob;
	delete[] ovl_prob;
	delete[] val;
	return best;
}

double GeneticAlgo::Solve()
{
	generateColony();
	double ans = DBL_MIN;
	for (int Gen = 0; Gen < Gmax; Gen++)
	{
		ans = std::max(ans, _step());
	}
	return ans;
}

GeneticAlgo::~GeneticAlgo()
{
	delete Colony;
}


double f(double x1, double x2) //target function
{
	return 21.5+x1*sin(4*pi*x1)+x2*sin(20*pi*x2);
}
Func F = &f;
int main()
{
	std::fstream f("performance.txt",std::fstream::app);
	//This example solves a simple maximum problem, the target function lies below
	GeneticAlgo* GA;
	BitConvertTools a(-3.0, 12.1, 0.0001);
	BitConvertTools b(4.1, 5.8, 0.0001);
	for (double pc = 0.1; pc <= 0.9; pc += 0.1)
	{
		for (double pm = 0.01; pm <= 0.09; pm += 0.01)
		{
			double ans = 0;
			for (int test = 0; test < 30; test++)
			{
				GA = new GeneticAlgo(&F, a, b, 100, pc, pm, 2000);
				ans += GA->Solve();
				delete GA;
			}
			f <<std::setprecision(10) << pc << " " << pm << " " << ans / 30 << std::endl;
			printf("%.1lf %.2lf %.10lf\n", pc, pm, ans / 30);
		}
	}
	return 0;
}