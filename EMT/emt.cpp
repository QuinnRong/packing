#include <cmath>
#include <fstream>

using namespace std;

double Maxwell(double k1, double k2, double mu)
{	/*
	k1: matrix
	k2: filler
	mu: volume fraction of filler
	*/
	return k1*(2*k1+k2-2*(k1-k2)*mu)/(2*k1+k2+(k1-k2)*mu);
}

double Reconstruction(double k1, double k2, double mu)
{	/*
	k1: matrix
	k2: filler
	mu: volume fraction of filler
	*/
	double temp = (3*mu -1)*k2+(3*(1-mu)-1)*k1;
	double result = temp + sqrt(temp*temp+8*k1*k2);
	return result/4;
}

double eval(double t, double x, double f)
{
	double temp = (1 - t) * (1 - x) / (f - x);
	return (temp*temp*temp*f - 1);
}

double Bruggeman(double k1, double k2, double mu)
{	/*
	k1: matrix
	k2: filler
	mu: volume fraction of filler
	solve:
	1 = (1-t)^3((1-x)/(f-x))^3f
	*/
	double x = k2 / k1;
	double t = mu;
	double min = k1;
	double max = k2;
	double f = (min + max) / 2;
	double err = eval(t, x, f);
	while (abs(err) > 0.001)
	{
		if (err > 0)
			max = f;
		else
			min = f;
		f = (min + max) / 2;
		err = eval(t, x, f);
	}
	return f;
}

void dump()
{
	double k1 = 1, k2 = 10;
	ofstream out1("./output/EMT_Maxwell.txt");
	ofstream out2("./output/EMT_Reconstruction.txt");
	ofstream out3("./output/EMT_Bruggeman.txt");
	for (int i = 0; i < 100; ++i)
	{
		double mu = i * 1.0 / 100 * 0.5 + 0.1;
		out1 << mu << " " << Maxwell(k1, k2, mu) << endl;
		out2 << mu << " " << Reconstruction(k1, k2, mu) << endl;
		out3 << mu << " " << Bruggeman(k1, k2, mu) << endl;
	}
	out1.close();
	out2.close();
	out3.close();
}

int main(int argc, char *argv[])
{
	dump();
}