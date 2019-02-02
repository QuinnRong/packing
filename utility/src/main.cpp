#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

#include "discretize.h"

using namespace std;

void run_once()
{
	int num = 125;
	double radius = 0.08743590505618876;
	int resolution = 100;
	int idx = 3;
	double dist[3] = {0.0, 0.0, 0.0};

	string input_dir = "../spheres/output/";
	string output_dir = "./output/";
	string ids = to_string(idx);
	string input_file = input_dir + "struct_" + ids + ".dat";
	Box box(input_file, num, radius);
	box.get_fenics_input_single_file(output_dir, idx, resolution);
	box.get_structure(output_dir + ids + ".txt", resolution);
	box.get_section(output_dir + ids, "x", dist[0], resolution);
	box.get_section(output_dir + ids, "y", dist[1], resolution);
	box.get_section(output_dir + ids, "z", dist[2], resolution);
}

void run_dataset()
{
	string root = "..";
	string input_dir  = root + SPT + "spheres" + SPT + "output" + SPT;
	string output_dir = root + SPT + "utility" + SPT + "output" + SPT;
	make_directory(output_dir);
	int resolution    = 100;

	ifstream ifs(input_dir + "summary.txt");

	string line; getline(ifs, line);
	string file; int num; double radius, pf;
	int idx = 0;
	while (ifs >> file >> num >> radius >> pf)
	{
		cout << file << " " << num << " " << radius << " " << pf << endl;
		Box box(input_dir + file, num, radius);

		box.get_structure(output_dir + file, resolution);

		// box.get_fenics_input(output_dir, idx++, resolution);
		// box.get_fenics_input_single_file(output_dir, idx++, resolution);

		// box.get_section(output_dir + file, "x", 0, resolution);
		// box.get_section(output_dir + file, "y", 0, resolution);
		// box.get_section(output_dir + file, "z", 0, resolution);
	}

	ifs.close();
}

void get_std()
{
	string root = "..";
	string input_dir  = root + SPT + "spheres" + SPT + "output" + SPT;
	string output_dir = root + SPT + "utility" + SPT + "output" + SPT;
	make_directory(output_dir);
	int resolution    = 100;

	ifstream ifs(input_dir + "summary.txt");

	string line; getline(ifs, line);
	string file; int num; double radius, pf;
	int idx = 0;
	double std_all = 0;
	int iter = resolution / 4;
	while (ifs >> file >> num >> radius >> pf)
	{
		cout << file << " " << num << " " << radius << " " << pf << endl;
		Box box(input_dir + file, num, radius);

		double std = 0;
		for (int i = 0; i < iter; ++i)
		{
			double section_pf = 0;
			section_pf += box.get_section_pf("x", i * 1.0 / resolution, resolution);
			section_pf += box.get_section_pf("y", i * 1.0 / resolution, resolution);
			section_pf += box.get_section_pf("x", i * 1.0 / resolution + 0.25, resolution);
			section_pf += box.get_section_pf("y", i * 1.0 / resolution + 0.25, resolution);
			section_pf += box.get_section_pf("x", i * 1.0 / resolution + 0.5, resolution);
			section_pf += box.get_section_pf("y", i * 1.0 / resolution + 0.5, resolution);
			section_pf += box.get_section_pf("x", i * 1.0 / resolution + 0.75, resolution);
			section_pf += box.get_section_pf("y", i * 1.0 / resolution + 0.75, resolution);
			section_pf /= 8; // cout << section_pf << " ";
			std += (pf - section_pf) * (pf - section_pf);
		}
		std_all += std;
		++idx;
		cout << sqrt(std / iter) / pf<< endl;
	}
	cout << endl << sqrt(std_all / iter / idx) / pf << endl;

	ifs.close();
}

int main()
{
	run_dataset();
	run_once();
	// get_std();

	return 0;
}
