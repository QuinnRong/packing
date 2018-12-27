#include <string>
#include <iostream>
#include <fstream>

#include "discretize.h"

using namespace std;

void run_once()
{
	Box box("../spheres/output/struct_0.dat", 100, 0.094413);
	box.get_structure("0.txt", 100);
	box.get_section("0.txt", "x", 0, 100);
	box.get_section("0.txt", "y", 0, 100);
	box.get_section("0.txt", "z", 0, 100);
}

int main()
{
	string root = "../..";
	string input_dir  = root + SPT + "spheres" + SPT + "run_2_train/output" + SPT;
	string output_dir = root + SPT + "utility" + SPT + "run_2_train/output" + SPT;
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

		// box.get_structure(output_dir + file, resolution);

		// box.get_fenics_input(output_dir, idx++, resolution);
		box.get_fenics_input_single_file(output_dir, idx++, resolution);

		// box.get_section(output_dir + file, "x", 0, resolution);
		// box.get_section(output_dir + file, "y", 0, resolution);
		// box.get_section(output_dir + file, "z", 0, resolution);
	}

	ifs.close();

	return 0;
}
