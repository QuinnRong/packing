#include <string>
#include <iostream>

#include "format_change.h"

using namespace std;

int main()
{
	double resolution = 100;
	
	double radius = 0.106;
	int num = 100;
	string file = "../spheres/output/struct.dat";
	int skip = 15;

	Box box(radius, num);
	box.get_field(file, skip);
	// box.print_field();
	box.get_structure(resolution);

	// string direction = "y";
	// double distance = 0;
	// double resolution = 50;
	// box.get_section(direction, distance, resolution);
}