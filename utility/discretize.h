#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <string>

#define SPT "\\"		// for windows
// #define SPT "/"		// for linux

class GridField
{
public:
	GridField(int num, double radius);
	~GridField();

	int get_cell(double x, double y, double z);
	void add_sphere(int id, double x, double y, double z);
	void print_field();
	void print_data();

	int ngrid, ncell, nsphere;
	int* cells;
	int* binlist;
};

class Box
{
public:
	Box(const std::string& filename, int n, double r);
	~Box();

	void print_field() { field.print_field(); };
	void print_data() { field.print_data(); };

	void get_section(const std::string& filename, std::string dir, double dis, int res);
	void get_structure(const std::string& filename, int res);

	void get_fenics_input(const std::string& path, int idx, int res);

private:
	void get_field(const std::string& filename);
	bool in_sphere(double x, double y, double z);

	int num;
	double radius;
	int skip = 15;		// # of lines before coords in struct file
	double** coords;
	GridField field;
};

#endif