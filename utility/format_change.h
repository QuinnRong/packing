#ifndef FORMAT_CHANGE_H
#define FORMAT_CHANGE_H

#include <string>

class GridField
{
public:
	GridField(double radius, int num);
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
	Box(double r, int n);
	Box(const std::string& filename, int n, int r);
	~Box();

	void get_field(const std::string& filename, int skip);
	void print_field() { field.print_field(); };
	void print_data() { field.print_data(); };

	void get_section(std::string dir, double dis, int res);
	void get_structure(int res);

private:
	bool in_sphere(double x, double y, double z);

	double radius;
	int num;
	double** coords;
	GridField field;
};

#endif