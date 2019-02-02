#include <fstream>
#include <iostream>
#include <unistd.h>		/* access */
#include <cstdlib>		/* system */

#include "discretize.h"

/*****
class GridField
*****/

GridField::GridField(int num, double radius)
{
	ngrid = 1 / radius;
	ncell = ngrid * ngrid * ngrid;
	cells = new int[ncell];
	for (int i = 0; i < ncell; ++i)
		cells[i] = -1;
	nsphere = num;
	binlist = new int[nsphere + 1];
	for (int i = 0; i <= nsphere; ++i)
		binlist[i] = -1;
}

GridField::~GridField()
{
	delete[] cells;
	delete[] binlist;
}

int GridField::get_cell(double x, double y, double z)
{
	int a = int(x * ngrid) % ngrid;
	int b = int(y * ngrid) % ngrid;
	int c = int(z * ngrid) % ngrid;
	return (a + b * ngrid + c * ngrid * ngrid);
}

void GridField::add_sphere(int id, double x, double y, double z)
{
	int cellID = get_cell(x, y, z);
	int prev = cells[cellID];
	if (prev != -1)
		binlist[id] = prev;
	cells[cellID] = id;
}

void GridField::print_field()
{
	int tot = 0;
	for (int i = 0; i < ncell; ++i)
	{
		if (cells[i] != -1)
		{
			std::cout << i << ": " << cells[i];
			++tot;
			int id = cells[i];
			while (binlist[id] != -1)
			{
				std::cout << " -> " << binlist[id];
				++tot;
				id = binlist[id];
			}
			std::cout << std::endl;
		}
	}
	std::cout << "tot = " << tot << std::endl;
}

void GridField::print_data()
{
	std::cout << "cells:" << std::endl;
	for (int i = 0; i < ncell; ++i)
	{
		std::cout << i << ": " << cells[i] << std::endl;
	}
	std::cout << "binlist:" << std::endl;
	for (int i = 1; i <= nsphere; ++i)
	{
		std::cout << i << ": " << binlist[i] << std::endl;
	}
}

/*****
class Box
*****/

Box::Box(const std::string& filename, int n, double r)
: num(n), radius(r), field(n, r)
{
	coords = new double*[num + 1];
	for (int i = 0; i <= num; ++i)
		coords[i] = new double[3];

	get_field(filename);
}

Box::~Box()
{
	delete[] coords;
}

void Box::get_field(const std::string& filename)
{
	std::ifstream ifs(filename);
	std::string line;
	for (int i = 0; i < skip; ++i)
		std::getline(ifs, line);
	int id, type;
	double x, y, z;
	for (int i = 0; i < num; ++i)
	{
		ifs >> id >> type >> x >> y >> z;
		coords[id][0] = x;
		coords[id][1] = y;
		coords[id][2] = z;
		field.add_sphere(id, x, y, z);
	}
	ifs.close();
}

void SaveFile(const std::string& file, int* type, int dimX, int dimY, int dimZ)
{
	FILE *out = fopen(file.c_str(), "w");

	fprintf(out, "\n%d atoms\n\n", dimX*dimY*dimZ);
	fprintf(out, "2 atom types\n\n");
	fprintf(out, "0 %d xlo xhi\n", dimX);
	fprintf(out, "0 %d ylo yhi\n", dimY);
	fprintf(out, "0 %d zlo zhi\n\n", dimZ);
	fprintf(out, "Masses\n\n");
	fprintf(out, "1 1\n");
	fprintf(out, "2 2\n\n");
	fprintf(out, "Atoms\n\n");

	int atomID = 0;
	for (int z = 0; z < dimZ; ++z)
	{
		for (int y = 0; y < dimY; ++y)
		{
			for (int x = 0; x < dimX; ++x)
			{
				fprintf(out, "%d %d %.1f %.1f %.1f\n", atomID + 1, type[atomID] + 1, x + 0.5, y + 0.5, z + 0.5);
				++atomID;
			}
		}
	}

	fclose(out);
}

double distance_square(double* curr, double* coor, double* pboffset)
{
	double res = 0;
	for (int i = 0; i < 3; ++i)
	{
		double dis = curr[i] - (coor[i] + pboffset[i]);
		res += dis * dis;
	}
	return res;
}

bool Box::in_sphere(double x, double y, double z)
{
    // now iterate through nearest neighbors
    double delta = 1.0 / field.ngrid;
    double curr[3] = {x, y, z};

    bool flag = false;
    double pboffset[3];        // nearest image offset
    int grid[3] = {-1, -1, -1};
    while (1)
    {
        for(int k = 0; k < 3; k++)
        { 
            if (curr[k] + grid[k] * delta < 0)              // out of bounds to left
                pboffset[k] = -1;
            else if (curr[k] + grid[k] * delta >= 1)   // out of bounds to right
                pboffset[k] = 1;
            else
                pboffset[k] = 0;
        }
        int cellID = field.get_cell(x + grid[0] * delta + 1, y + grid[1] * delta + 1, z + grid[2] * delta + 1);
        int j = field.cells[cellID];
        while (j != -1)
        {
        	if (distance_square(curr, coords[j], pboffset) <= radius*radius)
        	{
        		flag = true;
        		break;
        	}
            j = field.binlist[j];
        }
        if (flag)
        {
        	return true;
        }

        // A. Donev:     
        // This code makes this loop dimension-independent
        // It is basically a flattened-out loop nest of depth DIM
        // (-1, -1) (0, -1) (1, -1) (-1, 0) (0, 0) (1, 0) (-1, 1) (0, 1) (1, 1)
        int ii;
        for(ii = 0; ii < 3; ii++)
        {
            grid[ii] += 1;
            if(grid[ii] <= 1) break;
            grid[ii] = -1;
        }
        if(ii >= 3) break;
    }
    return false;
}

void Box::get_section(const std::string& filename, std::string dir, double dis, int res)
{
	std::string file = filename;
	file += dir + "_" + std::to_string(dis) + "_" + std::to_string(res) + ".txt";
	int* type = new int[res*res];

	if (dir == "x")
	{
		int atomID = 0;
		for (int z = 0; z < res; ++z)
		{
			for (int y = 0; y < res; ++y)
			{
				type[atomID] = in_sphere(dis, y * 1.0 / res, z * 1.0 / res) ? 1 : 0;
				++atomID;
			}
		}
		SaveFile(file, type, 1, res, res);
	}
	else if (dir == "y")
	{
		int atomID = 0;
		for (int z = 0; z < res; ++z)
		{
			for (int x = 0; x < res; ++x)
			{
				type[atomID] = in_sphere(x * 1.0 / res, dis, z * 1.0 / res) ? 1 : 0;
				++atomID;
			}
		}
		SaveFile(file, type, res, 1, res);
	}
	else if (dir == "z")
	{
		int atomID = 0;
		for (int y = 0; y < res; ++y)
		{
			for (int x = 0; x < res; ++x)
			{
				type[atomID] = in_sphere(x * 1.0 / res, y * 1.0 / res, dis) ? 1 : 0;
				++atomID;
			}
		}
		SaveFile(file, type, res, res, 1);
	}

	delete[] type;
}

double Box::get_section_pf(std::string dir, double dis, int res)
{
	int sum = 0;

	if (dir == "x")
	{
		for (int z = 0; z < res; ++z)
		{
			for (int y = 0; y < res; ++y)
			{
				sum += in_sphere(dis, y * 1.0 / res, z * 1.0 / res) ? 1 : 0;
			}
		}
	}
	else if (dir == "y")
	{
		for (int z = 0; z < res; ++z)
		{
			for (int x = 0; x < res; ++x)
			{
				sum += in_sphere(x * 1.0 / res, dis, z * 1.0 / res) ? 1 : 0;
			}
		}
	}
	else if (dir == "z")
	{
		for (int y = 0; y < res; ++y)
		{
			for (int x = 0; x < res; ++x)
			{
				sum += in_sphere(x * 1.0 / res, y * 1.0 / res, dis) ? 1 : 0;
			}
		}
	}

	return sum * 1.0 / res / res;
}

void Box::get_structure(const std::string& filename, int res)
{
	std::string file = filename;
	int* type = new int[res*res*res];

	int atomID = 0;
	for (int z = 0; z < res; ++z)
	{
		for (int y = 0; y < res; ++y)
		{
			for (int x = 0; x < res; ++x)
			{
				type[atomID] = in_sphere(x * 1.0 / res, y * 1.0 / res, z * 1.0 / res) ? 1 : 0;
				++atomID;
			}
		}
	}
	SaveFile(file, type, res, res, res);

	delete[] type;
}

void make_directory(const std::string &path)
{
    // path not exist
    if (access(path.c_str(), 0) == -1)
    {
        std::string cmd = "mkdir " + path;
        system(cmd.c_str());
    }
}

void delete_directory(const std::string &path)
{
    // path exist
    if (access(path.c_str(), 0) == 0)
    {
        std::string cmd = "rm -rf " + path;
        system(cmd.c_str());
    }
}

void Box::get_fenics_input(const std::string& root, int idx, int res)
{
    // mkdir
    std::string path = root + std::to_string(idx);
    delete_directory(path);
    make_directory(path);
    // save to file
    for (int z = 0; z < res; ++z)
    {
        std::string filename = path + "/3D_" + std::to_string(idx) + "_" + std::to_string(z) + ".dat";
        std::ofstream out(filename);
        for (int y = 0; y < res; ++y)
        {
            for (int x = 0; x < res; ++x)
            {
                out << (in_sphere(x * 1.0 / res, y * 1.0 / res, z * 1.0 / res) ? 1 : 0) << " ";
            }
            out << std::endl;
        }
        out.close();
    }
}

void Box::get_fenics_input_single_file(const std::string& root, int idx, int res)
{
    // save to file
    std::string filename = root + "3D_" + std::to_string(idx) + ".dat";
    std::ofstream out(filename);
    for (int z = 0; z < res; ++z)
    {
        for (int y = 0; y < res; ++y)
        {
            for (int x = 0; x < res; ++x)
            {
                out << (in_sphere(x * 1.0 / res, y * 1.0 / res, z * 1.0 / res) ? 1 : 0) << " ";
            }
            out << std::endl;
        }
    }
    out.close();
}