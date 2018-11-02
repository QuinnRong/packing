#include <fstream>
#include <iostream>

#include "format_change.h"

GridField::GridField(double radius, int num)
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

Box::Box(double r, int n): radius(r), num(n), field(r, n)
{
	coords = new double*[num + 1];
	for (int i = 0; i <= num; ++i)
		coords[i] = new double[3];
}

Box::~Box()
{
	delete[] coords;
}

void Box::get_field(const std::string& filename, int skip)
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
	freopen(file.c_str(), "w", stdout);

	printf("\n%d atoms\n\n", dimX*dimY*dimZ);
	printf("2 atom types\n\n");
	printf("0 %d xlo xhi\n", dimX);
	printf("0 %d ylo yhi\n", dimY);
	printf("0 %d zlo zhi\n\n", dimZ);
	printf("Masses\n\n");
	printf("1 1\n");
	printf("2 2\n\n");
	printf("Atoms\n\n");

	int atomID = 0;
	for (int z = 0; z < dimZ; ++z)
	{
		for (int y = 0; y < dimY; ++y)
		{
			for (int x = 0; x < dimX; ++x)
			{
				printf("%d %d %d %d %d\n", atomID + 1, type[atomID] + 1, x, y, z);
				++atomID;
			}
		}
	}

	fclose(stdout);
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

void Box::get_section(std::string dir, double dis, int res)
{
	std::string file = dir;
	file += "_" + std::to_string(dis) + "_" + std::to_string(res) + ".txt";
	int* type = new int[res*res];

	if (dir == "y")
	{
		int atomID = 0;
		for (int z = 0; z < res; ++z)
		{
			for (int x = 0; x < res; ++x)
			{
				type[atomID] = in_sphere(x * 1.0 / res, 0, z * 1.0 / res) ? 1 : 0;
				++atomID;
			}
		}
		SaveFile(file, type, res, 1, res);
	}

	delete[] type;
}

void Box::get_structure(int res)
{
	std::string file = std::to_string(res) + ".txt";
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
		SaveFile(file, type, res, res, res);
	}

	delete[] type;
}
