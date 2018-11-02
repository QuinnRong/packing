//===========================================================
//===========================================================
//===========================================================
//
//  Molecular dynamics simulation of hardspheres
//
//===========================================================
//===========================================================
//===========================================================

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>

#include "box.h"
#include "sphere.h"
#include "event.h"
#include "heap.h"
#include "read_input.h"


int main(int argc, char **argv)
{
    read_input input;
    int error = input.read(argc, argv);
    if (error) return error;

    std::vector<double> pf = {0.3, 0.4, 0.5, 0.6};

    for (int i = 0; i < pf.size(); ++i)
    {
        input.maxpf = pf[i];
        std::string writefile = "./output/struct_" + std::to_string(i) + ".dat";
        strcpy(input.writefile, writefile.c_str());
        std::string datafile =  "./output/status_" + std::to_string(i) + ".dat";
        strcpy(input.datafile, datafile.c_str());

        double r = pow(input.initialpf*pow(SIZE, DIM)/(input.N*VOLUMESPHERE), 1.0/((double)(DIM)));

        box b(input.N, r, input.growthrate, input.maxpf);

        std::cout << "ngrids = " << b.ngrids << std::endl;
        std::cout << "DIM = " << DIM << std::endl;

        std::cout << "Creating new positions of spheres" << std::endl;
        b.CreateSpheres(input.temp);

        std::ofstream output(input.datafile);
        output.precision(16);  

        output << "step packing-fraction pressure energy-change total-events" << std::endl;
        int step = 0;
        while ((b.pf < input.maxpf) && (b.pressure < input.maxpressure)) 
        {
            std::cout << "step = " << step << " pf = " << b.pf << " pressure = " << b.pressure << std::endl;
            b.Process(input.eventspercycle*input.N);
            output << step++ << " " << b.pf << " " << b.pressure << " "
                   << b.energychange << " " << b.neventstot << " " << std::endl;
            b.Synchronize(true);
        }
        std::cout << "\nfinal radius: " << b.r << std::endl;
        output << "\nfinal radius: " << b.r << std::endl;

        output.close();
        b.WriteConfiguration(input.writefile);
    }

    return 0;
}
