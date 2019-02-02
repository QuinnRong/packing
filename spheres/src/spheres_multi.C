#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string.h>

#include "box.h"
#include "read_input.h"

double rand_range(double r1, double r2)
{
    double tmp = rand() * 1.0 / RAND_MAX;
    return r1 + tmp * (r2 -r1);
}

int main(int argc, char **argv)
{
    read_input input;
    input.read(argc, argv);

    std::ofstream summary("./output/summary.txt");
    summary << "filename    N    radius     pf" << std::endl;
    
    double m_pf[3] = {0.1, 0.35, 0.6};
    srand (time(NULL));
    for (int i = 0; i < 3; ++i)
    {
        input.maxpf = m_pf[i];

        sprintf(input.writefile, "./output/struct_%d.dat", i);
        sprintf(input.datafile, "./output/statis_%d.dat", i);

        double r = pow(input.initialpf*pow(SIZE, DIM)/(input.N*VOLUMESPHERE), 1.0/((double)(DIM)));

        box b(input.N, r, input.growthrate, input.maxpf);

        b.CreateSpheres(input.temp);

        std::ofstream output(input.datafile);
        output.precision(16);  

        output << "step packing-fraction pressure energy-change total-events" << std::endl;
        int step = 0;
        while ((b.pf < input.maxpf) && (b.pressure < input.maxpressure)) 
        {
            // printf("step = %4d, pf = %.4f, pressure = %.4f\n", step, b.pf, b.pressure);
            b.Process(input.eventspercycle*input.N);
            output << step++ << " " << b.pf << " " << b.pressure << " "
                   << b.energychange << " " << b.neventstot << " " << std::endl;
            b.Synchronize(true);
        }
        printf("\n%d: %.4f -> final radius: %.6f\n", i, input.maxpf, b.r);
        summary << "struct_" << i << ".dat " << input.N << " " <<  b.r 
                << " " << b.pf << " " << std::endl;

        output.close();
        b.WriteConfiguration(input.writefile);
    }

    summary.close();

    return 0;
}
