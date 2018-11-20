all:
	g++ -std=c++11 main.cpp discretize.cpp
	sbatch discretize.slurm

clean:
	rm -rf output/* *.err *.out
