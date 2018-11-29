g++ -Iinclude -std=c++11 ./src/main.cpp ./src/discretize.cpp -o discretize

run:
	sbatch discretize.slurm

clean:
	-rm -rf output/* *.err *.out discretize
