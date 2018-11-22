gpu:
	sbatch tf_gpu.slurm

cpu:
	sbatch tf_cpu.slurm

fat:
	sbatch tf_fat.slurm

test:
	sbatch tf_gpu_inference.slurm

.PHONY:clean, delete_log, delete_model
clean:
	-rm *.err *.out
clean_all:
	-rm *.err *.out *.txt
	-rm -rf Model
