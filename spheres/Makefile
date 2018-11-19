
CXX = g++

CXXFLAGS = -Iinclude

TARGET = spheres

SRCDIR = ./src/
SRCS = $(wildcard $(SRCDIR)*.C)

OBJDIR = ./bin/
OBJS = $(addprefix $(OBJDIR), $(notdir $(patsubst %.C, %.o, $(SRCS))))

ifeq ($(wildcard $(OBJDIR)), )
$(shell mkdir -p $(OBJDIR))
endif

$(TARGET):$(OBJS)
	$(CXX) $^ -o $@

$(OBJDIR)%.o:$(SRCDIR)%.C
	$(CXX) -c $< -o $@ $(CXXFLAGS)

OUTDIR = ./output/
ifeq ($(wildcard $(OUTDIR)), )
$(shell mkdir -p $(OUTDIR))
endif

run_windows:
	./spheres.exe ./input/input.txt

run_hpc:
	./spheres ./input/input.txt

run:
	sbatch spheres.slurm

.PHONY:clean
clean:
	-rm -rf $(TARGET) bin $(OUTDIR)*