
CXX = g++ -std=c++11

CXXFLAGS = -Iinclude

TARGET = discretize

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
	./discretize.exe

run_hpc:
	./discretize

run:
	sbatch discretize.slurm

.PHONY:clean
clean:
	-rm -rf $(TARGET) bin $(OUTDIR) *.err *.out