
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

run:
	./spheres.exe ./input/input.txt

.PHONY:clean
clean:
	-rm -rf $(TARGET) bin output/*