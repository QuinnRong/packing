
CXX = g++

TARGET = spheres

SRCDIR = ./
SRCS = $(wildcard $(SRCDIR)*.C)

OBJDIR = ./bin/
OBJS = $(addprefix $(OBJDIR), $(notdir $(patsubst %.C, %.o, $(SRCS))))

ifeq ($(wildcard $(OBJDIR)), )
$(shell mkdir -p $(OBJDIR))
endif

$(TARGET):$(OBJS)
	$(CXX) $^ -o $@  $(LIBS)

$(OBJDIR)%.o:$(SRCDIR)%.C
	$(CXX) -c $< -o $@ $(CXXFLAGS)

run:
	./spheres.exe input.txt

.PHONY:clean
clean:
	-rm -rf $(TARGET) bin