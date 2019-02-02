
TARGET = emt

g++ -std=c++11 emt.cpp -o $(TARGET)

OUTDIR = ./output/
ifeq ($(wildcard $(OUTDIR)), )
$(shell mkdir -p $(OUTDIR))
endif

run: $(TARGET)
	./$(TARGET)

.PHONY:clean
clean:
	-rm -rf $(TARGET) $(OUTDIR)