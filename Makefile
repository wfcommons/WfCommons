# Makefile 

objects = bin/cpu-benchmark.o
CXX= g++
CPPFLAGS= -std=c++11
execname = bin/cpu-benchmark

# compile
$(execname): $(objects)
	$(CXX) $(CPPFLAGS) -o $(execname) $(objects)

#clean Makefile
clean:
	rm -rf $(objects) $(execname)


# #define variables
# objects= gpu-benchmark.o kernels.o 
# NVCC= nvcc               #cuda c compiler
# CPPFLAGS= -std=c++11
# opt= -O2                 #optimization flag
# LIBS=  
# execname= gpu-benchmark

# .PHONY: clean

# #compile
# $(execname): $(objects)
# 	$(NVCC) $(CPPFLAGS) $(opt) -o $(execname) $(objects) $(LIBS) 

# kernels.o: kernels.cu
# 	$(NVCC) $(CPPFLAGS) $(opt)  -c kernels.cu
# gpu-benchmark.o: gpu-benchmark.cu
# 	$(NVCC) $(CPPFLAGS) $(opt)  -c gpu-benchmark.cu


# #clean Makefile
# clean:
# 	rm $(objects)
# #end of Makefile
