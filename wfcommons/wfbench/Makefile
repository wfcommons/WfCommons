#Makefile 
#define variables
objects= pi-cuda.o kernels.o 
NVCC= nvcc               #cuda c compiler
CPPFLAGS= -std=c++11
opt= -O2                 #optimization flag
LIBS=  
execname= pi-cuda

.PHONY: clean

#compile
$(execname): $(objects)
	$(NVCC) $(CPPFLAGS) $(opt) -o $(execname) $(objects) $(LIBS) 

kernels.o: kernels.cu
	$(NVCC) $(CPPFLAGS) $(opt)  -c kernels.cu
pi-cuda.o: pi-cuda.cu
	$(NVCC) $(CPPFLAGS) $(opt)  -c pi-cuda.cu


#clean Makefile
clean:
	rm $(objects)
#end of Makefile
