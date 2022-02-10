#include <iostream>
#include <stdlib.h>
#include <time.h>

/**This function randomly accesses positions in an array much bigger than the cache and sums 1 unit to it.*/
void mem_access(long work) {
    long max = 7000000000;
    long *arr = new long[max]; 
    long index;
    srand (time(NULL));
    for(long i=0; i<work; i++){
        arr[rand() % max]++; 
    }
}
    
int main(int argc, char **argv) {

  
    // Process command-line args
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <work (# M memory accesses)> \n";
        exit(1);
    }
    
    long work;
    try {
       work = std::stol(argv[1]);
    } catch (std::invalid_argument &e) {
        std::cerr << "Invalid argument: " << e.what() << "\n";
        exit(1);
    }
    mem_access(1000000*work);
    
    exit(0);

}
