#include <vector>
#include <thread>
#include <iostream>
#include <random>
#include <cmath>
#include <functional>

#define PRECISION 100000L

/**
 * This function computes pi using work trials for a Monte Carlo method. It's
 * GOOD because it uses a good number generator.  Unfortunately, this method
 * leads to extra memory references that cause extra RAM references, and thus
 * extra cache misses.
 */

void compute_good_pi(long work) {

    std::uniform_real_distribution<double> random_dist(-0.5, 0.5);
    std::mt19937 rng;
    rng.seed((long)&work);

    double good_pi = 0.0;
    double x,y;
    for (long sample=0; sample < work; sample++) {
        x = random_dist(rng);
        y = random_dist(rng);
      good_pi += (double)(std::sqrt(x*x + y*y) < 0.5);
    } 
//    std::cout << "good pi = " << (good_pi/(double)work)/(0.5*0.5) << "\n";
         
}         


/**
 * This function computes pi using work trials for a Monte Carlo method. It's
 * TERRIBLE because it uses a custom, bad, number generator, which has too 
 * much bias to compute a good value a PI. The reason for using the generator
 * is that it does not cause extra memory references, and thus keeps this benchmark
 * 100% CPU intensive. 
 */
void compute_terrible_pi(long work) {

    long rng = (long)&work;
    double terrible_pi = 0.0;
    double x_value, y_value;
    for (long sample=0; sample < work; sample++) {
        rng = (((rng * 214013L + 2531011L) >> 16) & 32767);
        x_value = -0.5 + (rng % PRECISION) / (double)PRECISION;
        rng = (((rng * 214013L + 2531011L) >> 16) & 32767);
        y_value = -0.5 + (rng % PRECISION) / (double)PRECISION;
        terrible_pi += (double)(std::sqrt(x_value*x_value+ y_value*y_value) < 0.5);
    }
 //   std::cout << "terrible pi = " << (terrible_pi/(double)work)/ (0.5*0.5) << "\n";
}

int main(int argc, char **argv) {

  
    // Process command-line args
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <work (# 1M samples)>\n";
        exit(1);
    }
    
    long work;
    int num_threads;
    try {
       work = std::stol(argv[1]);
    //    num_threads = std::stoi(argv[2]);
    } catch (std::invalid_argument &e) {
        std::cerr << "Invalid argument: " << e.what() << "\n";
        exit(1);
    }

    // Create all threads
    // std::vector<std::thread> workers;
    // for (unsigned int i = 0; i < num_threads; i++) {
    //    auto t = std::thread([work]() {
    compute_terrible_pi(1000000*work);
    std::cout<<"Pi computed!"<<std::endl;
    //compute_good_pi(1000000*work);
        //   });
    //    workers.push_back(std::move(t));
    // }

    // Wait for all threads
    // for (unsigned int i = 0; i < num_threads; i++) {
    //     workers.at(i).join();
    // }

    exit(0);

}
