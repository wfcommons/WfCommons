/**
 * Copyright (c) 2021 The WfCommons Team.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include <algorithm>
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

#define PRECISION 100000L

/**
 * This function computes pi using work trials for a Monte Carlo method. It's
 * GOOD because it uses a good number generator.  Unfortunately, this method
 * leads to extra memory references that cause extra RAM references, and thus
 * extra cache misses.
 */
void computeGoodPi(long work) {
    std::uniform_real_distribution<double> random_dist(-0.5, 0.5);
    std::mt19937 rng;
    rng.seed((long) &work);

    double good_pi = 0.0;
    double x, y;
    for (long sample = 0; sample < work; sample++) {
        x = random_dist(rng);
        y = random_dist(rng);
        good_pi += (double) (std::sqrt(x * x + y * y) < 0.5);
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
void computeTerriblePi(unsigned long work) {
    unsigned long rng = (unsigned long) &work;
    double terrible_pi = 0.0;
    double x_value, y_value;
    for (unsigned long sample = 0; sample < work; sample++) {
        rng = (((rng * 214013L + 2531011L) >> 16) & 32767);
        x_value = -0.5 + (rng % PRECISION) / (double) PRECISION;
        rng = (((rng * 214013L + 2531011L) >> 16) & 32767);
        y_value = -0.5 + (rng % PRECISION) / (double) PRECISION;
        terrible_pi += (double) (std::sqrt(x_value * x_value + y_value * y_value) < 0.5);
    }
    //   std::cout << "terrible pi = " << (terrible_pi/(double)work)/ (0.5*0.5) << "\n";
}

/**This function randomly accesses positions in an array much bigger than the cache and sums 1 unit to it.*/
void computeMemAccess(unsigned long work) {
    //unsigned long max = 7000000000;
    unsigned long max = 300000000;
    unsigned long *arr = new unsigned long[max];
    srandom(time(nullptr));
    for (unsigned long i = 0; i < work; i++) {
        arr[random() % max]++;
    }
}

long getCmdOption(char **begin, char **end, const std::string &option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return std::stol(*itr);
    }
    return 0;
}

int main(int argc, char **argv) {

    auto start = std::chrono::system_clock::now();

    // process command-line args
    long cpu_work = getCmdOption(argv, argv + argc, "--cpu-work");
    long mem_work = getCmdOption(argv, argv + argc, "--mem-work");

    // sanity check
    if (cpu_work > 0 && mem_work > 0) {
        std::cerr << "ERROR: Cannot run CPU and Memory benchmark at the same time.\n"
                  << "Usage: " << argv[0] << " [--cpu-work|--mem-work] <work [#1M samples|#K memory accesses]>"
                  << std::endl;
        exit(1);

    } else if (cpu_work == 0 && mem_work == 0) {
        std::cerr << "ERROR: Work should be provided to either CPU or Memory.\n"
                  << "Usage: " << argv[0] << " [--cpu-work|--mem-work] <work [#1M samples|#K memory accesses]>"
                  << std::endl;
        exit(1);

    } else if (cpu_work < 0 || mem_work < 0) {
        std::cerr << "ERROR: Work should be an integer larger than 0.\n"
                  << "Usage: " << argv[0] << " [--cpu-work|--mem-work] <work [#1M samples|#K memory accesses]>"
                  << std::endl;
        exit(1);
    }

    // compute benchmark
    if (cpu_work > 0) {
        computeTerriblePi(1000000 * cpu_work);
    } else {
        computeMemAccess(1000 * mem_work);
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::string type = cpu_work > 0 ? "CPU" : "Memory";

    std::cout << type << ": " << elapsed_seconds.count() << std::endl;

    exit(0);
}
