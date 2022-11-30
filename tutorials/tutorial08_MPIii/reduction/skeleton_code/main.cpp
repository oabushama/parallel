#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

inline long exact(const long N){
    return N*(N+1)/2.0;
}

void reduce_mpi(const int rank, long& sum){
    // TODO: Perform the reduction using blocking collectives.
}

void reduce_manual(int rank, int size, long& sum){
    // TODO: Implement a tree based reduction using blocking point-to-point communication.
    // Size is a power of 2 for simplicity
}


int main(int argc, char** argv){
    const long N = 1000000;

    long sum = 0;

    sum = exact(N);

    std::cout << std::left << std::setw(25) << "Final result (exact): " << sum << std::endl;

    return 0;
}
