#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <mpi.h>

inline long exact(const long N){
    return N*(N+1)/2;
}

void reduce_mpi(const int rank, long& mysum){
    double sum;
    MPI_Reduce(&mysum, &sum, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD); //collect partial results from each rank to rank 0.
}

void reduce_manual(int rank, int size, long& sum){
    // Size is a power of 2 for simplicity
    const int TAG = 1337;

    for(int send_rec_border = size/2; send_rec_border >= 1; send_rec_border/=2)
    {
        if(rank < send_rec_border)
        {
            long buffer;
            MPI_Recv(&buffer, 1, MPI_LONG, rank+send_rec_border, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += buffer;
        } else if(rank <= send_rec_border*2){
            MPI_Send(&sum, 1, MPI_LONG, rank-send_rec_border, TAG, MPI_COMM_WORLD);
        }
    }
}


int main(int argc, char** argv){
    const long N = 1000000;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    // -------------------------
    // Case 1: naive reduction
    // -------------------------
/*
    long sum = 0;
    reduce_mpi(rank, sum);
*/

 
    // -------------------------
    // Case 2: distributed reduction
    // -------------------------

    // Determine work load per rank
    long N_per_rank = N / size;

    // Remark: Start summation at 1!!! (1+2+3+...+N)
    long N_start = rank * N_per_rank + 1;
    long N_end = (rank+1) * N_per_rank;
    // the last rank has to do a few additions more if size does not divide N
    if(rank == size-1){
        N_end += N % size;
    }

    // Compute partial sum: N_start + (N_start+1) + ... + (N_start+N_per_rank-1)
    long sum = 0;
    for(long i = N_start; i <= N_end; ++i){
        sum += i;
    }

    // Reduction
    reduce_manual(rank, size, sum);

 
    // -------------------------
    // Print the result
    // -------------------------
    if(rank == 0){
        std::cout << std::left << std::setw(25) << "Final result (exact): " << exact(N) << std::endl;
        std::cout << std::left << std::setw(25) << "Final result (MPI): " << sum << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();
    
    return 0;
}
