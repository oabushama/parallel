#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <fstream>
#include <string>
#include <mpi.h>


unsigned overlapMC(const double x2, const double R1, const double R2, size_t n, int rank, int procs)
{
	unsigned pts_inside = 0;

	std::default_random_engine g(42 + rank);  // random generator with seed 42
	std::uniform_real_distribution<double> u; // uniform distribution in [0, 1]

	// TODO_b: split the amount of work as equally as possible for each process.
	size_t n_local_size;
	size_t n_start = 0;
	size_t n_end   = n;

	for (size_t i = n_start; i < n_end; ++i)
	{
		// TODO_a: implement the MC integration part here!

	}

	return  pts_inside;
}




int main(int argc, char *argv[])
{
	// TODO_b: Start-up the MPI environment and determine this process' rank ID as
    // well as the total number of processes (=procs) involved in the
    // communicator
    int rank, procs = -1;


	const double R1 = 5.0;		// Radius of first circle
	const double R2 = 10.0;		// Radius of second circle
	const double x2 = 12.0;		// x2 coordinate of second circle center

	// TODO_a: calculate the rectangle area for which you uniformly sample x & y
	const double area_rectangle;

	size_t n = 1e9 + 1;// default number of MC samples

	double ref = 17.0097776; // reference solution

	double t0 = MPI_Wtime(); 

	unsigned local_sum = overlapMC(x2, R1, R2, n, rank, procs);

	unsigned global_sum = (procs == -1) ? local_sum : 0;

	// TODO_b: Sum up and average all the local_sums to the master ranks global_sum using MPI_Reduce


	double area = global_sum / double(n) * area_rectangle;

	double t1 = MPI_Wtime();

	if(rank == 0){
		double error = std::abs(area - ref);
		if(error > 1e-2){
			printf("ERROR: you should get pi, but got: %.20f\n", area);
		}
		else{
			printf("result:  %.20f\nref: %.20e\ntime: %.20f\n", area, ref, t1 - t0);

			std::string file_name = "out/";
			file_name += std::to_string(procs);
			file_name += ".txt";
			//output time in a file
			std::ofstream myfile;
			myfile.open (file_name);
			myfile << procs << " " << t1 - t0 << "\n";
			myfile.close();
		}
	}

	// TODO_b: Shutdown the MPI environment


	return 0;
}
