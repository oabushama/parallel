#include <chrono>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <vector>



struct Diagnostics {
    double time;
    double concentration;
    std::vector<int> histogram;
    
    Diagnostics(double time, double concentration, std::vector<int> &histogram)
        : time(time), concentration(concentration), histogram(histogram) {}
};



struct Diffusion {
    const double D;   // diffusion constant
    const double L;   // domain length
    size_t N;         // grid points per direction (whole grid is NxN)
    double h, dt;     // grid spacing and timestep
    double aux;       // auxiliary variable (prefactor)
    
    std::vector<double> c;     // solution vector
    std::vector<double> c_tmp; // temporary solution vector
    std::vector<Diagnostics> diag; // vector to store concentration values

    
    Diffusion(double D, double L, size_t N) : D(D), L(L), N(N) {
        h = L / (N - 1);
        // TODO: 1b replace constant value with dt_{max}
        dt = h * h / (4.0 * D);
        
        c.resize((N + 2), 0.0); // +2 for the ghost cells
        c_tmp.resize((N + 2), 0.0);
        
        aux = dt * D / (h * h);
        initialize_density();
    }


    void advance() {
    /* Central differences in space, forward Euler in time, Dirichlet BCs */
    // TODO: 1c Implement central difference in space, forward Euler in time
    // TODO: 1e Parallelize diffusion with OpenMP
#pragma omp parallel for 
        for (size_t j = 1; j <= N; ++j)
            c_tmp[j] = c[j] + aux * ( c[j+1] - 2*c[j] + c[j-1] );

        std::swap(c_tmp, c);
    }


    double compute_diagnostics(const double t) {
        double amount = 0.0;
        /* Integration to compute total concentration */
        // TODO: 1f Parallelize integration with OpenMP
#pragma omp parallel for reduction(+ : amount)
        for (size_t j = 1; j <= N; ++j)
            amount += c[j];

        amount *= h;  // Concentration = Mass / Volume

        printf("t = %lf amount = %lf\n", t, amount);
        return amount;
    }


    void compute_histogram(std::vector<int> &hist) {
        /* number of bins */
        const size_t M = hist.size();

        /* Initialize max and min concentration */
        double max_c, min_c;
        max_c = c[1];
        min_c = c[1];

        /* Find max and min concentration values */
        // TODO: 1f Parallelize max_c and min_c initialization
#pragma omp parallel for reduction(min : min_c) reduction(max : max_c)
        for (size_t j = 1; j <= N; ++j) {
            double c0 = c[j];
            if (c0 > max_c)
                max_c = c0;
            if (c0 < min_c)
                min_c = c0;
        }

        const double epsilon = 1e-8;
        double dc = (max_c - min_c + epsilon) / M;

        /* Accumulate equispaced bins */
        // TODO: 1f Parallelize bin accumulation
#pragma omp parallel
        {
            // local calculation of bins
            std::vector<int> local_hist(M, 0);
#pragma omp for nowait // nowait optional
            for (size_t j = 1; j <= N; ++j) {
                size_t bin = (c[j] - min_c) / dc;
                local_hist[bin]++;
            }
        // global aggregation of local histograms
#pragma omp critical
        for (size_t i = 0; i < M; ++i)
            hist[i] += local_hist[i];
        }


        /* ALTERNATIVE NAIVE SOLUTION (1/2pts, very slow)
         *
         * Accumulate equispaced bins
        // TODO: 1f Parallelize bin accumulation
#pragma omp parallel
        for (size_t j = 1; j <= N; ++j) {
            size_t bin = (c[j] - min_c) / dc;
#pragma omp atomic
            hist[bin]++;
        }
        */
    }


    void initialize_density() {
        const double bound = 0.25 * L;

        for (size_t j = 0; j < N; ++j) {
          if (std::abs(j * h - 0.5 * L) < bound)
            c[j+1] = 1;
          else
            c[j+1] = 0;
        }
    }


    void write_diagnostics(const std::string &filename) const {
        // write header
        std::ofstream out_file(filename, std::ios::out);
        out_file << "t concentration";
        for (size_t i = 0; i < diag[0].histogram.size(); ++i)
            out_file << " bin" << i;
        out_file << '\n';

        // write data
        for (const Diagnostics &d : diag) {
            out_file << d.time << ' ' << d.concentration;
            for (auto h : d.histogram)
                out_file << ' ' << h;
            out_file << '\n';
        }
        out_file.close();
  }

}; // end struct Diffusion



int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s N out\n", argv[0]);
        return 1;
    }

    /* Set up diffusion parameters */
    //
    /* Init diffusion constant */
    const double D = 1;

    /* Init domain size */
    const double L = 1;

    /* Init simulation length */
    const double T = 1;

    /* Init num grid points */
    const size_t N = std::stoul(argv[1]);

    /* Produce output */
    const size_t out = std::stoul(argv[2]); // YES == 1

    int numThreads = omp_get_num_threads();
    printf("Running Diffusion 2D on a %zu x %zu grid with %d threads.\n", N, N,
           numThreads);


    /* Initialize diffusion system */
    Diffusion system(D, L, N);


    /* System evolution in time */
    size_t numSteps = (T / system.dt + 1);
    auto tstart = std::chrono::steady_clock::now();
    for (size_t step = 0; step <= numSteps; ++step) {
        double t = system.dt * step;
        double amount = system.compute_diagnostics(t);

        std::vector<int> hist(10, 0); // histogram with 10 bins
        system.compute_histogram(hist);

        system.diag.push_back(Diagnostics(t, amount, hist));

        system.advance();
    }
    auto tend = std::chrono::steady_clock::now();
    double ms = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();


    /* Write diagnostics */
    if (out == 1)
        system.write_diagnostics("diagnostics.dat");

    printf("time: %lf\n", ms);

    return 0;
}
