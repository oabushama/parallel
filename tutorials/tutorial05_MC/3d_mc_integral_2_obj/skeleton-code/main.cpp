//#include <stdio.h>
//#include <stdlib.h>
#include <random>
#include <omp.h>

//global variables

//cube
double x1 = 1.5;
double y1_ = 2.5;
double z1 = 0.5;
double cube_l = 3.0*0.5;

//blue circle
double x2 = 4.0;
double y2 = 2.0;
double z2 = -0.5;
double blue_r = 2.0;


bool inside_cube(double x, double y, double z){
	//TODO: return true if x,y,z are inside the cube
	bool x_dim = x1 - cube_l < x && x < x1 + cube_l;
	bool y_dim = y1_ - cube_l < y && y < y1_ + cube_l;
	bool z_dim = z1 - cube_l < z && z < z1 + cube_l;

	return x_dim && y_dim && z_dim;
}

bool inside_blue_circle(double x, double y, double z){
	//TODO: return true if x,y,z are inside the blue circle
	return (x - x2)*(x - x2) + (y - y2)*(y - y2) + (z - z2)*(z - z2) < blue_r*blue_r;
}


unsigned overlapMC(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, size_t n)
{
	unsigned pts_inside = 0;

	#pragma omp parallel reduction(+:pts_inside)
	{
		int t = omp_get_thread_num();
		std::default_random_engine g(42 + t);  // random generator with seed 42
		
		//TODO: generate random uniform random sample generators here
		std::uniform_real_distribution<double> u_x(x_min, x_max);
		std::uniform_real_distribution<double> u_y(y_min, y_max);
		std::uniform_real_distribution<double> u_z(z_min, z_max);

		#pragma omp for
		for (size_t i = 0; i < n; ++i)
		{
			// TODO: implement the MC integration part here!
			double x = u_x(g);
			double y = u_y(g);
			double z = u_z(g);

			if(inside_cube(x,y,z) && inside_blue_circle(x,y,z)){
				++pts_inside;
			}
		}
	}

	return  pts_inside;
}




int main(){

	//TODO: compute the min and max distance of all three coordinates
	// these will be used to determine the sample area
	double x_min = std::min(x1 - cube_l, x2 - blue_r);
	double x_max = std::max(x1 + cube_l, x2 + blue_r);

	double y_min = std::min(y1_ - cube_l, y2 - blue_r);
	double y_max = std::max(y1_ + cube_l, y2 + blue_r);

	double z_min = std::min(z1 - cube_l, z2 - blue_r);
	double z_max = std::max(z1 + cube_l, z2 + blue_r);

	// TODO: compute the cuboid area for which you uniformly sample
	const double area_cuboid = (x_max - x_min)*(y_max - y_min)*(z_max - z_min);

	size_t n = 1e8;// default number of MC samples

	double ref = 3.4; // reference solution

	double t0 = omp_get_wtime();

	unsigned sum = overlapMC(x_min, x_max, y_min, y_max, z_min, z_max, n);

	//TODO: compute the integrated area
	double area = sum/double(n) * area_cuboid;

	double t1 = omp_get_wtime();

	double error = std::abs(area - ref);
	if(error > 1e-2){
		printf("ERROR: you should get %.20f, but got: %.20f\n", ref, area);
	}
	printf("result:  %.20f\nref: %.20e\ntime: %.20f\n", area, ref, t1 - t0);

	return 0;
}
