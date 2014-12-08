/* second version of N body simulation using CUDA */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <ctime>

using namespace std;


/* define the global constants */
const double G = 6.67 * pow(10, -11);
const double e = 0.00001;
const double period = 1;


/* define the structure of particle */
struct particle
{
	double m;
	double pos_x;
	double pos_y;
	double v_x;
	double v_y;
	double a_x;
	double a_y;

	particle(double m = 0, double pos_x = 0, double pos_y = 0, 
			double v_x = 0, double v_y = 0, double a_x = 0, double a_y = 0)
	{
		this->m			= m;
		this->pos_x		= pos_x;
		this->pos_y		= pos_y;
		this->v_x		= v_x;
		this->v_y		= v_y;
		this->a_x		= a_x;
		this->a_y		= a_y;
	}
};


struct my_double2
{
	double x, y;
	
	__device__ my_double2(double x = 0, double y = 0)
	{
		this->x = x;
		this->y = y;
	}
};



/* define the global data */
int g_N;									// number of particles
int g_P;									// number of particles in a tile
thrust::host_vector<particle> g_pv;			// particle vector


void setUp();



/* calculate the interaction between two bodies */
__device__ my_double2 bodyBodyAcceleration(double G, double e, particle b1, particle b2, my_double2 acceleration)
{	
	double r_2 = pow((b1.pos_x - b2.pos_x),2) + pow((b1.pos_y - b2.pos_y),2);
	b1.a_x = (-1) * G * b2.m * (b1.pos_x - b2.pos_x) / (pow(r_2 + e, 1.5));
	b1.a_y = (-1) * G * b2.m * (b1.pos_y - b2.pos_y) / (pow(r_2 + e, 1.5));
	
	acceleration.x += b1.a_x;
	acceleration.y += b1.a_y;
	return acceleration;
}


/* calculate the interaction inside a P*P block */
__device__ my_double2 tileAcceleration(double G, double e, particle b, my_double2 acceleration)
{
	extern __shared__ particle shParticles[];
	
	
	for(int i = 0; i < blockDim.x; ++i)
	{
		acceleration = bodyBodyAcceleration(G, e, b, shParticles[i], acceleration);
	}
	
	return acceleration;
}


/* update the position */
__device__ void updatePosition(double period, particle* particle_arr)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	/* compute the velocity */
	particle_arr[idx].v_x += particle_arr[idx].a_x * period;
	particle_arr[idx].v_y += particle_arr[idx].a_y * period;
	
	/* compute the new position */

	particle_arr[idx].pos_x += particle_arr[idx].v_x * period;
	particle_arr[idx].pos_y += particle_arr[idx].v_y * period;
}


/* calculate the whole acceleration */
__global__ void updateScene(int N, int P, double G, double e, double period, particle* particle_arr)
{
	extern __shared__ particle shParticles[];
	
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	particle ptc = particle_arr[id];
	
	my_double2 acceleration;
	acceleration.x = ptc.a_x;
	acceleration.y = ptc.a_y;
	
	int i, tile;
	for(i = 0, tile = 0; i < N; i += P, ++tile)
	{
		/* fill in the shared memory */
		int idx = tile * blockDim.x + threadIdx.x;
		shParticles[threadIdx.x] = particle_arr[idx];
		__syncthreads();
		
		/* calculate the acceleration with a tile */
		acceleration = tileAcceleration(G, e, ptc, acceleration);
		__syncthreads();	
	}
	
	ptc.a_x = acceleration.x;
	ptc.a_y = acceleration.y;
	
	updatePosition(period, particle_arr);
}



int main(int argc, char ** argv) {

	setUp();
	g_P = static_cast<int>(sqrt(g_N)) + 1;
	
	/* device copy of particle array */
	thrust::device_vector<particle> d_particle_arr = g_pv;
	
	/* get the raw pointer of particle array */
	particle *particle_arr = thrust::raw_pointer_cast(d_particle_arr.data());
	
	clock_t start, finish;
	start = clock();
	
	int time = 0;
	while(time < 100000)
	{	
		
		updateScene<<<g_P,g_P,g_P*sizeof(particle)>>>(g_N, g_P, G, e, period, particle_arr);
		/*
		g_pv = d_particle_arr;
		
		
		for ( int i = 0; i < g_N; ++i )
		{
			cout << "particle: " << i << " pos_x: " << g_pv[i].pos_x << " pos_y: " << g_pv[i].pos_y << endl;
		}
		*/
		time++;
	}
	
	
	finish = clock();
	cout << "Execution Time: " << (double)(finish-start)/CLOCKS_PER_SEC << endl;
	
	return 0;
}



/* read the input data */
void setUp()
{
	ifstream inFile;
	inFile.open("input.txt");
	
	inFile >> g_N;
	g_pv.resize(g_N);
	for ( int i = 0; i < g_N; ++i )
	{
		inFile >> g_pv[i].m >> g_pv[i].pos_x >> g_pv[i].pos_y
			   >> g_pv[i].v_x >> g_pv[i].v_y >> g_pv[i].a_x >> g_pv[i].a_y; 
	}
	
	inFile.close();
}

