/* openmp version of N body simulation */
/* Distributed and Parallel Computing 2nd project */
/* Author: Yangkai Zhou, Yimeng Li */

# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <vector>
# include <cstdlib>
# include <omp.h>
 
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

/* define the structure for force */
struct force
{
	double f_x;
	double f_y;

	force(double f_x = 0, double f_y = 0)
	{
		this->f_x = f_x;
		this->f_y = f_y;
	}
};

/* define the global data */
int g_N;						// number of particles
vector<force> g_F;				// net force for each particle
vector<particle> g_pv;			// particle vector

void setUp();
void update();

int main(int argc, char ** argv) {
	
	/* initialize */
	setUp();

	int time = 0;
	
	clock_t start, finish;
	start = clock();
	while (time < 200) {

		update();
		++time;
	}
	finish = clock();
	cout << "Execution Time: " << (double)(finish-start)/CLOCKS_PER_SEC << endl;
	
	return 0;
}

/* read the input data */
void setUp()
{	
	ifstream inFile;
	inFile.open("input_1.txt");
	
	inFile >> g_N;
	g_pv.resize(g_N);
	g_F.resize(g_N);
	for ( int i = 0; i < g_N; ++i )
	{
		inFile >> g_pv[i].m >> g_pv[i].pos_x >> g_pv[i].pos_y
			   >> g_pv[i].v_x >> g_pv[i].v_y >> g_pv[i].a_x >> g_pv[i].a_y; 
	}
	inFile.close();
}

/* update one frame */
void update()
{
	int i, j;
	omp_set_num_threads(16);
	# pragma omp parallel private(j)
	{
		/* compute the force */
		# pragma omp for schedule(dynamic)
			for ( i = 0; i < g_N; ++i ) {
				for ( j = 0; j < g_N; ++j ) {
					if (i == j)
						continue;
					double r_2 = pow((g_pv[i].pos_x - g_pv[j].pos_x),2) + pow((g_pv[i].pos_y - g_pv[j].pos_y),2);
				g_F[i].f_x += (-1) * G * g_pv[i].m * g_pv[j].m * (g_pv[i].pos_x - g_pv[j].pos_x) / (pow(r_2 + e,1.5));
				g_F[i].f_y += (-1) * G * g_pv[i].m * g_pv[j].m * (g_pv[i].pos_y - g_pv[j].pos_y) / (pow(r_2 + e,1.5));
				} 
			}
			
		# pragma omp for schedule(dynamic)
			for (int i = 0; i < g_N; ++i) {
				g_pv[i].a_x = g_F[i].f_x / g_pv[i].m;
				g_pv[i].a_y = g_F[i].f_y / g_pv[i].m;
			}	

		/* compute the velocity */
		# pragma omp for schedule(dynamic)
			for ( int i = 0; i < g_N; ++i ) {
				g_pv[i].v_x += g_pv[i].a_x * period;
				g_pv[i].v_y += g_pv[i].a_y * period;
			}
		# pragma omp for schedule(dynamic)
			for ( int i = 0; i < g_N; ++i ) {
				g_pv[i].pos_x += g_pv[i].v_x * period;
				g_pv[i].pos_y += g_pv[i].v_y * period;
			}
	}
}
