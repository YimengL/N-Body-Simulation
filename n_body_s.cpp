/* sequential version of N body simulation */
/* Distributed and Parallel Computing 2nd project */
/* Author: Yangkai Zhou, Yimeng Li */

# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <vector>
# include <cstdlib>

# include <GL/glew.h>
# include <GL/freeglut.h>
# include <GL/gl.h>

using namespace std;

/* the size of the opengl window */
# define WIDTH 400
# define HEIGHT 400

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

/* prototype of several helper function */
void setUp();
void update();
void display(void);
void timer(int value);

int main(int argc, char ** argv) {
	
	/* set up the window */
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (WIDTH, HEIGHT);
	glutInitWindowPosition (100, 100);
	glutCreateWindow ("Example");

	/* initialize */
	setUp();
	
	/* These three functions will ensure that display is called periodically, so that the particles are redrawn
	 * in their current position */
	glutDisplayFunc(display);
	glutTimerFunc(10, timer, 1);
	glutMainLoop();

	return 0;
}

void timer(int value) {
	update();
  	glutPostRedisplay();
  	glutTimerFunc(30, timer, value);
}

/* read the input data */
void setUp()
{
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	// Sets a 2d projection matrix
	// (0,0) is the lower left corner (WIDTH, HEIGHT) is the upper right
	glOrtho (0, WIDTH, 0, HEIGHT, 0, 1);
	glDisable(GL_DEPTH_TEST);
	
	ifstream inFile;
	inFile.open("input.txt");
	
	inFile >> g_N;
	g_pv.resize(g_N);
	g_F.resize(g_N);
	for (int i = 0; i < g_N; ++i) {
		inFile >> g_pv[i].m >> g_pv[i].pos_x >> g_pv[i].pos_y
			   >> g_pv[i].v_x >> g_pv[i].v_y >> g_pv[i].a_x >> g_pv[i].a_y; 
	}	
	inFile.close();
}

void display(void){
	glClear(GL_COLOR_BUFFER_BIT);
  	for(int i = 0; i < g_pv.size(); ++i) {
	    /* Get the ith particle */
	    particle p = g_pv[i];

	    /* Draw the particle as a little square. */
	    glBegin(GL_QUADS);
	    glColor3f (1.0, 1.0, 1.0);
	    glVertex2f(p.pos_x + 1.5, p.pos_y + 1.5);
	    glVertex2f(p.pos_x - 1.5, p.pos_y + 1.5);
	    glVertex2f(p.pos_x - 1.5, p.pos_y - 1.5);
	    glVertex2f(p.pos_x + 1.5, p.pos_y - 1.5);
	    glEnd();
  	}
  glFlush ();
}

/* update one frame */
void update()
{
	/* compute the force */
	for (int i = 0; i < g_N; ++i) {
		for (int j = 0; j < g_N; ++j) {
			if (i == j)
				continue;
			else {
				double r_2 = pow((g_pv[i].pos_x - g_pv[j].pos_x),2) + pow((g_pv[i].pos_y - g_pv[j].pos_y),2);
				g_F[i].f_x += (-1) * G * g_pv[i].m * g_pv[j].m * (g_pv[i].pos_x - g_pv[j].pos_x) / (pow(r_2 + e,1.5));
				g_F[i].f_y += (-1) * G * g_pv[i].m * g_pv[j].m * (g_pv[i].pos_y - g_pv[j].pos_y) / (pow(r_2 + e,1.5));
			}
		} 
	}

	/* using force to compute acceleration */
	for (int i = 0; i < g_N; ++i) {
		g_pv[i].a_x = g_F[i].f_x / g_pv[i].m;
		g_pv[i].a_y = g_F[i].f_y / g_pv[i].m;
	}

	/* compute the velocity */
	for (int i = 0; i < g_N; ++i) {
		g_pv[i].v_x += g_pv[i].a_x * period;
		g_pv[i].v_y += g_pv[i].a_y * period;
	}

	/* compute the new position */
	for (int i = 0; i < g_N; ++i) {
		g_pv[i].pos_x += g_pv[i].v_x * period;
		g_pv[i].pos_y += g_pv[i].v_y * period;
	}
}