/* sequential version of N body simulation */

# include <iostream>
# include <cmath>

using namespace std;

double G = 6.67 * pow(10, -11);
//# define G	6.67
double e = 0.00001;
//# define e 	1
# define period	100.0

void computeForce();
void compute_a();
void computerVelocity();
void computePosition();

float Force[2][2] = {{0, 0}, {0, 0}};

/* 0: m; (1, 2): (px, py); (3, 4): (vx, vy); (5, 6): (ax, ay) */
float particle[2][7] = {{5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5}, {2.2, 1.5, 1.0, 0.0, 0.0, 0.0, 0.0}};
// double particle[2][7] = {{1.0, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0}, {1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0}};

int main(int argc, char ** argv) {
	int time = 0;
	while (time < 10) {
		cout << Force[0][0] << "(" << particle[0][1] << "," << particle[0][2] << ")" << endl;
		cout << Force[1][0] << "(" << particle[1][1] << "," << particle[1][2] << ")" << endl;
		
		computeForce();
		compute_a();
		computerVelocity();
		computePosition();
		
		++time;
	}
	return 0;
}

/* compute the force */
void computeForce() {
	int i = 0;
	int j = 0;
	for (i = 0; i < 2; ++i) {
		for (j = 0; j < 2; ++j) {
			if (i == j)
				continue;
			else {
				double r_2 = pow((particle[i][1] - particle[j][1]), 2) + pow((particle[i][2] - particle[j][2]), 2);
				Force[i][0] += (-1) * G * particle[i][0] * particle[j][0] * (particle[i][1] - particle[j][1]) / (pow(r_2 + e, 1.5));	
				Force[i][1] += (-1) * G * particle[i][0] * particle[j][0] * (particle[i][2] - particle[j][2]) / (pow(r_2 + e, 1.5));
			}
		} 
	}
}

/* using force to compute acceleration */
void compute_a() {
	int i;
	for (i = 0; i < 2; ++i) {
		particle[i][5] = 1.0 * Force[i][0] / particle[i][0];
		particle[i][6] = 1.0 * Force[i][1] / particle[i][0];
	}
}

/* compute the velocity */
void computerVelocity() {
	int i;
	for (i = 0; i < 2; ++i) {
		particle[i][3] += 1.0 * particle[i][5] * period; 
		particle[i][4] += 1.0 * particle[i][6] * period;
	}
}

/* compute the new positon */
void computePosition() {
	int i;
	for (i = 0; i < 2; ++i) {
		particle[i][1] += particle[i][3] * period;
		cout << particle[i][1] << endl;
		particle[i][2] += particle[i][4] * period;
		cout << particle[i][2] << endl;
	}
}
