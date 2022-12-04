#include <stdio.h>
#include <math.h>
#include "utils.c"

void angle_axis2mat(double r[3][3], double axis[3], double angle) {
	double c = cos(angle);
	double s = sin(angle);

	r[0][0] = c + axis[0] * axis[0] * (1 - c);
	r[0][1] =     axis[0] * axis[1] * (1 - c) - axis[2] * s;
	r[0][2] =     axis[0] * axis[2] * (1 - c) + axis[1] * s;
	r[1][0] =     axis[1] * axis[0] * (1 - c) + axis[2] * s;
	r[1][1] = c + axis[1] * axis[1] * (1 - c);
	r[1][2] =     axis[1] * axis[2] * (1 - c) - axis[0] * s;
	r[2][0] =     axis[2] * axis[0] * (1 - c) - axis[1] * s;
	r[2][1] =     axis[2] * axis[1] * (1 - c) + axis[0] * s;
	r[2][2] = c + axis[2] * axis[2] * (1 - c);
}

void point_camera_at(cam_t* cam, double point[3], double angle, int direction) {
	int i;
	double x[3], y[3], z[3];
	double normx, normy, normz;
	double rot[3][3];

	// Compute z vector pointing from camera translation to the look-at-point
	z[0] = point[0] - cam->t[0];
	z[1] = point[1] - cam->t[1];
	z[2] = point[2] - cam->t[2];
	normz = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
	z[0] /= normz; z[1] /= normz; z[2] /= normz;

	// Construct matrix that rotates vectors around z by angle
	angle_axis2mat(rot, z, angle);

	// Compute potential x vectors in the xy-plane normal to z vector
	x[2] = 0;
	if (fabs(z[0]) < fabs(z[1])) {
		x[0] = direction;
		x[1] = - direction * z[0] / z[1];
	} else {
		x[0] = - direction * z[1] / z[0];
		x[1] = direction;
	}
	normx = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	x[0] /= normx; x[1] /= normx;

	// Rotate potential x vectors around z by angle
	x[0] = rot[0][0]*x[0] + rot[0][1]*x[1] + rot[0][2]*x[2];
	x[1] = rot[1][0]*x[0] + rot[1][1]*x[1] + rot[1][2]*x[2];
	x[2] = rot[2][0]*x[0] + rot[2][1]*x[1] + rot[2][2]*x[2];

	// Compute y vector normal to determined x and z vectors
	y[0] = z[1] * x[2] - z[2] * x[1];
	y[1] = z[2] * x[0] - z[0] * x[2];
	y[2] = z[0] * x[1] - z[1] * x[0];

	// Set rotation matrix
	cam->R[0][0] = -y[0];
	cam->R[1][0] = -y[1];
	cam->R[2][0] = -y[2];
	cam->R[0][1] =  x[0];
	cam->R[1][1] =  x[1];
	cam->R[2][1] =  x[2];
	cam->R[0][2] =  z[0];
	cam->R[1][2] =  z[1];
	cam->R[2][2] =  z[2];
}

void rotate_camera_around_axis(cam_t* cam, double axis[3], double angle) {
	int i, j, k;
	double r[3][3];
	double res[3][3];
	angle_axis2mat(r, axis, angle);

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			res[i][j] = 0;
			for (k = 0; k < 3; k++) {
				res[i][j] += r[i][k] * cam->R[k][j];
			}
		}
	}

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			cam->R[i][j] = res[i][j];
}
