#ifndef INTERSECTION
#define INTERSECTION
#include "utils.c"
#include "object_utils.c"

void sphere_collision(intersection_t* intersection, sphere_t s, double *x1, double *x2, double max_time) {
	double a, b, c, d, e, t;
	double f[3], y[3];
	for (int i = 0; i < 3; i++) {
		f[i] = x1[i+1] - s.p.t[i];
		y[i] = x2[i+1] - x1[i+1];
	}

	a = y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
	b = 2 * (y[0]*f[0] + y[1]*f[1] + y[2]*f[2]);
	c = f[0]*f[0] + f[1]*f[1] + f[2]*f[2] - s.sq_r;
	d = b / (2 * a);
	e = d * d - c / a;

	if (e < 0) {
		return;
	}

	t = -d - sqrt(e);
	if (t < intersection->t && 0 < t && t < max_time) {
		y[0] = x1[1] + t * (x2[1] - x1[1]) - s.p.t[0];
		y[1] = x1[2] + t * (x2[2] - x1[2]) - s.p.t[1];
		y[2] = x1[3] + t * (x2[3] - x1[3]) - s.p.t[2];
		f[0] = y[0]*s.p.tangent[0] + 	y[1]*s.p.tangent[1] + 	y[2]*s.p.tangent[2];
		f[1] = y[0]*s.p.binormal[0] + 	y[1]*s.p.binormal[1] + 	y[2]*s.p.binormal[2];
		f[2] = y[0]*s.p.n[0] + 		y[1]*s.p.n[1] +		y[2]*s.p.n[2];

		if (fabs(f[0]) < 1e-16) {
			a = 0;
		} else {
			a = atan2(f[1], f[0]);
		}

		if (fabs(f[2]) < 1e-16) {
			b = 0;
		} else {
			b = atan2(sqrt(f[0]*f[0] + f[1]*f[1]), f[2]);
		}

		a = (a + pi) / (2 * pi);
		b = b / pi;

		intersection->t = t;
		intersection->collided = true;
		grab_texture(intersection, s.tex, a, b);
	}
}

void disk_collision(intersection_t* intersection, disk_t d, double *x1, double *x2, double max_time) {
	double x, y, r, t;
	double a[3];

	t = -(d.p.d + d.p.n[0]*x1[1] + d.p.n[1]*x1[2] + d.p.n[2]*x1[3]) / (d.p.n[0]*(x2[1] - x1[1]) + d.p.n[1]*(x2[2] - x1[2]) + d.p.n[2]*(x2[3] - x1[3]));

	if (t < intersection->t && 0 < t && t < max_time) {
		a[0] = x1[1] + t * (x2[1] - x1[1]) - d.p.t[0];
		a[1] = x1[2] + t * (x2[2] - x1[2]) - d.p.t[1];
		a[2] = x1[3] + t * (x2[3] - x1[3]) - d.p.t[2];

		x = a[0] * d.p.tangent[0] + a[1] * d.p.tangent[1] + a[2] * d.p.tangent[2];
		y = a[0] * d.p.binormal[0] + a[1] * d.p.binormal[1] + a[2] * d.p.binormal[2];
		r = x * x + y * y;

		if (d.sq_r1 <= r && r < d.sq_r2) {
			x = (x + d.r2) / (2 * d.r2);
			y = (y + d.r2) / (2 * d.r2);

			intersection->t = t;
			intersection->collided = true;
			grab_texture(intersection, d.tex, x, y);
		}
	}
}

void plane_collision(intersection_t* intersection, plane_t p, double *x1, double *x2, double max_time) {
	double x, y, t;
	double a[3];

	t = -(p.p.d + p.p.n[0]*x1[1] + p.p.n[1]*x1[2] + p.p.n[2]*x1[3]) / (p.p.n[0]*(x2[1] - x1[1]) + p.p.n[1]*(x2[2] - x1[2]) + p.p.n[2]*(x2[3] - x1[3]));

	if (t < intersection->t && 0 < t && t < max_time) {
		a[0] = x1[1] + t * (x2[1] - x1[1]) - p.p.t[0];
		a[1] = x1[2] + t * (x2[2] - x1[2]) - p.p.t[1];
		a[2] = x1[3] + t * (x2[3] - x1[3]) - p.p.t[2];

		x = a[0] * p.p.tangent[0] + a[1] * p.p.tangent[1] + a[2] * p.p.tangent[2];
		y = a[0] * p.p.binormal[0] + a[1] * p.p.binormal[1] + a[2] * p.p.binormal[2];

		if (-p.L < x && x < p.L && -p.L < y && y < p.L) {
			x = (x + p.L) / (2 * p.L);
			y = (y + p.L) / (2 * p.L);

			intersection->t = t;
			intersection->collided = true;
			grab_texture(intersection, p.tex, x, y);
		}
	}
}

void check_collision(intersection_t* intersection, object_list_t objects, double* x1, double* x2, double max_time) {
	for (int i = 0; i < objects.n; i++) {
		if (objects.id[i] == sphere) {
			sphere_collision(intersection, *((sphere_t*) objects.obj[i]), x1, x2, max_time);
		} else if (objects.id[i] == disk) {
			disk_collision(intersection, *((disk_t*) objects.obj[i]), x1, x2, max_time);
		} else if (objects.id[i] == plane) {
			plane_collision(intersection, *((plane_t*) objects.obj[i]), x1, x2, max_time);
		}
	}
}
#endif
