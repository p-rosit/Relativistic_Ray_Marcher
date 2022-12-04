#include <math.h>
#include "utils.c"
#include "params.c"

void cam_coords2world_coords(cam_t cam, double *x) {
	double a, b, c, d;
	c = x[0];
	d = x[1];

	a = (c - cam.cy) * cam.f;
	b = (d - cam.cx) * cam.f;
	c = cam.a;

	x[0] = 0;
	x[1] = cam.R[0][0] * a + cam.R[0][1] * b + cam.R[0][2] * c + cam.t[0];
	x[2] = cam.R[1][0] * a + cam.R[1][1] * b + cam.R[1][2] * c + cam.t[1];
	x[3] = cam.R[2][0] * a + cam.R[2][1] * b + cam.R[2][2] * c + cam.t[2];
}

void world_coords2ray_dirs(cam_t cam, double *x, double *y) {
	double a, b, c, r;

	a = x[1] - cam.t[0];
	b = x[2] - cam.t[1];
	c = x[3] - cam.t[2];
	r = 1 / sqrt(a*a + b*b + c*c);

	y[0] = 1;
	y[1] = a * r;
	y[2] = b * r;
	y[3] = c * r;
}

void generate_ray(cam_t cam, img_t img, int i, int j, int *img_coords, double *x, double *y) {
	double dx = img.height / ((double) img.super_sample * img.height);
	double dy = img.width / ((double) img.super_sample * img.width);

	double p = dx * i + dx / 2;
	double q = dy * j + dy / 2;
	img_coords[0] = floor(p);
	img_coords[1] = floor(q);

	x[0] = p;
	x[1] = q;

	cam_coords2world_coords(cam, x);
	world_coords2ray_dirs(cam, x, y);
}
