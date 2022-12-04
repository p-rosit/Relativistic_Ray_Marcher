#ifndef OBJECT_UTILS
#define OBJECT_UTILS
#include <math.h>
#include "utils.c"

void prepare_position(pos_t* p) {
	// Normalize normal vector
        double norm = sqrt(p->n[0]*p->n[0] + p->n[1]*p->n[1] + p->n[2]*p->n[2]);
        p->n[0] /= norm;
        p->n[1] /= norm;
        p->n[2] /= norm;
        p->d = -(p->t[0]*p->n[0] + p->t[1]*p->n[1] + p->t[2]*p->n[2]);

	// Make sure tangent is actually normal to normal
	norm = p->tangent[0]*p->n[0] + p->tangent[1]*p->n[1] + p->tangent[2]*p->n[2];
	p->tangent[0] -= norm * p->n[0];
	p->tangent[1] -= norm * p->n[1];
	p->tangent[2] -= norm * p->n[2];

	// Normalize tangent
        norm = sqrt(p->tangent[0]*p->tangent[0] + p->tangent[1]*p->tangent[1] + p->tangent[2]*p->tangent[2]);
        p->tangent[0] /= norm;
        p->tangent[1] /= norm;
        p->tangent[2] /= norm;

	// Compute binormal
        p->binormal[0] = p->n[1]*p->tangent[2] - p->n[2]*p->tangent[1];
        p->binormal[1] = p->n[2]*p->tangent[0] - p->n[0]*p->tangent[2];
        p->binormal[2] = p->n[0]*p->tangent[1] - p->n[1]*p->tangent[0];
}

void prepare_sphere(sphere_t* s, char* file_name) {
	prepare_position(&s->p);
        s->sq_r = (s->r) * (s->r);
	load_texture(&s->tex, file_name);
}

void prepare_disk(disk_t* d, char* file_name) {
	prepare_position(&d->p);
        d->sq_r1 = (d->r1) * (d->r1);
        d->sq_r2 = (d->r2) * (d->r2);
        load_texture(&d->tex, file_name);
}

void prepare_plane(plane_t* p, char* file_name) {
        prepare_position(&p->p);
	load_texture(&p->tex, file_name);
}

void grab_texture(intersection_t* intersection, texture_t tex, double x, double y) {
	double  u = x * tex.m;
	double  v = y * tex.n;
	int u1 = floor(u);
	int u2 = ceil(u);
	int v1 = floor(v);
	int v2 = ceil(v);

	if (u1 < 0) {u1 = 0; u = 0; u2 = 0;}
	if (u2 >= tex.m) {u1 = tex.m - 1; u = tex.m - 1; u2 = tex.m - 1;}
	if (v1 < 0) {v1 = 0; v = 0; v2 = 0;}
	if (v2 >= tex.n) {v1 = tex.n - 1; v = tex.n - 1; v2 = tex.n - 1;}

	for (int i = 0; i < 3; i++)
		intersection->col[i] = ((u2 - u ) * (v2 - v ) * tex.tex[u1][v1][i] +
			  		(u  - u1) * (v2 - v ) * tex.tex[u2][v1][i] +
			  		(u2 - u ) * (v  - v1) * tex.tex[u1][v2][i] +
			  		(u  - u1) * (v  - v1) * tex.tex[u2][v2][i]);
}

void add_col(int img_coords[2], img_t img, double col[3]) {
        img.img[img_coords[0]][img_coords[1]][0] += col[0];
        img.img[img_coords[0]][img_coords[1]][1] += col[1];
        img.img[img_coords[0]][img_coords[1]][2] += col[2];
}
#endif
