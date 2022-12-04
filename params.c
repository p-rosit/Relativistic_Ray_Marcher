#ifndef PARAMS
#define PARAMS
#include "utils.c"
#include "cam_setup.c"
// gcc -fopenmp -O3 ray_tracing -lm -lpng

// What to render

// 4k, supersampling 5x5
// real		1460s
// debug	1456s
// euclidean r	63s
// euclidean d	56s

// 1080p, normal view, supersampling 5x5
// real
// debug
// euclidean r
// euclidean d

// 720p, rotation in xy-plane, supersampling 3x3
// real		200 x 93s
// debug
// mix real

// 720p, rotation in xz-plane, supersampling 3x3
// real
// debug


// Render parameters

const int width		= 1280;
const int height	= 720;
const int super_sample	= 3;
const int upsample	= 1;

const int ren		= relativistic;
const int tex		= real;
const bool progress_bar = true;

const bool include_disk = true;

const bool single_img = false;
const bool xy_rotation = true;
const char* image_name = "img";
const int n_image = 200;

// Camera parameters

cam_t cam = {	.t	= {-2, 0.0, 0.1},
		.f	= 1 / ((double) 2 * width),
		.a	= 0.3,
		.cx	= (double) width / 2,
		.cy	= (double) height / 2};

double cam_look_at[3] = {0, 0, 0};//{0, 0.5, 0.5};
double angle = -0.25;

// Ray tracing parameters

params_t ray_params = {	.ray_time	= 100.0,
			.min_tol	= 1e-3,
			.tol		= 1e-2,
			.max_tol	= 1e-1,
			.min_dt		= 1e-6,
			.max_dt		= 10,
			.max_trials	= 100};

// Parameters for object in scene

const double RS = 1.0;
const double EVENT_HORIZON = RS / 4.0;

const double sphere_radius = EVENT_HORIZON;
const double disk_inner_radius = 3 * EVENT_HORIZON;
const double disk_outer_radius = 6 * EVENT_HORIZON;

const double plane_length = 7;
#endif
