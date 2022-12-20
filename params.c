#ifndef PARAMS
#define PARAMS
#include "utils.c"
#include "cam_setup.c"
// gcc -fopenmp -O3 ray_tracing.c -lm -lpng
// convert -delay 5 -loop 0 *.png name.gif

// What to render

// 4k, supersampling 5x5
// stars	1460s
// debug	1456s
// euclidean s	63s
// euclidean d	56s

// 1080p, normal view, supersampling 5x5
// stars	580s
// debug	587s
// euclidean r	13s
// euclidean d	10s
// no disk rel	555s
// no disk euc	13s
// no disk s r	794
// no disk s e	12s

// 540p, rotation in xy-plane, supersampling 3x3
// stars	250 x 52s
// debug	250 x 52s
// euclidean s	250 x 2s
// euclidean d	250 x 2s
// mix stars
// mix debug	250 x 37s
// no disk rel	250 * 52s
// no disk euc  250 * 1s

// 720p, rotation in xz-plane, supersampling 3x3
// real
// debug	250 x 51s
// euclidean s	250 x 1s

// Render parameters

const int width		= 960;//1920;//1280; 960
const int height	= 540;//1080;//720;  540
const int super_sample	= 3;
const int upsample	= 1;

const int ren		= relativistic;
const int tex		= debug;
const bool progress_bar = true;

const bool include_disk = true;

const bool single_img = true;
const bool xy_rotation = false;
const char* image_name = "img";
const int n_image = 250;

// Camera parameters

cam_t cam = {	.t	= {-4, 0.0, 0.3},
		.f	= 1 / ((double) 2 * width),
		.a	= 0.3,
		.cx	= (double) width / 2,
		.cy	= (double) height / 2};

double cam_look_at[3] = {0, 0, 0};//{0, 0.5, 0.5};
double angle = 0;

// Ray tracing parameters

params_t ray_params = {.ray_time	= 100.0,
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
