#include <time.h>
#include <stdio.h>
#include <float.h>
#include <omp.h>

#include "params.c"
#include "utils.c"
#include "cam_utils.c"
#include "object_utils.c"
#include "numerical_integration.c"
#include "intersection.c"

void make_image(img_t img, object_list_t objects, char* image_name) {
	int i, j, k;
	for (i = 0; i < img.height; i++)
		for (j = 0; j < img.width; j++)
			for (k = 0; k < 3; k++)
				img.img[i][j][k] = 0;

	int m = img.super_sample * img.height;
	int n = img.super_sample * img.width;
	double start, end;
	start = omp_get_wtime();
	for (i = 0; i < m; i++) {
		// Print progress bar between each line of rays
		if (progress_bar) {
			end = omp_get_wtime();
			print_progress(i / ((double) m), end - start);
		}

		#pragma omp parallel for
		for (j = 0; j < n; j++) {
			int* img_coords = calloc(2, sizeof(double));
			double* x1 = calloc(4, sizeof(double));
			double* y1 = calloc(4, sizeof(double));
			double* x2 = calloc(4, sizeof(double));
			double* y2 = calloc(4, sizeof(double));
			double* ch = calloc(4, sizeof(double));
			double* k1 = calloc(4, sizeof(double));
			double* k2 = calloc(4, sizeof(double));
			double* k3 = calloc(4, sizeof(double));
			double* k3_alt = calloc(4, sizeof(double));
			double* k4 = calloc(4, sizeof(double));
			intersection_t intersection;
			intersection.col = calloc(3, sizeof(double));
			double max_time;

			if (ren == relativistic) {
				max_time = 1.0;
			} else if (ren == euclidean) {
				max_time = DBL_MAX;
			} else if (ren == mix) {
				if (m - i > j * ((double) img.height) / img.width) {
					max_time = 1.0;
				} else {
					max_time = DBL_MAX;
				}
			}

			generate_ray(cam, img, i, j, img_coords, x1, y1);
			trace_adaptive_ray(ray_params, &intersection, objects, x1, y1, x2, y2, ch, k1, k2, k3, k3_alt, k4, max_time);

			if (intersection.collided) {
				add_col(img_coords, img, intersection.col);
			}

			free(img_coords);
			free(x1);
			free(x2);
			free(y1);
			free(y2);
			free(ch);
			free(k1);
			free(k2);
			free(k3);
			free(k3_alt);
			free(k4);
			free(intersection.col);
		}
	}
	// Print final progress bar
	if (progress_bar) {
		end = omp_get_wtime();
		print_progress(1.0, end - start);
		printf("\n");
	}

	// Normalize pixels by amount of rays per pixel
	for (i = 0; i < img.height; i++)
		for (j = 0; j < img.width; j++)
			for (k = 0; k < 3; k++)
				img.img[i][j][k] = img.img[i][j][k] / (img.super_sample * img.super_sample);

	upsample_image(img);

	// Save image
//	save_image(img, image_name);
	save_png(img, image_name);
}
