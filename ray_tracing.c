#include <time.h>
#include <stdio.h>
#include <float.h>

#include "params.c"
#include "image_generation.c"

void main() {
	// Construct objects in scene
	sphere_t sphere = {	.p = {	.t = {0.0, 0.0, 0.0},
					.n = {0.0, 0.0, 1.0},
					.tangent = {0.0, 1.0, 0.0}},
				.r = sphere_radius};

	disk_t disk = {	.p = {	.t = {0.0, 0.0, 0.0},
				.n = {0.0, 0.0, 1.0},
				.tangent = {0.0, 1.0, 0.0}},
			.r1 = disk_inner_radius,
        	       	.r2 = disk_outer_radius};

	plane_t plane1 = {	.p = {	.t = {plane_length, 0.0, 0.0},
	                        	.n = {1.0, 0.0, 0.0},
	                        	.tangent = {0.0, 0.0, 1.0}},
	                  	.L = 1.1 * plane_length};
	plane_t plane2 = {	.p = {	.t = {-plane_length, 0.0, 0.0},
	                        	.n = {-1.0, 0.0, 0.0},
	                        	.tangent = {0.0, 0.0, -1.0}},
	                  	.L = 1.1 * plane_length};
	plane_t plane3 = {	.p = {	.t = {0.0, plane_length, 0.0},
	                        	.n = {0.0, 1.0, 0.0},
	                        	.tangent = {0.0, 0.0, 1.0}},
	                  	.L = 1.1 * plane_length};
	plane_t plane4 = {	.p = {	.t = {0.0, -plane_length, 0.0},
	                        	.n = {0.0, -1.0, 0.0},
	                        	.tangent = {0.0, 0.0, 1.0}},
	                  	.L = 1.1 * plane_length};
	plane_t plane5 = {	.p = {	.t = {0.0, 0.0, plane_length},
	                        	.n = {0.0, 0.0, 1.0},
	                        	.tangent = {1.0, 0.0, 0.0}},
	                  	.L = 1.1 * plane_length};
	plane_t plane6 = {	.p = {	.t = {0.0, 0.0, -plane_length},
	                        	.n = {0.0, 0.0, -1.0},
	                        	.tangent = {1.0, 0.0, 0.0}},
	                  	.L = 1.1 * plane_length};

	// Prepare objects in scene
	if (tex == real) {
		prepare_sphere(&sphere, "textures/black_hole.ppm");
		prepare_disk(&disk, "textures/disk.ppm");
		prepare_plane(&plane1, "textures/stars.ppm");
		prepare_plane(&plane2, "textures/stars.ppm");
		prepare_plane(&plane3, "textures/stars.ppm");
		prepare_plane(&plane4, "textures/stars.ppm");
		prepare_plane(&plane5, "textures/stars.ppm");
		prepare_plane(&plane6, "textures/stars.ppm");
	} else if (tex == debug) {
		prepare_sphere(&sphere, "textures/sphere_debug.ppm");
		prepare_disk(&disk, "textures/disk_debug.ppm");
		prepare_plane(&plane1, "textures/plane_debug1.ppm");
		prepare_plane(&plane2, "textures/plane_debug2.ppm");
		prepare_plane(&plane3, "textures/plane_debug3.ppm");
		prepare_plane(&plane4, "textures/plane_debug4.ppm");
		prepare_plane(&plane5, "textures/plane_debug5.ppm");
		prepare_plane(&plane6, "textures/plane_debug6.ppm");
	} else {
		printf("Invalid texture loading choice. Unreachable error.\n");
	}

	// Add objects to object list
	object_list_t objects;
	make_list(8, &objects);
	add_sphere(&objects, &sphere);
	if (include_disk) {
		add_disk(  &objects, &disk);
	}
	add_plane( &objects, &plane1);
	add_plane( &objects, &plane2);
	add_plane( &objects, &plane3);
	add_plane( &objects, &plane4);
	add_plane( &objects, &plane5);
	add_plane( &objects, &plane6);

	img_t img;
	make_img(&img, height, width, super_sample, upsample);

	char img_name[60];
	if (single_img) {
		// Render single image
		point_camera_at(&cam, cam_look_at, angle, 1);
		sprintf(img_name, "%s.png", image_name);
		make_image(img, objects, img_name);
	} else {
		// Render series of images
		// Make camera look at some point
		double axis[3];
		if (xy_rotation) {
			axis[0] = 0.0;
			axis[1] = 0.0;
			axis[2] = 1.0;
		} else {
			axis[0] = 0.0;
			axis[1] = 1.0;
			axis[2] = 0.0;
		}
		double cam_radius = 4.5;

		double start_angle = 0;
		double end_angle = 2 * pi;
		double da = (end_angle - start_angle) / n_image;

		cam.t[0] = cam_radius;
		cam.t[1] = 0;
		if (xy_rotation) {
			cam.t[2] = -0.3;
		} else {
			cam.t[2] = 0;
		}
		point_camera_at(&cam, cam_look_at, 0, 1);
		for (int i = 0; i < n_image; i++) {
			if (xy_rotation) {
				cam.t[0] = cam_radius * cos(i * da);
				cam.t[1] = cam_radius * sin(i * da);
				cam.t[2] = -0.3;
			} else {
				cam.t[0] = cam_radius * cos(i * da);
				cam.t[1] = 0;
				cam.t[2] = cam_radius * sin(-i * da);
			}

			sprintf(img_name, "%s%03d.png", image_name, i);
			printf("Rendering image: %03d\n", i);
			make_image(img, objects, img_name);
			rotate_camera_around_axis(&cam, axis, da);
		}
	}

	free_img(img);

	free(objects.id);
	free(objects.obj);

	free_texture(sphere.tex);
	free_texture(disk.tex);
	free_texture(plane1.tex);
	free_texture(plane2.tex);
	free_texture(plane3.tex);
	free_texture(plane4.tex);
	free_texture(plane5.tex);
	free_texture(plane6.tex);
}
