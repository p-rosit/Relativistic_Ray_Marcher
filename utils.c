#ifndef UTILS
#define UTILS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <png.h>

#define RGB_COMPONENT_COLOUR (255)
#define PBSTR   ("###################################################")
#define PBWIDTH (50)

double pi = 3.1415918281828459045;
typedef int bool;
const bool true = 1;
const bool false = 0;

typedef struct img_t		img_t;
typedef struct cam_t		cam_t;
typedef struct params_t		params_t;

typedef struct pos_t		pos_t;
typedef struct intersection_t 	intersection_t;
typedef struct texture_t	texture_t;
typedef struct sphere_t		sphere_t;
typedef struct disk_t		disk_t;
typedef struct triangle_t	triangle_t;
typedef struct plane_t		plane_t;

typedef struct object_list_t object_list_t;

enum render {relativistic, euclidean, mix};
enum texture_choice {debug, real};

struct img_t {
	double*** img;
	double*** upsample_img;
	int height;
	int width;
	int super_sample;
	int upsample;
};

struct cam_t {
	double t[3];
	double R[3][3];
	double f;
	double a;
	double cx;
	double cy;
};

struct params_t {
	double ray_time;
	double min_tol;
	double tol;
	double max_tol;
	double min_dt;
	double max_dt;
	double max_trials;
};

struct pos_t {
        double t[3];
        double n[3];
        double tangent[3];
        double binormal[3];
        double d;
};

struct intersection_t {
	double t;
	int collided;
	double* col;
};

struct texture_t {
        int m;
        int n;
        double*** tex;
};

struct material_t {
	double reflectivity;
	double refraction_index;
	double transparency;
};

enum object {sphere, disk, plane};

struct sphere_t {
	pos_t p;
        double r;
	double sq_r;
	texture_t tex;
};

struct disk_t {
	pos_t p;
        double r1;
        double r2;
	double sq_r1;
	double sq_r2;
        texture_t tex;
};

struct plane_t {
        pos_t p;
        double L;
        texture_t tex;
};

struct object_list_t {
	int n;
	int* id;
	void** obj;
};

void make_img(img_t* img, int height, int width, int super_sample, int upsample) {
	int i, j, k;

	img->img = calloc(height, sizeof(double**));
	for (i = 0; i < height; i++) {
		img->img[i] = calloc(width, sizeof(double*));
		for (j = 0; j < width; j++) {
			img->img[i][j] = calloc(3, sizeof(double));
		}
	}
	img->upsample_img = calloc(upsample * height, sizeof(double**));
	for (i = 0; i < upsample * height; i++) {
		img->upsample_img[i] = calloc(upsample * width, sizeof(double*));
		for (j = 0; j < upsample * width; j++) {
			img->upsample_img[i][j] = calloc(3, sizeof(double));
		}
	}

	img->height = height;
	img->width = width;
	img->super_sample = super_sample;
	img->upsample = upsample;
}

void free_img(img_t img) {
	int i, j;
	for (i = 0; i < img.height; i++) {
		for (j = 0; j < img.width; j++) {
			free(img.img[i][j]);
		}
		free(img.img[i]);
	}
	free(img.img);
	for (i = 0; i < img.upsample * img.height; i++) {
		for (j = 0; j < img.upsample * img.width; j++) {
			free(img.upsample_img[i][j]);
		}
		free(img.upsample_img[i]);
	}
	free(img.upsample_img);
}

double max(double a, double b) {
	if (a >= b) {
		return a;
	} else {
		return b;
	}
}

double min(double a, double b) {
	if (a <= b) {
		return a;
	} else {
		return b;
	}
}

void make_list(int n_objects, object_list_t* list) {
	list->n = 0;
	list->id = calloc(n_objects, sizeof(int));
	list->obj = calloc(n_objects, sizeof(void*));
}

void add_sphere(object_list_t *list, sphere_t* s) {
	list->id[list->n] = sphere;
	list->obj[list->n] = s;
	list->n += 1;
}

void add_disk(object_list_t *list, disk_t* d) {
	list->id[list->n] = disk;
	list->obj[list->n] = d;
	list->n += 1;
}

void add_plane(object_list_t *list, plane_t* p) {
	list->id[list->n] = plane;
	list->obj[list->n] = p;
	list->n += 1;
}

void free_list(object_list_t list) {
	free(list.id);
	free(list.obj);
}

void print_progress(double percentage, int time) {
	int val = (int) (100 * percentage);
	int lpad = (int) (percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;

	printf("\r%3d%% [%.*s%*s] Time taken: %ds", val, lpad, PBSTR, rpad, "", time);
	fflush(stdout);
}

void upsample_image(img_t img) {
	int i, j, k, l, m;
	double x, y;
	int x1, x2, y1, y2;
	double dp = 1.0 / img.upsample;

	if (img.upsample > 1) {
		for (i = 0; i < img.height; i++) {
			for (j = 0; j < img.width; j++) {
				for (k = 0; k < img.upsample; k++) {
					for (l = 0; l < img.upsample; l++) {
						x = i + k * dp + 0.5 * dp - 0.5;
						y = j + l * dp + 0.5 * dp - 0.5;
						x1 = floor(x);
						x2 = ceil(x);
						y1 = floor(y);
						y2 = ceil(y);

						if (x1 < 0) {x1 = 0; x = 0;}
						if (x2 >= img.height - 1) {x2 = img.height - 1; x = img.height - 1;}
						if (y1 < 0) {y1 = 0; y = 0;}
						if (y2 >= img.width - 1) {y2 = img.width - 1; y = img.width - 1;}

						for (m = 0; m < 3; m++) {
							img.upsample_img[i*img.upsample + k][j*img.upsample + l][m] = (	(x2 - x ) * (y2 - y ) * img.img[x1][y1][m] +
															(x  - x1) * (y2 - y ) * img.img[x2][y1][m] +
															(x2 - x ) * (y  - y1) * img.img[x1][y2][m] +
															(x  - x1) * (y  - y1) * img.img[x2][y2][m]);
						}
					}
				}
			}
		}
	} else {
		for (i = 0; i < img.height; i++)
			for (j = 0; j < img.width; j++)
				for (k = 0; k < 3; k++)
					img.upsample_img[i][j][k] = img.img[i][j][k];
	}
}

void save_image(img_t img, char* image_name) {
        FILE *f = fopen(image_name, "wb");
        fprintf(f, "P6\n%i %i 255\n", img.upsample * img.width, img.upsample * img.height);

        for (int i = 0; i < img.upsample * img.height; i++) {
            for (int j = 0; j < img.upsample * img.width; j++) {
                fputc((int) 255 * img.upsample_img[i][j][0], f);
                fputc((int) 255 * img.upsample_img[i][j][1], f);
                fputc((int) 255 * img.upsample_img[i][j][2], f);
            }
        }
        fclose(f);
}

void save_png(img_t img, char* image_name) {
	int y;
	FILE *fp = fopen(image_name, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open file '%s'\n", image_name);
		exit(1);
	}

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (!png) {
		fprintf(stderr, "Unable to make png write struct for '%s'\n", image_name);
		exit(1);
	}

	png_infop info = png_create_info_struct(png);
	if (!info) {
		fprintf(stderr, "Unable to make png info struct for '%s'\n", image_name);
		exit(1);
	}

	if (setjmp(png_jmpbuf(png))) {
		fprintf(stderr, "No clue: '%s'\n", image_name);
		exit(1);
	}

	png_init_io(png, fp);

	int m = img.upsample * img.height;
	int n = img.upsample * img.width;

	png_set_IHDR(
		png,
		info,
		n,
		m,
		8,
		PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT
	);
	png_write_info(png, info);

	int i, j, k;
	unsigned char** png_img = calloc(m, sizeof(unsigned char*));
	for (i = 0; i < m; i++) {
		png_img[i] = calloc(n * 3, sizeof(unsigned char));
		for (j = 0; j < n; j++) {
			for (k = 0; k < 3; k++) {
				png_img[i][3 * j + k] = 255 * img.upsample_img[i][j][k];
			}
		}
	}

	png_write_image(png, png_img);
	png_write_end(png, NULL);

	fclose(fp);
	png_destroy_write_struct(&png, &info);

	for (i = 0; i < m; i++)
		free(png_img[i]);
	free(png_img);
}

void free_texture(texture_t texture) {
	for (int i = 0; i < texture.m; i++) {
		for (int j = 0; j < texture.n; j++)
			free(texture.tex[i][j]);
		free(texture.tex[i]);
	}
	free(texture.tex);
}

void load_texture(texture_t* texture, const char* file_name) {
	char buff[16];
	int c, rgb_comp_colour, pix;

	// Open PPM file for reading
        FILE *f = fopen(file_name, "rb");
        if (!f) {
		fprintf(stderr, "Unable to open file '%s'\n", file_name);
		exit(1);
        }

        // Read image format
        if (!fgets(buff, sizeof(buff), f)) {
        	perror(file_name);
        	exit(1);
        }

	//check the image format
	if (buff[0] != 'P' || buff[1] != '6') {
		fprintf(stderr, "Invalid image format (must be 'P6')\n");
		exit(1);
	}

	//check for comments
	c = getc(f);
	while (c == '#') {
		while (getc(f) != '\n') ;
		c = getc(f);
	}

	ungetc(c, f);
	//read image size information
	if (fscanf(f, "%d %d", &texture->m, &texture->n) != 2) {
		fprintf(stderr, "Invalid image size (error loading '%s')\n", file_name);
		exit(1);
	}

	//read rgb component
	if (fscanf(f, "%d", &rgb_comp_colour) != 1) {
		fprintf(stderr, "Invalid rgb component (error loading '%s')\n", file_name);
		exit(1);
	}

	//check rgb component depth
	if (rgb_comp_colour != RGB_COMPONENT_COLOUR) {
		fprintf(stderr, "'%s' does not have 8-bits components\n", file_name);
		exit(1);
	}

	while (fgetc(f) != '\n') ;

	//memory allocation for pixel data
	texture->tex = calloc(texture->m, sizeof(double**));
	for (int i = 0; i < texture->m; i++) {
		texture->tex[i] = calloc(texture->n, sizeof(double*));
		for (int j = 0; j < texture->n; j++)
			texture->tex[i][j] = calloc(3, sizeof(double));
	}
	unsigned char* data = calloc(texture->m * texture->n * 3, sizeof(unsigned char));

	if (!data) {
		fprintf(stderr, "Unable to allocate memory\n");
		exit(1);
	}

	//read pixel data from file
	if (fread(data, 3 * texture->n, texture->m, f) != texture->n) {
		fprintf(stderr, "Error loading image '%s'\n", file_name);
		exit(1);
	}

	//read pixel data from file
	for (int i = 0; i < texture->m; i++) {
		for (int j = 0; j < texture->m; j++) {
			for (int k = 0; k < 3; k++) {
				texture->tex[i][j][k] = ((double) data[i * texture->n * 3 + j * 3 + k]) / 255.0;
			}
		}
	}

	free(data);
	fclose(f);
}


#endif
