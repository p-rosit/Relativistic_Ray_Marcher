#include <math.h>
#include "params.c"
#include "utils.c"
#include "intersection.c"

void christoffel(double *ch, double *x) {
	/*
		Calculate christoffel symbols for each point in the array x,
		the entries of the array are the x, y, z coordinates of the point
		at which the christoffel symbols are to be computed. Saved as

			x[0] = t_coord
			x[1] = x_coord
			x[2] = y_coord
			x[3] = z_coord

		Values in christoffel symbols are saved as

			ch[0] = f_dot
			ch[1] = gf_dot
			ch[1] = gg_dot

		where d/dx f = x * f_dot and the same holds for g. The christoffel
		symbols at a point can be recovered as lam^k_{ij} = chr[i][j][k]:

			chr[:,:,t] = [[  0, ggx, ggy, ggy],
			 	      [ggx,   0,   0,   0],
			 	      [ggy,   0,   0,   0],
				      [ggz,   0,   0,   0]]
			chr[:,:,x] = [[gfx,   0,   0,   0],
				      [  0,  fx,  fy,  fz],
				      [  0,  fy, -fx,   0],
				      [  0,  fz,   0, -fx]]
			chr[:,:,y] = [[gfy,   0,   0,   0],
				      [  0, -fy,  fx,   0],
				      [  0,  fx,  fy,  fz],
				      [  0,   0,  fz, -fy]]
			chr[:,:,z] = [[gfz,   0,   0,   0],
				      [  0, -fz,   0,  fx],
				      [  0,   0, -fz,  fy],
				      [  0,  fx,  fy,  fz]]
	*/
	double Rinv, rsRinv, one_plus_rinv, one_minus_rinv, plus_squared, plus_cubed, minus_squared, Rcubed, g_dot;
	double finv, ginv;

	Rinv = 1 / sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);

	Rcubed = RS * Rinv * Rinv * Rinv;

	rsRinv = 0.25 * RS * Rinv;
	one_plus_rinv = 1 + rsRinv;
	one_minus_rinv = 1 - rsRinv;

	plus_squared = one_plus_rinv * one_plus_rinv;
	plus_cubed = one_plus_rinv * plus_squared;
	minus_squared = one_minus_rinv * one_minus_rinv;

	finv = 0.5 / (plus_squared * plus_squared);
	ginv = 0.5 * plus_squared / minus_squared;

	g_dot = Rcubed * one_minus_rinv / plus_cubed;

	ch[0] = -finv * Rcubed * plus_cubed;
	ch[1] = -finv * g_dot;
	ch[2] = -ginv * g_dot;
}

void geodesic_deriv(double *x, double *y, double *d, double *ch) {
	/*
		Calculate acceleration of geodesic according to the equation

			d2/ds2 p^k + chr[i][j][k] * d/ds p^i * d/ds p^j = 0
	*/
	double fx, fy, fz, gfx, gfy, gfz, ggx, ggy, ggz;
	double tt, tx, ty, tz, xx, xy, xz, yy, yz, zz;

	// Relevant derivatives
	fx = x[1] * ch[0];
	fy = x[2] * ch[0];
	fz = x[3] * ch[0];

	gfx = x[1] * ch[1];
	gfy = x[2] * ch[1];
	gfz = x[3] * ch[1];

	ggx = x[1] * ch[2];
	ggy = x[2] * ch[2];
	ggz = x[3] * ch[2];

	// Multiplication of entries in derivative vector
	tt = y[0] * y[0];
	tx = y[0] * y[1];
	ty = y[0] * y[2];
	tz = y[0] * y[3];

	xx = y[1] * y[1];
	xy = y[1] * y[2];
	xz = y[1] * y[3];

	yy = y[2] * y[2];
	yz = y[2] * y[3];

	zz = y[3] * y[3];

	d[0] = -2 * (ggx * tx + ggy * ty + ggz * tz);
	d[1] = -(gfx * tt + fx * ( xx - yy - zz) + 2 * (          fy * xy + fz * xz));
	d[2] = -(gfy * tt + fy * (-xx + yy - zz) + 2 * (fx * xy           + fz * yz));
	d[3] = -(gfz * tt + fz * (-xx - yy + zz) + 2 * (fx * xz + fy * yz          ));
}
/*
void rk4_step(double dt, double *restrict x, double *restrict y, double *restrict new_x, double *restrict new_y, double *restrict ch, double *restrict k1, double *restrict k2, double *restrict k3, double *restrict k4) {
	/*
		Computes one runge-kutta step according to the following formula

			k1 = geodesic_deriv(x           , y            )
			k2 = geodesic_deriv(x + dt/2 * y, y + dt/2 * k1)
			k3 = geodesic_deriv(x + dt/2 * y, y + dt/2 * k2)
			k4 = geodesic_deriv(x +   dt * y, y +   dt * k3)

			x_new = x + dt * (1 + 2 * dt / 3) * y
			y_new = y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

	christoffel(ch, x);
	geodesic_deriv(x, y, k1, ch);

	new_x[0] = x[0] + .5 * dt * y[0];
	new_x[1] = x[1] + .5 * dt * y[1];
	new_x[2] = x[2] + .5 * dt * y[2];
	new_x[3] = x[3] + .5 * dt * y[3];

	new_y[0] = y[0] + .5 * dt * k1[0];
	new_y[1] = y[1] + .5 * dt * k1[1];
	new_y[2] = y[2] + .5 * dt * k1[2];
	new_y[3] = y[3] + .5 * dt * k1[3];

	christoffel(ch, new_x);
	geodesic_deriv(new_x, new_y, k2, ch);

	new_y[0] = y[0] + .5 * dt * k2[0];
	new_y[1] = y[1] + .5 * dt * k2[1];
	new_y[2] = y[2] + .5 * dt * k2[2];
	new_y[3] = y[3] + .5 * dt * k2[3];

	geodesic_deriv(new_x, new_y, k3, ch);

	new_x[0] = x[0] + dt * y[0];
	new_x[1] = x[1] + dt * y[1];
	new_x[2] = x[2] + dt * y[2];
	new_x[3] = x[3] + dt * y[3];

	new_y[0] = y[0] + dt * k3[0];
	new_y[1] = y[1] + dt * k3[1];
	new_y[2] = y[2] + dt * k3[2];
	new_y[3] = y[3] + dt * k3[3];

	christoffel(ch, new_x);
	geodesic_deriv(new_x, new_y, k4, ch);

	double step = dt * (1 + dt / 2);
	new_x[0] = x[0] + step * y[0];
	new_x[1] = x[1] + step * y[1];
	new_x[2] = x[2] + step * y[2];
	new_x[3] = x[3] + step * y[3];

	new_y[0] = y[0] + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
	new_y[1] = y[1] + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
	new_y[2] = y[2] + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6;
	new_y[3] = y[3] + dt * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6;
}
*/
double rk34_step(double dt, double *x, double *y, double *new_x, double *new_y, double *ch, double *k1, double *k2, double *k3, double *k3_alt, double *k4) {
	/*
		Computes one runge-kutta step according to the following formula

			k1 = geodesic_deriv(x              , y               )
			k2 = geodesic_deriv(x + .5 * dt * y, y + .5 * dt * k1)
			k3 = geodesic_deriv(x + .5 * dt * y, y + .5 * dt * k2)
			k4 = geodesic_deriv(x +      dt * y, y +      dt * k3)

			new_x = x + dt * (1 + 2 * dt / 3) * y
			new_y = y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

		an embedded RK3 method is computed according to

		 	k3_alt = geodesic_deriv(x + .75 * dt * y, y + 0.75 * dt * k2)

			z_x = x + dt * (1 + dt / 2) * y
			z_y = y + dt * (2 * k1 + 3 * k2 + 4 * k3_alt) / 9

		and the error is estimated according to

			error = norm([new_x - z_x, new_y - z_y])

			new_x - z_x = y * dt * dt / 6
			new_y - z_y = dt * (-k1/6 + k3/3 - 2*k3_alt/3 + k4/6)
	*/
	int i;
	christoffel(ch, x);
	geodesic_deriv(x, y, k1, ch);

	for (i = 0; i < 4; i++) {
		new_x[i] = x[i] + .5 * dt * y[i];
		new_y[i] = y[i] + .5 * dt * k1[i];
	}

	christoffel(ch, new_x);
	geodesic_deriv(new_x, new_y, k2, ch);

	for (i = 0; i < 4; i++)
		new_y[i] = y[i] + .5 * dt * k2[i];

	geodesic_deriv(new_x, new_y, k3, ch);

	for (i = 0; i < 4; i++) {
		new_x[0] = x[i] + .75 * dt * y[i];
		new_y[0] = y[i] + .75 * dt * k2[i];
	}

	christoffel(ch, new_x);
	geodesic_deriv(new_x, new_y, k3_alt, ch);

	for (i = 0; i < 4; i++) {
		new_x[i] = x[i] + dt * y[i];
		new_y[i] = y[i] + dt * k3[i];
	}

	christoffel(ch, new_x);
	geodesic_deriv(new_x, new_y, k4, ch);

	double step = dt * (1 + dt / 2);

	for (i = 0; i < 4; i++) {
		new_x[i] = x[i] + step * y[i];
		new_y[i] = y[i] + dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
	}

	for (i = 0; i < 4; i++) {
		k1[i] = y[i] * dt * dt / 6;
		k2[i] = dt * (-k1[i] / 6 + k3[i] / 3 - 2 * k3_alt[i] / 3 + k4[i] / 6);
	}

	return sqrt(k1[0]*k1[0] + k1[1]*k1[1] + k1[2]*k1[2] + k1[3]*k1[3] + k2[0]*k2[0] + k2[1]*k2[1] + k2[2]*k2[2] + k2[3]*k2[3]);
}

double update_step(params_t ray_params, double dt, double error) {
	return max(min(0.9 * dt * sqrt(ray_params.tol / (2 * error)), ray_params.max_dt), ray_params.min_dt);
}

void trace_adaptive_ray(params_t ray_params, intersection_t* intersection, object_list_t objects, double* x1, double* y1, double* x2, double* y2, double* ch, double* k1, double* k2, double* k3, double* k3_alt, double* k4, double max_time) {
	double dt = 20;
	double new_dt = -1;
	double T = 0;
	double* temp, error;
	bool accepted = false;
	intersection->t = DBL_MAX;
	intersection->collided = false;

	int trials;
	while (T < ray_params.ray_time) {
		new_dt = min(dt, ray_params.ray_time - T);

		// Compute new step
		error = rk34_step(new_dt, x1, y1, x2, y2, ch, k1, k2, k3, k3_alt, k4);

		trials = 0;
		while (error < ray_params.min_tol || error > ray_params.max_tol) {
			// Compute new step and estimate error
			error = rk34_step(new_dt, x1, y1, x2, y2, ch, k1, k2, k3, k3_alt, k4);
			trials++;


			// If new step size was clipped trying to update step size will not help
			if (new_dt == ray_params.min_dt || new_dt == ray_params.max_dt) {
				break;
			}

			// If error is in tolerance accept step
			if (ray_params.min_tol < error && error < ray_params.max_tol) {
				break;
			}

			// If step is not accepted, estimate a new step size
			if ((error < ray_params.min_tol || error > ray_params.max_tol) && (new_dt != ray_params.min_dt || new_dt != ray_params.max_dt)) {
				new_dt = update_step(ray_params, new_dt, error);
			}

			// Safeguard against infinite loops
			if (trials > ray_params.max_trials) {
				break;
			}
		}

		dt = new_dt;
		T += dt;
		temp = x2; x2 = x1; x1 = temp;
		temp = y2; y2 = y1; y1 = temp;

		check_collision(intersection, objects, x2, x1, max_time);

		if (intersection->collided) {
			break;
		}
	}
}
