#include "auxiliarytool.h"

int Kdelta ( int i, int j )
{
	if (i == j) return (1);

	else return (0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void generate_perpendicular_vector ( double * n3, double * res )
{
	double p;


	res [0] = (double) rand () / RAND_MAX; p = (double) rand () / RAND_MAX; if ( p > 0.5 ) res [0] = - res [0];

	res [1] = (double) rand () / RAND_MAX; p = (double) rand () / RAND_MAX; if ( p > 0.5 ) res [1] = - res [1];

	res [2] = (double) rand () / RAND_MAX; p = (double) rand () / RAND_MAX; if ( p > 0.5 ) res [2] = - res [2];


	p = res [0] * n3 [0] + res [1] * n3 [1] + res [2] * n3 [2];

	res [0] = res [0] - p * n3 [0]; res [1] = res [1] - p * n3 [1]; res [2] = res [2] - p * n3 [2];


	p = sqrt (res [0] * res [0] + res [1] * res [1] + res [2] * res [2]);

	res [0] = res [0] / p; res [1] = res [1] / p; res [2] = res [2] / p;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cross_product ( double * v1, double * v2, double * res )
{
	res [0] = v1 [1] * v2 [2] - v1 [2] * v2 [1];       /////////////////// i, j, K cross product: i=j*k, j=k*i, k=i*j

	res [1] = v1 [2] * v2 [0] - v1 [0] * v2 [2];

	res [2] = v1 [0] * v2 [1] - v1 [1] * v2 [0];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void maxwell_1 ( double * array, int array_size )
{
	double r1, r2, tmp1, tmp2;

	int i, j;

	if ( (array_size % 2) == 0 ) j = array_size;

	else
	{
		j = array_size - 1;

gen1:

		r1 = (double) rand () / RAND_MAX; if ( r1 == 0 ) goto gen1;

		r2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( - 2.0 * log (r1) ); tmp2 = 2.0 * pi * r2;

		array [ array_size ] = tmp1 * cos (tmp2);
	}

	for (i = 1; i < j + 1; i+= 2)
	{

gen2:

	r1 = (double) rand () / RAND_MAX; if ( r1 == 0 ) goto gen2;

	r2 = (double) rand () / RAND_MAX;

	tmp1 = sqrt ( - 2.0 * log (r1) ); tmp2 = 2.0 * pi * r2;

	array [i] = tmp1 * cos (tmp2); array [i + 1] = tmp1 * sin (tmp2);

	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void generate_atom_velocity ( int i, int j )
{
	int k, l;

	double a1, a2, tmp1, tmp2;

	/////////////////////////////////////////////////////////////////

	k = atom_key [i][j];

	l = part_key [k];

	/////////////////////////////////////////////////////////////////

genvx:

	a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvx;

	a2 = (double) rand () / RAND_MAX;

	tmp1 = sqrt ( T / MASS [l] ) * sqrt ( -2.0 * log (a1) );

	tmp2 = 2.0 * pi * a2;

	vx [k] = tmp1 * cos ( tmp2 );

genvy:

	a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvy;

	a2 = (double) rand () / RAND_MAX;

	tmp1 = sqrt ( T / MASS [l] ) * sqrt ( -2.0 * log (a1) );

	tmp2 = 2.0 * pi * a2;

	vy [k] = tmp1 * cos ( tmp2 );

genvz:

	a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvz;

	a2 = (double) rand () / RAND_MAX;

	tmp1 = sqrt ( T / MASS [l] ) * sqrt ( -2.0 * log (a1) );

	tmp2 = 2.0 * pi * a2;

	vz [k] = tmp1 * cos ( tmp2 );

	/////////////////////////////////////////////////////////////////

	x_old [k] = x [k] - vx [k] * h;

	y_old [k] = y [k] - vy [k] * h;

	z_old [k] = z [k] - vz [k] * h;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void move_atom ( int i, int j )
{
	int k, l;

	double F1, F2, F3;

	/////////////////////////////////////////////////////////////////

	k = atom_key [i][j];

	l = part_key [k];

	/////////////////////////////////////////////////////////////////

	F1 = ( fx [k] - VISC [l] * vx [k] ) / MASS [l];

	F2 = ( fy [k] - VISC [l] * vy [k] ) / MASS [l];

	F3 = ( fz [k] - VISC [l] * vz [k] ) / MASS [l];

	/////////////////////////////////////////////////////////////////

	x_old [k] = x [k];

	y_old [k] = y [k];

	z_old [k] = z [k];

	/////////////////////////////////////////////////////////////////

	x [k] += h * vx [k] + hsq2 * F1;

	y [k] += h * vy [k] + hsq2 * F2;

	z [k] += h * vz [k] + hsq2 * F3;

	/////////////////////////////////////////////////////////////////

	vx [k] += h * F1;

	vy [k] += h * F2;

	vz [k] += h * F3;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void leap_frog_move_atom ( int i, int j )
{
	int k, l;

	double x_tmp, y_tmp, z_tmp;

	double vx_tmp, vy_tmp, vz_tmp;

	/////////////////////////////////////////////////////////////////

	k = atom_key [i][j];

	l = part_key [k];

	/////////////////////////////////////////////////////////////////

	x_tmp = x [k];

	y_tmp = y [k];

	z_tmp = z [k];

	/////////////////////////////////////////////////////////////////

	vx_tmp = vx [k];

	vy_tmp = vy [k];

	vz_tmp = vz [k];

	/////////////////////////////////////////////////////////////////

	vx [k] = K1 [l] * vx_tmp + K2 [l] * fx [k];

	vy [k] = K1 [l] * vy_tmp + K2 [l] * fy [k];

	vz [k] = K1 [l] * vz_tmp + K2 [l] * fz [k];

	/////////////////////////////////////////////////////////////////

	x [k] = 2.0 * x_tmp - x_old [k] + ( vx [k] - vx_tmp ) * h;

	y [k] = 2.0 * y_tmp - y_old [k] + ( vy [k] - vy_tmp ) * h;

	z [k] = 2.0 * z_tmp - z_old [k] + ( vz [k] - vz_tmp ) * h;

	/////////////////////////////////////////////////////////////////

	x_old [k] = x_tmp;

	y_old [k] = y_tmp;

	z_old [k] = z_tmp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void half_shift ( double * xx )
{
	if ( xx [0] > side_hX ) xx [0] = xx [0] - side_X; if ( xx [0] < - side_hX) xx [0] = xx [0] + side_X;

	if ( xx [1] > side_hY ) xx [1] = xx [1] - side_Y; if ( xx [1] < - side_hY) xx [1] = xx [1] + side_Y;

	if ( xx [2] > side_hZ ) xx [2] = xx [2] - side_Z; if ( xx [2] < - side_hZ) xx [2] = xx [2] + side_Z;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CC_vector ( int i1, int i2, double * xx )
{
	xx [0] = x [i1] - x [i2];

	xx [1] = y [i1] - y [i2];

	xx [2] = z [i1] - z [i2];

	half_shift (xx);
}

double bond_length ( int i1, int i2 )
{
	double v1 [3], r1;

	/////////////////////////////////////////////////////////////////

	v1 [0] = x [i2] - x [i1];

	v1 [1] = y [i2] - y [i1];

	v1 [2] = z [i2] - z [i1];

	/////////////////////////////////////////////////////////////////

	r1 = sqrt ( v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2] );

	/////////////////////////////////////////////////////////////////

	return ( r1 );
}

double valence_angle ( int i1, int i2, int i3 )
{
	double v1 [3], v2 [3], v1_norm, v2_norm, v12_norm, kosinus, theta;

	/////////////////////////////////////////////////////////////////

	v1 [0] = x [i1] - x [i2];

	v1 [1] = y [i1] - y [i2];

	v1 [2] = z [i1] - z [i2];

	/////////////////////////////////////////////////////////////////

	v2 [0] = x [i3] - x [i2];

	v2 [1] = y [i3] - y [i2];

	v2 [2] = z [i3] - z [i2];

	/////////////////////////////////////////////////////////////////

	v1_norm = v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2];

	/////////////////////////////////////////////////////////////////

	v2_norm = v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2];

	/////////////////////////////////////////////////////////////////

	v12_norm = sqrt ( v1_norm * v2_norm );

	/////////////////////////////////////////////////////////////////

	kosinus = ( v1 [0] * v2 [0] + v1 [1] * v2 [1] + v1 [2] * v2 [2] ) / v12_norm;

	theta = acos ( kosinus );

	/////////////////////////////////////////////////////////////////

	return ( theta );
}

double dihedral_angle ( int i1, int i2, int i3, int i4 )
{
	double v1 [3], v2 [3], v3 [3], m [3], n [3], kosinus, sinus, psi, tmp1;

	/////////////////////////////////////////////////////////////////

	v1 [0] = x [i2] - x [i1];

	v1 [1] = y [i2] - y [i1];

	v1 [2] = z [i2] - z [i1];

	/////////////////////////////////////////////////////////////////

	v2 [0] = x [i3] - x [i2];

	v2 [1] = y [i3] - y [i2];

	v2 [2] = z [i3] - z [i2];

	/////////////////////////////////////////////////////////////////

	v3 [0] = x [i4] - x [i3];

	v3 [1] = y [i4] - y [i3];

	v3 [2] = z [i4] - z [i3];

	/////////////////////////////////////////////////////////////////

	cross_product ( v1, v2, m );

	cross_product ( v2, v3, n );

	/////////////////////////////////////////////////////////////////

	tmp1 = sqrt ( v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2] );

	kosinus = m [0] * n [0] + m [1] * n [1] + m [2] * n [2];

	sinus = ( v1 [0] * n [0] + v1 [1] * n [1] + v1 [2] * n [2] ) * tmp1;

	psi = atan2 ( sinus, kosinus );

	/////////////////////////////////////////////////////////////////

	return ( psi );
}
