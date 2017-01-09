#include "boxoperation.h"

int box_init()
{
  BORDER_MIN = 1.5 * ( 1.2 * D [4] );
	BORDER_MAX = BORDER_MIN + 2.0;
  return 1;
}

int set_box ()
{
	int i, J1, J2;

	double x_min, y_min, z_min;

	double x_max, y_max, z_max;

	double A, B, C, F, G, H;

	double a0, a1, a2;

	double * ep1 = NULL, * ep2 = NULL, * ep3 = NULL;

	double n1 [3], n2 [3], n3 [3];

	///////////Protein dimensions///////////
	////////////////////////////////////////

	a0 = 0; a1 = 0; a2 = 0;

	for (i = 1; i < NTP_atom + 1; i++)
	{
		a0 += x [i];

		a1 += y [i];

		a2 += z [i];
	}

	a0 /= NTP_atom;

	a1 /= NTP_atom;

	a2 /= NTP_atom;

	///////////////////////////////////////////////

	A = 0; B = 0; C = 0; F = 0; G = 0; H = 0;

	for (i = 1; i < NTP_atom + 1; i++)
	{
		x_min = x [i] - a0;

		y_min = y [i] - a1;

		z_min = z [i] - a2;

		////////////////////////////////////////

		A += y_min * y_min + z_min * z_min;

		B += x_min * x_min + z_min * z_min;

		C += x_min * x_min + y_min * y_min;

		F += y_min * z_min;

		G += x_min * z_min;

		H += x_min * y_min;
	}

	///////////////////////////////////////////////

	if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6  )
	{
		a0 = sqrt (B * B - 2.0 * A * B + A * A + 4.0 * H * H);

		/////////////////////////////////

		x_max = 0.5 * (B + A + a0);

		n1 [0] = 0.5 * (B - A - a0) / H;

		n1 [1] = 1;

		n1 [2] = 0;

		a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

		n1 [0] /= a1;

		n1 [1] /= a1;

		n1 [2] /= a1;

		/////////////////////////////////

		y_max = 0.5 * (B + A - a0);

		n2 [0] = 0.5 * (B - A + a0) / H;

		n2 [1] = 1;

		n2 [2] = 0;

		a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

		n2 [0] /= a1;

		n2 [1] /= a1;

		n2 [2] /= a1;

		/////////////////////////////////

		z_max = C;

		n3 [0] = 0;

		n3 [1] = 0;

		n3 [2] = 1;

		/////////////////////////////////

		if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }

		else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

		if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

		else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

		else { z_min = z_max; ep3 = n3; }

		/////////////////////////////////

		if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; }

		if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
	}

	///////////////////////////////////////////////

	else if ( fabs (F) < 1.0e-6 && fabs (H) < 1.0e-6  )
	{
		a0 = sqrt (C * C - 2.0 * A * C + A * A + 4.0 * G * G);

		/////////////////////////////////

		x_max = 0.5 * (C + A + a0);

		n1 [0] = 0.5 * (C - A - a0) / G;

		n1 [1] = 0;

		n1 [2] = 1;

		a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

		n1 [0] /= a1;

		n1 [1] /= a1;

		n1 [2] /= a1;

		/////////////////////////////////

		y_max = 0.5 * (C + A - a0);

		n2 [0] = 0.5 * (C - A + a0) / G;

		n2 [1] = 0;

		n2 [2] = 1;

		a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

		n2 [0] /= a1;

		n2 [1] /= a1;

		n2 [2] /= a1;

		/////////////////////////////////

		z_max = B;

		n3 [0] = 0;

		n3 [1] = 1;

		n3 [2] = 0;

		/////////////////////////////////

		if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }

		else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

		if ( z_max < x_min )
    { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

		else if ( z_max < y_min )
    { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

		else
    { z_min = z_max; ep3 = n3; }

		/////////////////////////////////

		if ( ep1 [0] < 0 )
    { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; }

		if ( ep2 [1] < 0 )
    { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
	}

	///////////////////////////////////////////////

	else if ( fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
	{
		a0 = sqrt (C * C - 2.0 * B * C + B * B + 4.0 * F * F);

		/////////////////////////////////

		x_max = 0.5 * (C + B + a0);

		n1 [0] = 0;

		n1 [1] = 0.5 * (C - B - a0) / F;

		n1 [2] = 1;

		a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

		n1 [0] /= a1;

		n1 [1] /= a1;

		n1 [2] /= a1;

		/////////////////////////////////

		y_max = 0.5 * (C + B - a0);

		n2 [0] = 0;

		n2 [1] = 0.5 * (C - B + a0) / F;

		n2 [2] = 1;

		a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

		n2 [0] /= a1;

		n2 [1] /= a1;

		n2 [2] /= a1;

		/////////////////////////////////

		z_max = A;

		n3 [0] = 1;

		n3 [1] = 0;

		n3 [2] = 0;

		/////////////////////////////////

		if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }

		else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

		if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

		else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

		else { z_min = z_max; ep3 = n3; }

		/////////////////////////////////

		if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; }

		if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
	}

	///////////////////////////////////////////////

	else if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
	{
		x_max = A;

		n1 [0] = 1;

		n1 [1] = 0;

		n1 [2] = 0;

		/////////////////////////////////

		y_max = B;

		n2 [0] = 0;

		n2 [1] = 1;

		n2 [2] = 0;

		/////////////////////////////////

		z_max = C;

		n3 [0] = 0;

		n3 [1] = 0;

		n3 [2] = 1;

		/////////////////////////////////

		if ( x_max < y_max ) { x_min = x_max; y_min = y_max; ep1 = n1; ep2 = n2; }

		else { x_min = y_max; y_min = x_max; ep1 = n2; ep2 = n1; }

		if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

		else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

		else { z_min = z_max; ep3 = n3; }

		/////////////////////////////////

		if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; }

		if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
	}

	///////////////////////////////////////////////

	else
	{
		a2 = - A - B - C;

		a1 = A * B + A * C + B * C - H * H - G * G - F * F;

		a0 = A * F * F + B * G * G + C * H * H + 2.0 * H * G * F - A * B * C;

		/////////////////////////////////

		x_min = - a2 / 3.0;

		x_max = ( 9.0 * a1 * a2 - 27 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;

		y_max = x_min * x_min - a1 / 3.0;

		y_min = 2.0 * sqrt ( y_max );

		z_min = acos ( x_max * pow ( y_max, - 1.5 ) ) / 3.0;

		a0 = pi / 3.0;

		/////////////////////////////////

		x_max = x_min + y_min * cos ( z_min );

		y_max = x_min - y_min * cos ( z_min - a0 );

		z_max = x_min - y_min * cos ( z_min + a0 );

		/////////////////////////////////

		if ( x_max < y_max ) { x_min = x_max; y_min = y_max; }

		else { x_min = y_max; y_min = x_max; }

		if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; }

		else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; }

		else z_min = z_max;

		/////////////////////////////////

		x_max = (A - x_min) * F + G * H;

		y_max = (B - x_min) * G + F * H;

		z_max = (C - x_min) * H + G * F;

		a0 = x_max / y_max;

		a1 = x_max / z_max;

		a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

		n1 [0] = 1.0 / a2;

		n1 [1] = a0 / a2;

		n1 [2] = a1 / a2;

		ep1 = n1;

		/////////////////////////////////

		x_max = (A - y_min) * F + G * H;

		y_max = (B - y_min) * G + F * H;

		z_max = (C - y_min) * H + G * F;

		a0 = y_max / x_max;

		a1 = y_max / z_max;

		a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

		n2 [0] = a0 / a2;

		n2 [1] = 1.0 / a2;

		n2 [2] = a1 / a2;

		ep2 = n2;

		/////////////////////////////////

		ep3 = n3;
	}

	///////////////////////////////////////////////

	cross_product ( ep1, ep2, ep3 );

	/////////////////////////////////
	/////////////////Set box////////////////
	////////////////////////////////////////

	for (i = 1; i < NTP_atom + 1; i++)
	{
		x_max = x [i] * ep1 [0] + y [i] * ep1 [1] + z [i] * ep1 [2];

		y_max = x [i] * ep2 [0] + y [i] * ep2 [1] + z [i] * ep2 [2];

		z_max = x [i] * ep3 [0] + y [i] * ep3 [1] + z [i] * ep3 [2];

		x [i] = x_max;

		y [i] = y_max;

		z [i] = z_max;
	}

	/////////////////////////////////////////////////

	x_min = x [1]; y_min = y [1]; z_min = z [1];

	x_max = x [1]; y_max = y [1]; z_max = z [1];

	for (i = 2; i < NTP_atom + 1; i++)
	{
		if ( x [i] < x_min ) x_min = x [i];

		if ( x [i] > x_max ) x_max = x [i];

		/////////////////////////////////////////////////

		if ( y [i] < y_min ) y_min = y [i];

		if ( y [i] > y_max ) y_max = y [i];

		/////////////////////////////////////////////////

		if ( z [i] < z_min ) z_min = z [i];

		if ( z [i] > z_max ) z_max = z [i];
	}

	/////////////////////////////////////////////////

	side_X = x_max - x_min + 2.0 * BORDER_MAX;

	side_Y = y_max - y_min + 2.0 * BORDER_MAX;

	side_Z = z_max - z_min + 2.0 * BORDER_MAX;

	side_XYZ = side_X * side_Y * side_Z;

	/////////////////////////////////////////////////

	side_hX = 0.5 * side_X;

	side_hY = 0.5 * side_Y;

	side_hZ = 0.5 * side_Z;

	/////////////////////////////////////////////////

	for (i = 1; i < NTP_atom + 1; i++)
	{
		x [i] += BORDER_MAX - x_min;

		y [i] += BORDER_MAX - y_min;

		z [i] += BORDER_MAX - z_min;

		/////////////////////////////////////////////////

		x_max = x_old [i] * ep1 [0] + y_old [i] * ep1 [1] + z_old [i] * ep1 [2] + BORDER_MAX - x_min;

		y_max = x_old [i] * ep2 [0] + y_old [i] * ep2 [1] + z_old [i] * ep2 [2] + BORDER_MAX - y_min;

		z_max = x_old [i] * ep3 [0] + y_old [i] * ep3 [1] + z_old [i] * ep3 [2] + BORDER_MAX - z_min;

		x_old [i] = x_max;

		y_old [i] = y_max;

		z_old [i] = z_max;

		/////////////////////////////////////////////////

		x_max = vx [i] * ep1 [0] + vy [i] * ep1 [1] + vz [i] * ep1 [2];

		y_max = vx [i] * ep2 [0] + vy [i] * ep2 [1] + vz [i] * ep2 [2];

		z_max = vx [i] * ep3 [0] + vy [i] * ep3 [1] + vz [i] * ep3 [2];

		vx [i] = x_max;

		vy [i] = y_max;

		vz [i] = z_max;
	}

	/////////////////////////////////////////////////

	for (i = 1; i < NBP_atom + 1; i++)
	{
		if ( RIGID_SET [i] > RIGID_END [i] ) continue;

		/////////////////////////////////////////////////

		J1 = 3 * (i - 1);

		J2 = 9 * (i - 1);

		/////////////////////////////////////////////////

		x_max = CMS [J1 + 1] * ep1 [0] + CMS [J1 + 2] * ep1 [1] + CMS [J1 + 3] * ep1 [2] + BORDER_MAX - x_min;

		y_max = CMS [J1 + 1] * ep2 [0] + CMS [J1 + 2] * ep2 [1] + CMS [J1 + 3] * ep2 [2] + BORDER_MAX - y_min;

		z_max = CMS [J1 + 1] * ep3 [0] + CMS [J1 + 2] * ep3 [1] + CMS [J1 + 3] * ep3 [2] + BORDER_MAX - z_min;

		CMS [J1 + 1] = x_max;

		CMS [J1 + 2] = y_max;

		CMS [J1 + 3] = z_max;

		/////////////////////////////////////////////////

		x_max = VCMS [J1 + 1] * ep1 [0] + VCMS [J1 + 2] * ep1 [1] + VCMS [J1 + 3] * ep1 [2];

		y_max = VCMS [J1 + 1] * ep2 [0] + VCMS [J1 + 2] * ep2 [1] + VCMS [J1 + 3] * ep2 [2];

		z_max = VCMS [J1 + 1] * ep3 [0] + VCMS [J1 + 2] * ep3 [1] + VCMS [J1 + 3] * ep3 [2];

		VCMS [J1 + 1] = x_max;

		VCMS [J1 + 2] = y_max;

		VCMS [J1 + 3] = z_max;

		/////////////////////////////////////////////////

		x_max = AXES [J2 + 1] * ep1 [0] + AXES [J2 + 2] * ep1 [1] + AXES [J2 + 3] * ep1 [2];

		y_max = AXES [J2 + 1] * ep2 [0] + AXES [J2 + 2] * ep2 [1] + AXES [J2 + 3] * ep2 [2];

		z_max = AXES [J2 + 1] * ep3 [0] + AXES [J2 + 2] * ep3 [1] + AXES [J2 + 3] * ep3 [2];

		AXES [J2 + 1] = x_max;

		AXES [J2 + 2] = y_max;

		AXES [J2 + 3] = z_max;

		/////////////////////////////////////////////////

		x_max = AXES [J2 + 4] * ep1 [0] + AXES [J2 + 5] * ep1 [1] + AXES [J2 + 6] * ep1 [2];

		y_max = AXES [J2 + 4] * ep2 [0] + AXES [J2 + 5] * ep2 [1] + AXES [J2 + 6] * ep2 [2];

		z_max = AXES [J2 + 4] * ep3 [0] + AXES [J2 + 5] * ep3 [1] + AXES [J2 + 6] * ep3 [2];

		AXES [J2 + 4] = x_max;

		AXES [J2 + 5] = y_max;

		AXES [J2 + 6] = z_max;

		/////////////////////////////////////////////////

		x_max = AXES [J2 + 7] * ep1 [0] + AXES [J2 + 8] * ep1 [1] + AXES [J2 + 9] * ep1 [2];

		y_max = AXES [J2 + 7] * ep2 [0] + AXES [J2 + 8] * ep2 [1] + AXES [J2 + 9] * ep2 [2];

		z_max = AXES [J2 + 7] * ep3 [0] + AXES [J2 + 8] * ep3 [1] + AXES [J2 + 9] * ep3 [2];

		AXES [J2 + 7] = x_max;

		AXES [J2 + 8] = y_max;

		AXES [J2 + 9] = z_max;

		/////////////////////////////////////////////////

		x_max = W [J1 + 1] * ep1 [0] + W [J1 + 2] * ep1 [1] + W [J1 + 3] * ep1 [2];

		y_max = W [J1 + 1] * ep2 [0] + W [J1 + 2] * ep2 [1] + W [J1 + 3] * ep2 [2];

		z_max = W [J1 + 1] * ep3 [0] + W [J1 + 2] * ep3 [1] + W [J1 + 3] * ep3 [2];

		W [J1 + 1] = x_max;

		W [J1 + 2] = y_max;

		W [J1 + 3] = z_max;
	}
  return 1;
}
