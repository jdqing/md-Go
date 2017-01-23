#include "boxoperation.h"

void add_crowder ( int tp )
{
	int i, j, k;

	int i2, j2, k2;

	int J1, J2;

	int In1, In2, In12;

	double tmp1, tmp2, tmp3, xx [3];

	double A, B, C, F, G, H;

	double a0, a1, a2;

	double * ep1 = NULL, * ep2 = NULL, * ep3 = NULL;

	double n1 [3], n2 [3], n3 [3];

	/////////////////////////////////////////////////////////////////

	NB_atom += 1;

	i = NB_atom;

	/////////////////////////////////////////////////////////////////

	crowder_mass [tp] += 1;

	crowder_content [tp] = (int *) realloc ( crowder_content [tp], (crowder_mass [tp] + 1) * sizeof(int) );

	crowder_content [tp][ crowder_mass [tp] ] = i;

	/////////////////////////////////////////////////////////////////

	J1 = 3 * (i - 1);

	J2 = 9 * (i - 1);

	/////////////////////////////////////////////////////////////////

	crowder_key = (int *) realloc ( crowder_key, (i + 1) * sizeof(int) );

	crowder_key [i] = tp;

	/////////////////////////////////////////////////////////////////

	NS_atom = (int *) realloc ( NS_atom, (i + 1) * sizeof(int) );

	NS_atom [i] = m_crwd [tp] - 1;

	/////////////////////////////////////////////////////////////////

	atom_key = (int **) realloc ( atom_key, (i + 1) * sizeof(int *) );

	atom_key [i] = NULL;

	atom_key [i] = (int *) calloc ( NS_atom [i] + 1, sizeof(int) );

	/////////////////////////////////////////////////////////////////

	RIGID_SET = (int *) realloc ( RIGID_SET, (i + 1) * sizeof(int) );

	RIGID_END = (int *) realloc ( RIGID_END, (i + 1) * sizeof(int) );

	RIGID_SET [i] = RIGID_SET_crwd [tp];

	RIGID_END [i] = RIGID_END_crwd [tp];

	/////////////////////////////////////////////////////////////////

	for (j = 0; j < NS_atom [i] + 1; j++)
	{
		NT_atom += 1;

		atom_key [i][j] = NT_atom;

		/////////////////////////////////////////////////////////////////

		part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

		part_key [ NT_atom ] = part_key_crwd [tp][j];

		/////////////////////////////////////////////////////////////////

		maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );

		maxi_key [ NT_atom ] = maxi_key_crwd [tp][j];

		/////////////////////////////////////////////////////////////////

		INDX = (int *) realloc ( INDX, (NT_atom + 1) * sizeof(int) );

		JNDX = (int *) realloc ( JNDX, (NT_atom + 1) * sizeof(int) );

		INDX [ NT_atom ] = i;

		JNDX [ NT_atom ] = j;
	}

	/////////////////////////////////////////////////////////////////

	IR1 = (double *) realloc ( IR1, (i + 1) * sizeof(double) );

	IR1 [i] = 0;

	IR2 = (double *) realloc ( IR2, (i + 1) * sizeof(double) );

	IR2 [i] = 0;

	IR3 = (double *) realloc ( IR3, (i + 1) * sizeof(double) );

	IR3 [i] = 0;

	/////////////////////////////////////////////////////////////////

	CMS = (double *) realloc ( CMS, (3 * i + 1) * sizeof(double) );

	CMS [J1 + 1] = 0;

	CMS [J1 + 2] = 0;

	CMS [J1 + 3] = 0;

	/////////////////////////////////////////////////////////////////

	VCMS = (double *) realloc ( VCMS, (3 * i + 1) * sizeof(double) );

	VCMS [J1 + 1] = 0;

	VCMS [J1 + 2] = 0;

	VCMS [J1 + 3] = 0;

	/////////////////////////////////////////////////////////////////

	AXES = (double *) realloc ( AXES, (9 * i + 1) * sizeof(double) );

	AXES [J2 + 1] = 0;

	AXES [J2 + 2] = 0;

	AXES [J2 + 3] = 0;

	AXES [J2 + 4] = 0;

	AXES [J2 + 5] = 0;

	AXES [J2 + 6] = 0;

	AXES [J2 + 7] = 0;

	AXES [J2 + 8] = 0;

	AXES [J2 + 9] = 0;

	/////////////////////////////////////////////////////////////////

	W = (double *) realloc ( W, (3 * i + 1) * sizeof(double) );

	W [J1 + 1] = 0;

	W [J1 + 2] = 0;

	W [J1 + 3] = 0;

	/////////////////////////////////////////////////////////////////

	x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );

	y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

	z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );

	x [ NT_atom ] = 0;

	y [ NT_atom ] = 0;

	z [ NT_atom ] = 0;

	/////////////////////////////////////////////////////////////////

	x_old = (double *) realloc ( x_old, (NT_atom + 1) * sizeof(double) );

	y_old = (double *) realloc ( y_old, (NT_atom + 1) * sizeof(double) );

	z_old = (double *) realloc ( z_old, (NT_atom + 1) * sizeof(double) );

	x_old [ NT_atom ] = 0;

	y_old [ NT_atom ] = 0;

	z_old [ NT_atom ] = 0;

	/////////////////////////////////////////////////////////////////

	vx = (double *) realloc ( vx, (NT_atom + 1) * sizeof(double) );

	vy = (double *) realloc ( vy, (NT_atom + 1) * sizeof(double) );

	vz = (double *) realloc ( vz, (NT_atom + 1) * sizeof(double) );

	vx [ NT_atom ] = 0;

	vy [ NT_atom ] = 0;

	vz [ NT_atom ] = 0;

	/////////////////////////////////////////////////////////////////

	fx = (double *) realloc ( fx, (NT_atom + 1) * sizeof(double) );

	fy = (double *) realloc ( fy, (NT_atom + 1) * sizeof(double) );

	fz = (double *) realloc ( fz, (NT_atom + 1) * sizeof(double) );

	fx [ NT_atom ] = 0;

	fy [ NT_atom ] = 0;

	fz [ NT_atom ] = 0;

	/////////////////////////////////////////////////////////////////

	RMASS = (double *) realloc ( RMASS, (i + 1) * sizeof(double) );

	RMASS [i] = 0;

	/////////////////////////////////////////////////////////////////

	RVISC = (double *) realloc ( RVISC, (i + 1) * sizeof(double) );

	RVISC [i] = 0;

	/////////////////////////////////////////////////////////////////

	RXYZ = (double *) realloc ( RXYZ, (3 * i + 1) * sizeof(double) );

	RXYZ [J1 + 1] = 0;

	RXYZ [J1 + 2] = 0;

	RXYZ [J1 + 3] = 0;

	/////////////////////////////////////////////////////////////////

	RX = (double **) realloc ( RX, (i + 1) * sizeof(double *) );

	RX [i] = NULL;

	RX [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	RY = (double **) realloc ( RY, (i + 1) * sizeof(double *) );

	RY [i] = NULL;

	RY [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	RZ = (double **) realloc ( RZ, (i + 1) * sizeof(double *) );

	RZ [i] = NULL;

	RZ [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	AV = (double *) realloc ( AV, (i + 1) * sizeof(double) );

	AV [i] = 0;

	BV = (double *) realloc ( BV, (i + 1) * sizeof(double) );

	BV [i] = 0;

	CV = (double *) realloc ( CV, (i + 1) * sizeof(double) );

	CV [i] = 0;

	FV = (double *) realloc ( FV, (i + 1) * sizeof(double) );

	FV [i] = 0;

	GV = (double *) realloc ( GV, (i + 1) * sizeof(double) );

	GV [i] = 0;

	HV = (double *) realloc ( HV, (i + 1) * sizeof(double) );

	HV [i] = 0;

	/////////////////////////////////////////////////////////////////

gencoords:

	a1 = (double) rand () / RAND_MAX;

	tmp1 = a1 * side_X;

	a1 = (double) rand () / RAND_MAX;

	tmp2 = a1 * side_Y;

	a1 = (double) rand () / RAND_MAX;

	tmp3 = a1 * side_Z;

	/////////////////////////////////////////////////////////////////

	n1 [0] = rand ();

	a1 = (double) rand () / RAND_MAX;

	if (a1 > 0.5) n1 [0] = - n1 [0];

	/////////////////////////////////////////////////

	n1 [1] = rand ();

	a1 = (double) rand () / RAND_MAX;

	if (a1 > 0.5) n1 [1] = - n1 [1];

	/////////////////////////////////////////////////

	n1 [2] = rand ();

	a1 = (double) rand () / RAND_MAX;

	if (a1 > 0.5) n1 [2] = - n1 [2];

	/////////////////////////////////////////////////

	a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

	n1 [0] /= a1;

	n1 [1] /= a1;

	n1 [2] /= a1;

	/////////////////////////////////////////////////

	generate_perpendicular_vector ( n1, n2 );

	/////////////////////////////////////////////////

	cross_product ( n1, n2, n3 );

	/////////////////////////////////////////////////

	for (j = 0; j < NS_atom [i] + 1; j++)
	{
		k = atom_key [i][j];

		In1 = part_key [k];

		/////////////////////////////////////////////////////////////////

		x [k] = tmp1 + n1 [0] * RX_crwd [tp][j] + n2 [0] * RY_crwd [tp][j] + n3 [0] * RZ_crwd [tp][j];

		y [k] = tmp2 + n1 [1] * RX_crwd [tp][j] + n2 [1] * RY_crwd [tp][j] + n3 [1] * RZ_crwd [tp][j];

		z [k] = tmp3 + n1 [2] * RX_crwd [tp][j] + n2 [2] * RY_crwd [tp][j] + n3 [2] * RZ_crwd [tp][j];

		/////////////////////////////////////////////////////////////////

		for (i2 = 1; i2 < i; i2++)
		{
			for (j2 = 0; j2 < NS_atom [i2] + 1; j2++)
			{
				k2 = atom_key [i2][j2];

				In2 = part_key [k2];

				if (In1 < In2) In12 = ns * In1 - In1 * (In1 + 1) / 2 + In2;

				else In12 = ns * In2 - In2 * (In2 + 1) / 2 + In1;

				/////////////////////////////////////////////////////////////////

				xx [0] = x [k] - x [k2];

				xx [1] = y [k] - y [k2];

				xx [2] = z [k] - z [k2];

				half_shift ( xx );

				a2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

				if ( a2 < 0.7 * PD_sq [ In12 ] ) goto gencoords;
			}
		}
	}

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////

	for (j = 0; j < RIGID_SET [i]; j++) generate_atom_velocity ( i, j );

	for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++) generate_atom_velocity ( i, j );

	/////////////////////////////////////////////////////////////////

	if ( RIGID_SET [i] <= RIGID_END [i] )
	{
		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			In1 = part_key [k];

			/////////////////////////////////////////////////////////////////

			RMASS [i] += MASS [ In1 ];

			RVISC [i] += VISC [ In1 ];

			/////////////////////////////////////////////////////////////////

			CMS [J1 + 1] += MASS [ In1 ] * x [k];

			CMS [J1 + 2] += MASS [ In1 ] * y [k];

			CMS [J1 + 3] += MASS [ In1 ] * z [k];
		}

		CMS [J1 + 1] /= RMASS [i];

		CMS [J1 + 2] /= RMASS [i];

		CMS [J1 + 3] /= RMASS [i];

		/////////////////////////////////////////////////////////////////

		A = 0; B = 0; C = 0; F = 0; G = 0; H = 0;

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			In1 = part_key [k];

			/////////////////////////////////////////////////////////////////

			xx [0] = x [k] - CMS [J1 + 1];

			xx [1] = y [k] - CMS [J1 + 2];

			xx [2] = z [k] - CMS [J1 + 3];

			/////////////////////////////////////////////////////////////////

			A += MASS [ In1 ] * ( xx [1] * xx [1] + xx [2] * xx [2] );

			B += MASS [ In1 ] * ( xx [0] * xx [0] + xx [2] * xx [2] );

			C += MASS [ In1 ] * ( xx [0] * xx [0] + xx [1] * xx [1] );

			F += MASS [ In1 ] * xx [1] * xx [2];

			G += MASS [ In1 ] * xx [0] * xx [2];

			H += MASS [ In1 ] * xx [0] * xx [1];
		}

		/////////////////////////////////////////////////////////////////

		if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6  )
		{
			a0 = sqrt (B * B - 2.0 * A * B + A * A + 4.0 * H * H);

			/////////////////////////////////

			xx [0] = 0.5 * (B + A + a0);

			n1 [0] = 0.5 * (B - A - a0) / H;

			n1 [1] = 1;

			n1 [2] = 0;

			a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

			n1 [0] /= a1;

			n1 [1] /= a1;

			n1 [2] /= a1;

			/////////////////////////////////

			xx [1] = 0.5 * (B + A - a0);

			n2 [0] = 0.5 * (B - A + a0) / H;

			n2 [1] = 1;

			n2 [2] = 0;

			a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

			n2 [0] /= a1;

			n2 [1] /= a1;

			n2 [2] /= a1;

			/////////////////////////////////

			xx [2] = C;

			n3 [0] = 0;

			n3 [1] = 0;

			n3 [2] = 1;

			/////////////////////////////////

			if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }

			else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

			if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

			else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

			else { tmp3 = xx [2]; ep3 = n3; }
		}

		///////////////////////////////////////////////

		else if ( fabs (F) < 1.0e-6 && fabs (H) < 1.0e-6  )
		{
			a0 = sqrt (C * C - 2.0 * A * C + A * A + 4.0 * G * G);

			/////////////////////////////////

			xx [0] = 0.5 * (C + A + a0);

			n1 [0] = 0.5 * (C - A - a0) / G;

			n1 [1] = 0;

			n1 [2] = 1;

			a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

			n1 [0] /= a1;

			n1 [1] /= a1;

			n1 [2] /= a1;

			/////////////////////////////////

			xx [1] = 0.5 * (C + A - a0);

			n2 [0] = 0.5 * (C - A + a0) / G;

			n2 [1] = 0;

			n2 [2] = 1;

			a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

			n2 [0] /= a1;

			n2 [1] /= a1;

			n2 [2] /= a1;

			/////////////////////////////////

			xx [2] = B;

			n3 [0] = 0;

			n3 [1] = 1;

			n3 [2] = 0;

			/////////////////////////////////

			if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }

			else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

			if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

			else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

			else { tmp3 = xx [2]; ep3 = n3; }
		}

		///////////////////////////////////////////////

		else if ( fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
		{
			a0 = sqrt (C * C - 2.0 * B * C + B * B + 4.0 * F * F);

			/////////////////////////////////

			xx [0] = 0.5 * (C + B + a0);

			n1 [0] = 0;

			n1 [1] = 0.5 * (C - B - a0) / F;

			n1 [2] = 1;

			a1 = sqrt ( n1 [0] * n1 [0] + n1 [1] * n1 [1] + n1 [2] * n1 [2] );

			n1 [0] /= a1;

			n1 [1] /= a1;

			n1 [2] /= a1;

			/////////////////////////////////

			xx [1] = 0.5 * (C + B - a0);

			n2 [0] = 0;

			n2 [1] = 0.5 * (C - B + a0) / F;

			n2 [2] = 1;

			a1 = sqrt ( n2 [0] * n2 [0] + n2 [1] * n2 [1] + n2 [2] * n2 [2] );

			n2 [0] /= a1;

			n2 [1] /= a1;

			n2 [2] /= a1;

			/////////////////////////////////

			xx [2] = A;

			n3 [0] = 1;

			n3 [1] = 0;

			n3 [2] = 0;

			/////////////////////////////////

			if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }

			else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

			if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

			else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

			else { tmp3 = xx [2]; ep3 = n3; }
		}

		///////////////////////////////////////////////

		else if ( fabs (F) < 1.0e-6 && fabs (G) < 1.0e-6 && fabs (H) < 1.0e-6  )
		{
			xx [0] = A;

			n1 [0] = 1;

			n1 [1] = 0;

			n1 [2] = 0;

			/////////////////////////////////

			xx [1] = B;

			n2 [0] = 0;

			n2 [1] = 1;

			n2 [2] = 0;

			/////////////////////////////////

			xx [2] = C;

			n3 [0] = 0;

			n3 [1] = 0;

			n3 [2] = 1;

			/////////////////////////////////

			if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; ep1 = n1; ep2 = n2; }

			else { tmp1 = xx [1]; tmp2 = xx [0]; ep1 = n2; ep2 = n1; }

			if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; ep3 = ep2; ep2 = ep1; ep1 = n3; }

			else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; ep3 = ep2; ep2 = n3; }

			else { tmp3 = xx [2]; ep3 = n3; }
		}

		///////////////////////////////////////////////

		else
		{
			a2 = - A - B - C;

			a1 = A * B + A * C + B * C - H * H - G * G - F * F;

			a0 = A * F * F + B * G * G + C * H * H + 2.0 * H * G * F - A * B * C;

			/////////////////////////////////

			tmp1 = - a2 / 3.0;

			xx [0] = ( 9.0 * a1 * a2 - 27 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;

			xx [1] = tmp1 * tmp1 - a1 / 3.0;

			tmp2 = 2.0 * sqrt ( xx [1] );

			tmp3 = acos ( xx [0] * pow ( xx [1], - 1.5 ) ) / 3.0;

			a0 = pi / 3.0;

			/////////////////////////////////

			xx [0] = tmp1 + tmp2 * cos ( tmp3 );

			xx [1] = tmp1 - tmp2 * cos ( tmp3 - a0 );

			xx [2] = tmp1 - tmp2 * cos ( tmp3 + a0 );

			/////////////////////////////////

			if ( xx [0] < xx [1] ) { tmp1 = xx [0]; tmp2 = xx [1]; }

			else { tmp1 = xx [1]; tmp2 = xx [0]; }

			if ( xx [2] < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = xx [2]; }

			else if ( xx [2] < tmp2 ) { tmp3 = tmp2; tmp2 = xx [2]; }

			else tmp3 = xx [2];

			/////////////////////////////////

			xx [0] = (A - tmp1) * F + G * H;

			xx [1] = (B - tmp1) * G + F * H;

			xx [2] = (C - tmp1) * H + G * F;

			a0 = xx [0] / xx [1];

			a1 = xx [0] / xx [2];

			a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

			n1 [0] = 1.0 / a2;

			n1 [1] = a0 / a2;

			n1 [2] = a1 / a2;

			ep1 = n1;

			/////////////////////////////////

			xx [0] = (A - tmp2) * F + G * H;

			xx [1] = (B - tmp2) * G + F * H;

			xx [2] = (C - tmp2) * H + G * F;

			a0 = xx [1] / xx [0];

			a1 = xx [1] / xx [2];

			a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

			n2 [0] = a0 / a2;

			n2 [1] = 1.0 / a2;

			n2 [2] = a1 / a2;

			ep2 = n2;

			/////////////////////////////////

			ep3 = n3;
		}

		/////////////////////////////////////////////////////////////////

		IR1 [i] = tmp1;

		IR2 [i] = tmp2;

		IR3 [i] = tmp3;

		/////////////////////////////////////////////////////////////////

		if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; }

		if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }

		cross_product ( ep1, ep2, ep3 );

		/////////////////////////////////////////////////////////////////

		AXES [J2 + 1] = ep1 [0]; AXES [J2 + 2] = ep1 [1]; AXES [J2 + 3] = ep1 [2];

		AXES [J2 + 4] = ep2 [0]; AXES [J2 + 5] = ep2 [1]; AXES [J2 + 6] = ep2 [2];

		AXES [J2 + 7] = ep3 [0]; AXES [J2 + 8] = ep3 [1]; AXES [J2 + 9] = ep3 [2];

		/////////////////////////////////////////////////////////////////

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			In1 = part_key [k];

			/////////////////////////////////////////////////////////////////

			xx [0] = x [k] - CMS [J1 + 1];

			xx [1] = y [k] - CMS [J1 + 2];

			xx [2] = z [k] - CMS [J1 + 3];

			/////////////////////////////////////////////////////////////////

			RX [i][j] = xx [0] * ep1 [0] + xx [1] * ep1 [1] + xx [2] * ep1 [2];

			RY [i][j] = xx [0] * ep2 [0] + xx [1] * ep2 [1] + xx [2] * ep2 [2];

			RZ [i][j] = xx [0] * ep3 [0] + xx [1] * ep3 [1] + xx [2] * ep3 [2];

			/////////////////////////////////////////////////////////////////

			RXYZ [J1 + 1] += VISC [ In1 ] * RX [i][j];

			RXYZ [J1 + 2] += VISC [ In1 ] * RY [i][j];

			RXYZ [J1 + 3] += VISC [ In1 ] * RZ [i][j];

			/////////////////////////////////////////////////////////////////

			AV [i] += VISC [ In1 ]  * ( RY [i][j] * RY [i][j] + RZ [i][j] * RZ [i][j] );

			BV [i] += VISC [ In1 ] * ( RX [i][j] * RX [i][j] + RZ [i][j] * RZ [i][j] );

			CV [i] += VISC [ In1 ] * ( RX [i][j] * RX [i][j] + RY [i][j] * RY [i][j] );

			FV [i] += VISC [ In1 ] * RY [i][j] * RZ [i][j];

			GV [i] += VISC [ In1 ] * RX [i][j] * RZ [i][j];

			HV [i] += VISC [ In1 ] * RX [i][j] * RY [i][j];
		}

		/////////////////////////////////////////////////////////////////

genw1:
		if ( IR1 [i] == 0 ) goto genw2;

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw1;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR1 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 1] = tmp1 * cos ( tmp2 );

genw2:
		if ( IR2 [i] == 0 ) goto genw3;

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw2;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR2 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 2] = tmp1 * cos ( tmp2 );

genw3:
		if ( IR3 [i] == 0 ) goto genvcms1;

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw3;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR3 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 3] = tmp1 * cos ( tmp2 );

genvcms1:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms1;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 1] = tmp1 * cos ( tmp2 );

genvcms2:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms2;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 2] = tmp1 * cos ( tmp2 );

genvcms3:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms3;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 3] = tmp1 * cos ( tmp2 );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void remove_crowder ( int i )
{
	int j, k, k2;

	int alist1_mass = 0;

	int * alist1_content = NULL, * alist2_content = NULL;

	/////////////////////////////////////////////////////////////////

	alist1_mass = NS_atom [i] + 1;

	alist1_content = (int *) calloc ( alist1_mass, sizeof(int) );

	alist2_content = (int *) calloc ( alist1_mass, sizeof(int) );

	/////////////////////////////////////////////////////////////////

	for (j = 0; j < alist1_mass; j++)
	{
		alist1_content [j] = atom_key [i][j];

		alist2_content [j] = NT_atom - j;
	}

	/////////////////////////////////////////////////////////////////

check_lists:

	for (j = 0; j < alist1_mass; j++)
	{
		for (k = 0; k < alist1_mass; k++)
		{
			if ( alist1_content [j] == alist2_content [k] )
			{
				alist1_content [j] = alist1_content [ alist1_mass - 1 ];

				alist2_content [k] = alist2_content [ alist1_mass - 1 ];

				alist1_mass -= 1;

				goto check_lists;
			}
		}
	}

	/////////////////////////////////////////////////////////////////

	for (j = 0; j < alist1_mass; j++)
	{
		k = alist1_content [j];

		k2 = alist2_content [j];

		/////////////////////////////////////////////////////////////////

		atom_key [ INDX [k2] ][ JNDX [k2] ] = k;

		part_key [k] = part_key [k2];

		maxi_key [k] = maxi_key [k2];

		INDX [k] = INDX [k2];

		JNDX [k] = JNDX [k2];

		/////////////////////////////////////////////////////////////////

		x [k] = x [k2];

		y [k] = y [k2];

		z [k] = z [k2];

		/////////////////////////////////////////////////////////////////

		x_old [k] = x_old [k2];

		y_old [k] = y_old [k2];

		z_old [k] = z_old [k2];

		/////////////////////////////////////////////////////////////////

		vx [k] = vx [k2];

		vy [k] = vy [k2];

		vz [k] = vz [k2];

		/////////////////////////////////////////////////////////////////

		fx [k] = fx [k2];

		fy [k] = fy [k2];

		fz [k] = fz [k2];
	}

	/////////////////////////////////////////////////////////////////

	free ( alist1_content );

	free ( alist2_content );

	/////////////////////////////////////////////////////////////////

	NT_atom -= (NS_atom [i] + 1);

	/////////////////////////////////////////////////////////////////

	part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

	maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );

	INDX = (int *) realloc ( INDX, (NT_atom + 1) * sizeof(int) );

	JNDX = (int *) realloc ( JNDX, (NT_atom + 1) * sizeof(int) );

	/////////////////////////////////////////////////////////////////

	x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );

	y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

	z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	x_old = (double *) realloc ( x_old, (NT_atom + 1) * sizeof(double) );

	y_old = (double *) realloc ( y_old, (NT_atom + 1) * sizeof(double) );

	z_old = (double *) realloc ( z_old, (NT_atom + 1) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	vx = (double *) realloc ( vx, (NT_atom + 1) * sizeof(double) );

	vy = (double *) realloc ( vy, (NT_atom + 1) * sizeof(double) );

	vz = (double *) realloc ( vz, (NT_atom + 1) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	fx = (double *) realloc ( fx, (NT_atom + 1) * sizeof(double) );

	fy = (double *) realloc ( fy, (NT_atom + 1) * sizeof(double) );

	fz = (double *) realloc ( fz, (NT_atom + 1) * sizeof(double) );

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////

	k = crowder_key [i];

	for (j = 1; j < crowder_mass [k] + 1; j++)
	{
		if ( crowder_content [k][j] == i ) crowder_content [k][j] = crowder_content [k][ crowder_mass [k] ];
	}

	crowder_mass [k] -= 1;

	crowder_content [k] = (int *) realloc ( crowder_content [k], (crowder_mass [k] + 1) * sizeof(int) );

	/////////////////////////////////////////////////////////////////

	k = crowder_key [ NB_atom ];

	for (j = 1; j < crowder_mass [k] + 1; j++)
	{
		if ( crowder_content [k][j] == NB_atom ) crowder_content [k][j] = i;
	}

	/////////////////////////////////////////////////////////////////

	crowder_key [i] = crowder_key [ NB_atom ];

	crowder_key = (int *) realloc ( crowder_key, NB_atom * sizeof(int) );

	/////////////////////////////////////////////////////////////////

	NS_atom [i] = NS_atom [ NB_atom ];

	NS_atom = (int *) realloc ( NS_atom, NB_atom * sizeof(int) );

	/////////////////////////////////////////////////////////////////

	RIGID_SET [i] = RIGID_SET [ NB_atom ];

	RIGID_SET = (int *) realloc ( RIGID_SET, NB_atom * sizeof(int) );

	/////////////////////////////////////////////////////////////////

	RIGID_END [i] = RIGID_END [ NB_atom ];

	RIGID_END = (int *) realloc ( RIGID_END, NB_atom * sizeof(int) );

	/////////////////////////////////////////////////////////////////

	IR1 [i] = IR1 [ NB_atom ];

	IR1 = (double *) realloc ( IR1, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	IR2 [i] = IR2 [ NB_atom ];

	IR2 = (double *) realloc ( IR2, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	IR3 [i] = IR3 [ NB_atom ];

	IR3 = (double *) realloc ( IR3, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	RMASS [i] = RMASS [ NB_atom ];

	RMASS = (double *) realloc ( RMASS, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	RVISC [i] = RVISC [ NB_atom ];

	RVISC = (double *) realloc ( RVISC, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	AV [i] = AV [ NB_atom ];

	AV = (double *) realloc ( AV, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	BV [i] = BV [ NB_atom ];

	BV = (double *) realloc ( BV, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	CV [i] = CV [ NB_atom ];

	CV = (double *) realloc ( CV, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	FV [i] = FV [ NB_atom ];

	FV = (double *) realloc ( FV, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	GV [i] = GV [ NB_atom ];

	GV = (double *) realloc ( GV, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	HV [i] = HV [ NB_atom ];

	HV = (double *) realloc ( HV, NB_atom * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	k = 3 * (i - 1);

	k2 = 3 * (NB_atom - 1);

	/////////////////////////////////////////////////////////////////

	CMS [k + 1] = CMS [k2 + 1];

	CMS [k + 2] = CMS [k2 + 2];

	CMS [k + 3] = CMS [k2 + 3];

	CMS = (double *) realloc ( CMS, (3 * NB_atom - 2) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	VCMS [k + 1] = VCMS [k2 + 1];

	VCMS [k + 2] = VCMS [k2 + 2];

	VCMS [k + 3] = VCMS [k2 + 3];

	VCMS = (double *) realloc ( VCMS, (3 * NB_atom - 2) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	W [k + 1] = W [k2 + 1];

	W [k + 2] = W [k2 + 2];

	W [k + 3] = W [k2 + 3];

	W = (double *) realloc ( W, (3 * NB_atom - 2) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	RXYZ [k + 1] = RXYZ [k2 + 1];

	RXYZ [k + 2] = RXYZ [k2 + 2];

	RXYZ [k + 3] = RXYZ [k2 + 3];

	RXYZ = (double *) realloc ( RXYZ, (3 * NB_atom - 2) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	k = 9 * (i - 1);

	k2 = 9 * (NB_atom - 1);

	/////////////////////////////////////////////////////////////////

	AXES [k + 1] = AXES [k2 + 1];

	AXES [k + 2] = AXES [k2 + 2];

	AXES [k + 3] = AXES [k2 + 3];

	AXES [k + 4] = AXES [k2 + 4];

	AXES [k + 5] = AXES [k2 + 5];

	AXES [k + 6] = AXES [k2 + 6];

	AXES [k + 7] = AXES [k2 + 7];

	AXES [k + 8] = AXES [k2 + 8];

	AXES [k + 9] = AXES [k2 + 9];

	AXES = (double *) realloc ( AXES, (9 * NB_atom - 8) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	free ( atom_key [i] );

	atom_key [i] = atom_key [ NB_atom ];

	atom_key = (int **) realloc ( atom_key, NB_atom * sizeof(int *) );

	/////////////////////////////////////////////////////////////////

	free ( RX [i] );

	RX [i] = RX [ NB_atom ];

	RX = (double **) realloc ( RX, NB_atom * sizeof(double *) );

	/////////////////////////////////////////////////////////////////

	free ( RY [i] );

	RY [i] = RY [ NB_atom ];

	RY = (double **) realloc ( RY, NB_atom * sizeof(double *) );

	/////////////////////////////////////////////////////////////////

	free ( RZ [i] );

	RZ [i] = RZ [ NB_atom ];

	RZ = (double **) realloc ( RZ, NB_atom * sizeof(double *) );

	/////////////////////////////////////////////////////////////////

	if (i < NB_atom)
	{
		for (j = 0; j < NS_atom [i] + 1; j++)
		{
			k = atom_key [i][j];

			INDX [k] = i;
		}
	}

	/////////////////////////////////////////////////////////////////

	NB_atom -= 1;
}



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

void treat_overlaps ( void )
{
	int i, j, k;

	int i2, j2, k2;

	int j0, k0;

	int m, n;

	int status, round;

	int In1, In2, In12;

	int J1, J2;

	int outlier_mass = 0, * outlier_content = NULL;

	double D_max, a0, a1, a2, xx [3], n1 [3];

	/////////////////////////////////////////////////////////////////

	D_max = D [0];

	for (i = 1; i < ns; i++) { if ( D [i] > D_max ) D_max = D [i]; }

	/////////////////////////////////////////////////

	for (i = NBP_atom + 1; i < NB_atom + 1; i++)
	{
		status = 0;

		/////////////////////////////////////////////////////////////////

		for (j = 0; j < NS_atom [i] + 1; j++)
		{
			k = atom_key [i][j];

			In1 = part_key [k];

			a2 = ( D [ In1 ] + D_max ) / 2;

			/////////////////////////////////////////////////////////////////

			if ( x [k] < a2 || x [k] > ( side_X - a2 ) ) { status = 1; break; }

			if ( y [k] < a2 || y [k] > ( side_Y - a2 ) ) { status = 1; break; }

			if ( z [k] < a2 || z [k] > ( side_Z - a2 ) ) { status = 1; break; }
		}

		/////////////////////////////////////////////////////////////////

		if ( status )
		{
			outlier_mass += 1;

			outlier_content = (int *) realloc ( outlier_content, (outlier_mass + 1) * sizeof(int) );

			outlier_content [ outlier_mass ] = i;
		}
	}

	/////////////////////////////////////////////////////////////////

	round = 1;

	status = 0;

	for (m = 1; m < outlier_mass + 1; m++)
	{
		i = outlier_content [m];

		J1 = 3 * (i - 1);

		/////////////////////////////////////////////////////////////////

		for (n = m + 1; n < outlier_mass + 1; n++)
		{
			i2 = outlier_content [n];

			J2 = 3 * (i2 - 1);

			/////////////////////////////////////////////////////////////////

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				In1 = part_key [k];

				/////////////////////////////////////////////////////////////////

				for (j2 = 0; j2 < NS_atom [i2] + 1; j2++)
				{
					k2 = atom_key [i2][j2];

					In2 = part_key [k2];

					/////////////////////////////////////////////////////////////////

					if (In1 < In2) In12 = ns * In1 - In1 * (In1 + 1) / 2 + In2;

					else In12 = ns * In2 - In2 * (In2 + 1) / 2 + In1;

					/////////////////////////////////////////////////////////////////

					xx [0] = x [k] - x [k2];

					xx [1] = y [k] - y [k2];

					xx [2] = z [k] - z [k2];

					half_shift ( xx );

					/////////////////////////////////////////////////////////////////

					a2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

					a1 = 0.65 * PD_sq [ In12 ] / a2;

					/////////////////////////////////////////////////////////////////

					if (a1 > 1)
					{
						status += 1;

						a0 = 0.5 * ( sqrt (a1) - 0.99 );

						//printf ( "%le\n", a2 / PD_sq [ In12 ] );

						/////////////////////////////////////////////////////////////////

						n1 [0] = xx [0] * a0;

						n1 [1] = xx [1] * a0;

						n1 [2] = xx [2] * a0;

						/////////////////////////////////////////////////////////////////

						if ( RIGID_SET [i] <= RIGID_END [i] )
						{
							CMS [J1 + 1] += n1 [0];

							CMS [J1 + 2] += n1 [1];

							CMS [J1 + 3] += n1 [2];
						}

						/////////////////////////////////////////////////////////////////

						for (j0 = 0; j0 < NS_atom [i] + 1; j0++)
						{
							k0 = atom_key [i][j0];

							/////////////////////////////////////////////////////////////////

							x [k0] += n1 [0];

							y [k0] += n1 [1];

							z [k0] += n1 [2];

							/////////////////////////////////////////////////////////////////

							x_old [k0] += n1 [0];

							y_old [k0] += n1 [1];

							z_old [k0] += n1 [2];
						}

						/////////////////////////////////////////////////////////////////

						if ( RIGID_SET [i2] <= RIGID_END [i2] )
						{
							CMS [J2 + 1] -= n1 [0];

							CMS [J2 + 2] -= n1 [1];

							CMS [J2 + 3] -= n1 [2];
						}

						/////////////////////////////////////////////////////////////////

						for (j0 = 0; j0 < NS_atom [i2] + 1; j0++)
						{
							k0 = atom_key [i2][j0];

							/////////////////////////////////////////////////////////////////

							x [k0] -= n1 [0];

							y [k0] -= n1 [1];

							z [k0] -= n1 [2];

							/////////////////////////////////////////////////////////////////

							x_old [k0] -= n1 [0];

							y_old [k0] -= n1 [1];

							z_old [k0] -= n1 [2];
						}
					}
				}
			}
		}
	}

	//printf ( "klok = %ld, round = %d, status = %d\n", klok, round, status );

	/////////////////////////////////////////////////////////////////

check_for_overlaps:

	round += 1;

	status = 0;

	for (m = 1; m < outlier_mass + 1; m++)
	{
		i = outlier_content [m];

		J1 = 3 * (i - 1);

		/////////////////////////////////////////////////////////////////

		for (i2 = 1; i2 < NB_atom + 1; i2++)
		{
			if (i2 == i) continue;

			/////////////////////////////////////////////////////////////////

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				In1 = part_key [k];

				/////////////////////////////////////////////////////////////////

				for (j2 = 0; j2 < NS_atom [i2] + 1; j2++)
				{
					k2 = atom_key [i2][j2];

					In2 = part_key [k2];

					/////////////////////////////////////////////////////////////////

					if (In1 < In2) In12 = ns * In1 - In1 * (In1 + 1) / 2 + In2;

					else In12 = ns * In2 - In2 * (In2 + 1) / 2 + In1;

					/////////////////////////////////////////////////////////////////

					xx [0] = x [k] - x [k2];

					xx [1] = y [k] - y [k2];

					xx [2] = z [k] - z [k2];

					half_shift ( xx );

					/////////////////////////////////////////////////////////////////

					a2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

					a1 = 0.65 * PD_sq [ In12 ] / a2;

					/////////////////////////////////////////////////////////////////

					if (a1 > 1)
					{
						status += 1;

						a0 = sqrt (a1) - 0.99;

						//printf ( "%le\n", a2 / PD_sq [ In12 ] );

						/////////////////////////////////////////////////////////////////

						generate_perpendicular_vector ( xx, n1 );

						n1 [0] = xx [0] * a0 + n1 [0] * 1.0;

						n1 [1] = xx [1] * a0 + n1 [1] * 1.0;

						n1 [2] = xx [2] * a0 + n1 [2] * 1.0;

						/////////////////////////////////////////////////////////////////

						if ( RIGID_SET [i] <= RIGID_END [i] )
						{
							CMS [J1 + 1] += n1 [0];

							CMS [J1 + 2] += n1 [1];

							CMS [J1 + 3] += n1 [2];
						}

						/////////////////////////////////////////////////////////////////

						for (j0 = 0; j0 < NS_atom [i] + 1; j0++)
						{
							k0 = atom_key [i][j0];

							/////////////////////////////////////////////////////////////////

							x [k0] += n1 [0];

							y [k0] += n1 [1];

							z [k0] += n1 [2];

							/////////////////////////////////////////////////////////////////

							x_old [k0] += n1 [0];

							y_old [k0] += n1 [1];

							z_old [k0] += n1 [2];
						}
					}
				}
			}
		}
	}

	//printf ( "klok = %ld, round = %d, status = %d\n", klok, round, status );

	if ( status && round < 100 ) goto check_for_overlaps;

	/////////////////////////////////////////////////////////////////

	for (m = 1; m < outlier_mass + 1; m++)
	{
		i = outlier_content [m];

		J1 = 3 * (i - 1);

		J2 = atom_key [i][0];

		////////////////////////////////////////////////////

		if ( x [J2] < 0 )
		{
			if ( RIGID_SET [i] <= RIGID_END [i] ) CMS [J1 + 1] += side_X;

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				x [k] += side_X;

				x_old [k] += side_X;
			}
		}

		////////////////////////////////////////////////////

		if ( x [J2] > side_X )
		{
			if ( RIGID_SET [i] <= RIGID_END [i] ) CMS [J1 + 1] -= side_X;

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				x [k] -= side_X;

				x_old [k] -= side_X;
			}
		}

		////////////////////////////////////////

		if ( y [J2] < 0 )
		{
			if ( RIGID_SET [i] <= RIGID_END [i] ) CMS [J1 + 2] += side_Y;

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				y [k] += side_Y;

				y_old [k] += side_Y;
			}
		}

		////////////////////////////////////////

		if ( y [J2] > side_Y )
		{
			if ( RIGID_SET [i] <= RIGID_END [i] ) CMS [J1 + 2] -= side_Y;

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				y [k] -= side_Y;

				y_old [k] -= side_Y;
			}
		}

		////////////////////////////////////////

		if ( z [J2] < 0 )
		{
			if ( RIGID_SET [i] <= RIGID_END [i] ) CMS [J1 + 3] += side_Z;

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				z [k] += side_Z;

				z_old [k] += side_Z;
			}
		}

		////////////////////////////////////////

		if ( z [J2] > side_Z )
		{
			if ( RIGID_SET [i] <= RIGID_END [i] ) CMS [J1 + 3] -= side_Z;

			for (j = 0; j < NS_atom [i] + 1; j++)
			{
				k = atom_key [i][j];

				z [k] -= side_Z;

				z_old [k] -= side_Z;
			}
		}
	}

	/////////////////////////////////////////////////////////////////

	if ( outlier_mass ) free ( outlier_content );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void adjust_box ( void )
{
	int i, k;

	int J1, J2;

	int In1, In2, In12, InInt;

	double x_min, y_min, z_min;

	double x_max, y_max, z_max;

	double A, B, C, F, G, H;

	double a0, a1, a2;

	double * ep1 = NULL, * ep2 = NULL, * ep3 = NULL;

	double n1 [3], n2 [3], n3 [3], xx [3];

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

	/////////////////////////////////

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

		if ( z_max < x_min ) { z_min = y_min; y_min = x_min; x_min = z_max; ep3 = ep2; ep2 = ep1; ep1 = n3; }

		else if ( z_max < y_min ) { z_min = y_min; y_min = z_max; ep3 = ep2; ep2 = n3; }

		else { z_min = z_max; ep3 = n3; }

		/////////////////////////////////

		if ( ep1 [0] < 0 ) { ep1 [0] = - ep1 [0]; ep1 [1] = - ep1 [1]; ep1 [2] = - ep1 [2]; }

		if ( ep2 [1] < 0 ) { ep2 [0] = - ep2 [0]; ep2 [1] = - ep2 [1]; ep2 [2] = - ep2 [2]; }
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

	///////////////Adjust box///////////////
	////////////////////////////////////////

	for (i = 1; i < NT_atom + 1; i++)
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

	for (i = 1; i < NT_atom + 1; i++)
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

	for (i = 1; i < NB_atom + 1; i++)
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

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////

	//printf ("\n\n");

	for (In1 = 0; In1 < n_crwd; In1++)
	{
		N_crwd [ In1 ] = (int) ( phi_crwd [ In1 ] * side_XYZ / V_crwd [ In1 ] );

		//printf ( "klok = %ld, crowder_mass [%d] = %d, N_crwd [%d] = %d\n", klok, In1, crowder_mass [ In1 ], In1, N_crwd [ In1 ]  );

		if ( crowder_mass [ In1 ] > N_crwd [ In1 ] )
		{
			do
			{
				a0 = 0;

				for (In2 = 1; In2 < crowder_mass [ In1 ] + 1; In2++)
				{
					i = crowder_content [ In1 ][ In2 ];

					k = atom_key [i][0];

					/////////////////////////////////////////////////////////////////

					xx [0] = x [k] - side_hX;

					xx [1] = y [k] - side_hY;

					xx [2] = z [k] - side_hZ;

					/////////////////////////////////////////////////////////////////

					a2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

					if (a2 > a0) { a0 = a2; J1 = i; }
				}

				remove_crowder ( J1 );

				//printf ( "Removed crowder %d\n", J1 );
			}
			while ( crowder_mass [ In1 ] > N_crwd [ In1 ] );
		}
	}

	/////////////////////////////////////////////////////////////////

	treat_overlaps ();

	/////////////////////////////////////////////////////////////////

	for (In1 = 0; In1 < n_crwd; In1++)
	{
		if ( crowder_mass [ In1 ] < N_crwd [ In1 ] )
		{
			do add_crowder ( In1 ); while ( crowder_mass [ In1 ] < N_crwd [ In1 ] );

			//printf ( "Inserted new crowders: crowder_mass [%d] = %d, N_crwd [%d] = %d\n", In1, crowder_mass [ In1 ], In1, N_crwd [ In1 ]  );
		}
	}

	/////////////////////////////////////////////////////////////////

	for (In1 = 0; In1 < ns; In1++)
	{
		for (In2 = In1; In2 < ns; In2++)
		{
			In12 = In1 * ns - In1 * (1 + In1) / 2 + In2;

			for (InInt = 0; InInt < npi; InInt++)
			{
				populate_list ( In1, In2, InInt );

				PROGRESS [ In12 ][ InInt ] = 0;

				//printf ("list_mass [%d][%d] = %d\n\n", In12, InInt, list_mass [In12][InInt] );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_overlaps ( void )
{
	int i, j, k;

	int i2, j2, k2;

	int i3;

	int In1, In2, In12;

	int outlier_mass = 0, * outlier_content = NULL;

	double a0, a1, xx [3];

	/////////////////////////////////////////////////

	 a1 = D [0];

	for (i = 1; i < ns; i++) { if ( D [i] > a1 ) a1 = D [i]; }

	/////////////////////////////////////////////////

	for (i = NBP_atom + 1; i < NB_atom + 1; i++)
	{
		i3 = 0;

		/////////////////////////////////////////////////////////////////

		for (j = 0; j < NS_atom [i] + 1; j++)
		{
			k = atom_key [i][j];

			In1 = part_key [k];

			a0 = ( D [ In1 ] + a1 ) / 2;

			/////////////////////////////////////////////////////////////////

			if ( x [k] < a0 || x [k] > ( side_X - a0 ) ) { i3 = 1; break; }

			if ( y [k] < a0 || y [k] > ( side_Y - a0 ) ) { i3 = 1; break; }

			if ( z [k] < a0 || z [k] > ( side_Z - a0 ) ) { i3 = 1; break; }
		}

		/////////////////////////////////////////////////////////////////

		if ( i3 )
		{
			outlier_mass += 1;

			outlier_content = (int *) realloc ( outlier_content, (outlier_mass + 1) * sizeof(int) );

			outlier_content [ outlier_mass ] = i;
		}
	}

	/////////////////////////////////////////////////////////////////

	a1 = 1.0;

	/////////////////////////////////////////////////////////////////

	for (i3 = 1; i3 < outlier_mass + 1; i3++)
	{
		i = outlier_content [i3];

		/////////////////////////////////////////////////////////////////

		for (j = 0; j < NS_atom [i] + 1; j++)
		{
			k = atom_key [i][j];

			In1 = part_key [k];

			for (i2 = 1; i2 < NB_atom + 1; i2++)
			{
				if (i2 == i) continue;

				for (j2 = 0; j2 < NS_atom [i2] + 1; j2++)
				{
					k2 = atom_key [i2][j2];

					In2 = part_key [k2];

					if (In1 < In2) In12 = ns * In1 - In1 * (In1 + 1) / 2 + In2;

					else In12 = ns * In2 - In2 * (In2 + 1) / 2 + In1;

					/////////////////////////////////////////////////////////////////

					xx [0] = x [k] - x [k2];

					xx [1] = y [k] - y [k2];

					xx [2] = z [k] - z [k2];

					half_shift ( xx );

					/////////////////////////////////////////////////////////////////

					a0 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

					/////////////////////////////////////////////////////////////////

					a0 = a0 / PD_sq [ In12 ];

					if ( a0 < a1 ) a1 = a0;
				}
			}
		}
	}

	printf ( "%le\n", a1 );

	/////////////////////////////////////////////////////////////////

	free ( outlier_content );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int check_box ( void )
{
	int i, status;

	/////////////////////////////////////

	status = 0;

	for (i = 1; i < NTP_atom + 1; i++)
	{
		if ( x [i] < BORDER_MIN ) { status = 1; break; }

		if ( x [i] > side_X - BORDER_MIN ) { status = 1; break; }

		if ( y [i] < BORDER_MIN ) { status = 1; break; }

		if ( y [i] > side_Y - BORDER_MIN ) { status = 1; break; }

		if ( z [i] < BORDER_MIN ) { status = 1; break; }

		if ( z [i] > side_Z - BORDER_MIN ) { status = 1; break; }
	}

	/////////////////////////////////////

	if ( status || (klok % 10000) == 0 )
	{
		adjust_box ();

		//printf ( "klok = %ld: box adjusted\n", klok );
	}

	return ( status );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Bfactor ( char * file_Bfactor )
{
	FILE * f1;

	long CNT;

	int i, k;

	double x_cm, y_cm, z_cm, r2;

	double * rB = NULL;

	double * rB2 = NULL;

	////////////////////////////////////////

	CNT = 0;

	x_cm = 0;

	y_cm = 0;

	z_cm = 0;

	for (i = 1; i < NTP_atom + 1; i++)
	{
		k = part_key [i];

		if (k != 3)
		{
			CNT += 1;

			x_cm += x [i];

			y_cm += y [i];

			z_cm += z [i];
		}
	}

	x_cm /= CNT;

	y_cm /= CNT;

	z_cm /= CNT;

	///////////////////////////////////////////////

	for (i = 1; i < NTP_atom + 1; i++)
	{
		k = part_key [i];

		if ( k != 3 )
		{
			r2 = (x [i] - x_cm) * (x [i] - x_cm) + (y [i] - y_cm) * (y [i] - y_cm) + (z [i] - z_cm) * (z [i] - z_cm);

			RB [i] += sqrt (r2);

			RB2 [i] += r2;
		}
	}

	///////////////////////////////////////////////

	COUNT_Bfactor += 1;

	///////////////////////////////////////////////

	if ( COUNT_Bfactor % 100000 == 0 )
	{
		rB = (double *) calloc ( NTP_atom + 1, sizeof(double) );

		rB2 = (double *) calloc ( NTP_atom + 1, sizeof(double) );

		////////////////////////////////////////

		f1 = fopen ( file_Bfactor, "r" );

		if (!f1) CNT = 0;

		else
		{
			fseek ( f1, 0L, SEEK_SET );

			fscanf ( f1, "%ld\n", &CNT );

			for (i = 1; i < NTP_atom + 1; i++)
			{
				k = part_key [i];

				if ( k != 3 ) fscanf ( f1, "%le %le\n", rB + i, rB2 + i );
			}

			fclose (f1);
		}

		///////////////////////////////////////////////

		f1 = fopen ( file_Bfactor, "w" );

		fseek ( f1, 0L, SEEK_SET );

		fprintf ( f1, "%ld\n", CNT + 1 );

		for (i = 1; i < NTP_atom + 1; i++)
		{
			k = part_key [i];

			if ( k != 3 )
			{
				rB [i] += RB [i] / COUNT_Bfactor;

				rB2 [i] += RB2 [i] / COUNT_Bfactor;

				//////////////////////////////////////////

				fprintf ( f1, "%.16le %.16le\n", rB [i], rB2 [i] );

				//////////////////////////////////////////

				RB [i] = 0;

				RB2 [i] = 0;
			}
		}

		fclose (f1);

		//////////////////////////////////////////

		COUNT_Bfactor = 0;

		//////////////////////////////////////////

		free ( rB );

		free ( rB2 );
	}
}
