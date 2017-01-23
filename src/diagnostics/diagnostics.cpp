#include "diagnostics.h"

void velocity_rescale ( double * sc )
{
	int i, j, k, l, J1, nV, nVCM, nW;

	double T0, T1, T2;

	////////////////////////////

	T0 = 0; nV = 0;

	T1 = 0; nVCM = 0;

	T2 = 0; nW = 0;

	for (i = 1; i < NB_atom + 1; i++)
	{
		for (j = 0; j < RIGID_SET [i]; j++)
		{
			nV += 3;

			k = atom_key [i][j];

			l = part_key [k];

			T0 += MASS [l] * ( vx [k] * vx [k] + vy [k] * vy [k] + vz [k] * vz [k] );
		}

		for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++)
		{
			nV += 3;

			k = atom_key [i][j];

			l = part_key [k];

			T0 += MASS [l] * ( vx [k] * vx [k] + vy [k] * vy [k] + vz [k] * vz [k] );
		}

		if ( RIGID_SET [i] > RIGID_END [i] ) continue;



		///////////////////////////////////////////////////////////////////////////////////////////////

		J1 = 3 * (i - 1);

		////////////////////////////

		nVCM += 3;

		T1 += RMASS [i] * VCMS [J1 + 1] * VCMS [J1 + 1];

		T1 += RMASS [i] * VCMS [J1 + 2] * VCMS [J1 + 2];

		T1 += RMASS [i] * VCMS [J1 + 3] * VCMS [J1 + 3];



		///////////////////////////////////////////////////////////////////////////////////////////////

		if ( IR1 [i] == 0 && IR2 [i] == 0 ) continue;

		else if ( IR1 [i] == 0 )
		{
			nW += 2;

			T2 += IR2 [i] * W [J1 + 2] * W [J1 + 2];

			T2 += IR3 [i] * W [J1 + 3] * W [J1 + 3];
		}

		else
		{
			nW += 3;

			T2 += IR1 [i] * W [J1 + 1] * W [J1 + 1];

			T2 += IR2 [i] * W [J1 + 2] * W [J1 + 2];

			T2 += IR3 [i] * W [J1 + 3] * W [J1 + 3];
		}
	}

	////////////////////////////

	sc [0] = sqrt ( nV * T / T0 );

	sc [1] = sqrt ( nVCM * T / T1 );

	sc [2] = sqrt ( nW * T / T2 );

	sc [3] = sqrt ( ( nV + nVCM + nW ) * T / ( T0 + T1 + T2 ) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double check_maxwell ( void )
{
	int i, j, k, l, nV, count;

	double vsq;

	////////////////////////////////

	nV = 0; count = 0;

	for (i = 1; i < NB_atom + 1; i++)
	{
		for (j = 0; j < RIGID_SET [i]; j++)
		{
			nV += 3;

			k = atom_key [i][j];

			l = part_key [k];

			////////////////////////////////

			vsq = MASS [l] * vx [k] * vx [k];

			if ( vsq > T ) count += 1;

			////////////////////////////////

			vsq = MASS [l] * vy [k] * vy [k];

			if ( vsq > T ) count += 1;

			////////////////////////////////

			vsq = MASS [l] * vz [k] * vz [k];

			if ( vsq > T ) count += 1;
		}

		////////////////////////////////

		for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++)
		{
			nV += 3;

			k = atom_key [i][j];

			l = part_key [k];

			////////////////////////////////

			vsq = MASS [l] * vx [k] * vx [k];

			if ( vsq > T ) count += 1;

			////////////////////////////////

			vsq = MASS [l] * vy [k] * vy [k];

			if ( vsq > T ) count += 1;

			////////////////////////////////

			vsq = MASS [l] * vz [k] * vz [k];

			if ( vsq > T ) count += 1;
		}
	}

	////////////////////////////////

	vsq = (double) count * 100 / nV;

	return ( vsq );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double check_cm_maxwell ( void )
{
	int i, J1, nVCM, count;

	double vsq;

	////////////////////////////////

	nVCM = 0; count = 0;

	for (i = 1; i < NB_atom + 1; i++)
	{
		if ( RIGID_SET [i] > RIGID_END [i] ) continue;

		////////////////////////////////

		J1 = 3 * (i - 1);

		////////////////////////////////

		nVCM += 3;

		////////////////////////////////

		vsq = RMASS [i] * VCMS [J1 + 1] * VCMS [J1 + 1];

		if ( vsq > T ) count += 1;

		////////////////////////////////

		vsq = RMASS [i] * VCMS [J1 + 2] * VCMS [J1 + 2];

		if ( vsq > T ) count += 1;

		////////////////////////////////

		vsq = RMASS [i] * VCMS [J1 + 3] * VCMS [J1 + 3];

		if ( vsq > T ) count += 1;
	}

	vsq = (double) count * 100 / nVCM;

	return ( vsq );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double check_angular_maxwell ( void )
{
	int i, J1, nW, count;

	double wsq;

	////////////////////////////////

	nW = 0; count = 0;

	for (i = 1; i < NB_atom + 1; i++)
	{
		if ( RIGID_SET [i] > RIGID_END [i] ) continue;

		////////////////////////////////

		J1 = 3 * (i - 1);

		////////////////////////////////

		if ( IR1 [i] == 0 && IR2 [i] == 0 ) continue;

		////////////////////////////////

		else if ( IR1 [i] == 0 )
		{
			nW += 2;

			////////////////////////////////

			wsq = IR2 [i] * W [J1 + 2] * W [J1 + 2];

			if ( wsq > T ) count += 1;

			////////////////////////////////

			wsq = IR3 [i] * W [J1 + 3] * W [J1 + 3];

			if ( wsq > T ) count += 1;
		}

		else
		{
			nW += 3;

			////////////////////////////////

			wsq = IR1 [i] * W [J1 + 1] * W [J1 + 1];

			if ( wsq > T ) count += 1;

			////////////////////////////////

			wsq = IR2 [i] * W [J1 + 2] * W [J1 + 2];

			if ( wsq > T ) count += 1;

			////////////////////////////////

			wsq = IR3 [i] * W [J1 + 3] * W [J1 + 3];

			if ( wsq > T ) count += 1;
		}
	}

	wsq = (double) count * 100 / nW;

	return ( wsq );
}
