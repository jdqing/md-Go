#include "move.h"

void move_rigid_units ( void )
{
	int i, j, k, J1, J2;

	double tmp1, tmp2, tmp3;

	double v1, v2, v3;

	double w1, w2, w3;

	double vec1 [3], vec2 [3], vec3 [3];

	double FC1, FC2, FC3;

	double F1, F2, F3;

	double TC1, TC2, TC3;

	double T1, T2, T3;

	//////////////////////////////////////////

	for (i = 1; i < NB_atom + 1; i++)
	{
		for (j = 0; j < RIGID_SET [i]; j++) leap_frog_move_atom ( i, j );

		for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++) leap_frog_move_atom ( i, j );

		if ( RIGID_SET [i] > RIGID_END [i] ) continue;

		/////////////////////////////////////////////////////////////////

		J1 = 3 * (i - 1);

		J2 = 9 * (i - 1);

		//////////////////////////////////////////

		v1 = VCMS [J1 + 1];

		v2 = VCMS [J1 + 2];

		v3 = VCMS [J1 + 3];

		/////////////////////////////////

		if ( IR1 [i] == 0 && IR2 [i] == 0 )
		{
			k = atom_key [i][0];

			//////////////////////////////////////////

			F1 = ( fx [k] - RVISC [i] * v1 ) / RMASS [i];

			F2 = ( fy [k] - RVISC [i] * v2 ) / RMASS [i];

			F3 = ( fz [k] - RVISC [i] * v3 ) / RMASS [i];

			//////////////////////////////////////////

			VCMS [J1 + 1] += h * F1;

			VCMS [J1 + 2] += h * F2;

			VCMS [J1 + 3] += h * F3;

			//////////////////////////////////////////

			CMS [J1 + 1] += h * v1 + hsq2 * F1;

			CMS [J1 + 2] += h * v2 + hsq2 * F2;

			CMS [J1 + 3] += h * v3 + hsq2 * F3;

			//////////////////////////////////////

			x_old [k] = x [k];

			y_old [k] = y [k];

			z_old [k] = z [k];

			//////////////////////////////////////

			x [k] = CMS [J1 + 1];

			y [k] = CMS [J1 + 2];

			z [k] = CMS [J1 + 3];
		}

		/////////////////////////////////

		else if ( IR1 [i] == 0 )
		{
			w2 = W [J1 + 2];

			w3 = W [J1 + 3];

			//////////////////////////////////////////

			vec1 [0] = AXES [J2 + 1]; vec1 [1] = AXES [J2 + 2]; vec1 [2] = AXES [J2 + 3];

			vec2 [0] = AXES [J2 + 4]; vec2 [1] = AXES [J2 + 5]; vec2 [2] = AXES [J2 + 6];

			vec3 [0] = AXES [J2 + 7]; vec3 [1] = AXES [J2 + 8]; vec3 [2] = AXES [J2 + 9];

			//Calculate conservative forces in the box coordinate system
			//and conservative torques in the inertia coordinate system

			FC1 = 0; FC2 = 0; FC3 = 0;

			TC1 = 0; TC2 = 0; TC3 = 0;

			for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
			{
				k = atom_key [i][j];

				/////////////////////////////////////

				FC1 += fx [k];

				FC2 += fy [k];

				FC3 += fz [k];

				/////////////////////////////////////

				TC2 -= RX [i][j] * ( fx [k] * vec3 [0] + fy [k] * vec3 [1] + fz [k] * vec3 [2] );

				TC3 += RX [i][j] * ( fx [k] * vec2 [0] + fy [k] * vec2 [1] + fz [k] * vec2 [2] );
			}

			//Step1

			//Calculate full forces in the box coordinate system

			F2 = w3 * RXYZ [J1 + 1];

			F3 = - w2 * RXYZ [J1 + 1];

			//////////////////////////////////////

			tmp1 = F2 * vec2 [0] + F3 * vec3 [0];

			tmp2 = F2 * vec2 [1] + F3 * vec3 [1];

			tmp3 = F2 * vec2 [2] + F3 * vec3 [2];

			//////////////////////////////////////

			F1 = ( FC1 - RVISC [i] * v1 - tmp1 ) / RMASS [i];

			F2 = ( FC2 - RVISC [i] * v2 - tmp2 ) / RMASS [i];

			F3 = ( FC3 - RVISC [i] * v3 - tmp3 ) / RMASS [i];

			//Calculate full torques in the inertia coordinate system

			tmp2 = v1 * vec2 [0] + v2 * vec2 [1] + v3 * vec2 [2];

			tmp3 = v1 * vec3 [0] + v2 * vec3 [1] + v3 * vec3 [2];

			//////////////////////////////////////

			T2 = ( TC2 + RXYZ [J1 + 1] * tmp3 - BV [i] * w2 ) / IR2 [i];

			T3 = ( TC3 - RXYZ [J1 + 1] * tmp2 - CV [i] * w3 ) / IR3 [i];

			//Step2///////////////////////////////

			//Find new velocities in the box coordinate system

			VCMS [J1 + 1] += h * F1;

			VCMS [J1 + 2] += h * F2;

			VCMS [J1 + 3] += h * F3;

			//Find new angular velocities in the inertia coordinate system

			W [J1 + 2] += h * T2;

			W [J1 + 3] += h * T3;

			//Step3///////////////////////////////

			//Calcultate new positions in the box coordinate system

			CMS [J1 + 1] += h * v1 + hsq2 * F1;

			CMS [J1 + 2] += h * v2 + hsq2 * F2;

			CMS [J1 + 3] += h * v3 + hsq2 * F3;

			//Calcultate new vectors of the inertia coordinate system

			F1 = w3 * vec2 [0] - w2 * vec3 [0];

			F2 = w3 * vec2 [1] - w2 * vec3 [1];

			F3 = w3 * vec2 [2] - w2 * vec3 [2];

			/////////////////////////////////////////

			tmp1 = - (w2 * w2 + w3 * w3);

			tmp2 = T3;

			tmp3 = - T2;

			//////////////////////////////////////////

			v1 = vec1 [0] * tmp1 + vec2 [0] * tmp2 + vec3 [0] * tmp3;

			v2 = vec1 [1] * tmp1 + vec2 [1] * tmp2 + vec3 [1] * tmp3;

			v3 = vec1 [2] * tmp1 + vec2 [2] * tmp2 + vec3 [2] * tmp3;

			////////////////////////////////////////////////////

			FC1 = vec1 [0] + h * F1 + hsq2 * v1;

			FC2 = vec1 [1] + h * F2 + hsq2 * v2;

			FC3 = vec1 [2] + h * F3 + hsq2 * v3;

			//////////////////////////////////////////
			//////////////////////////////////////////

			F1 = - w3 * vec1 [0];

			F2 = - w3 * vec1 [1];

			F3 = - w3 * vec1 [2];

			//////////////////////////////////////////

			tmp1 = - T3;

			tmp2 = - w3 * w3;

			tmp3 = w2 * w3;

			//////////////////////////////////////////

			v1 = vec1 [0] * tmp1 + vec2 [0] * tmp2 + vec3 [0] * tmp3;

			v2 = vec1 [1] * tmp1 + vec2 [1] * tmp2 + vec3 [1] * tmp3;

			v3 = vec1 [2] * tmp1 + vec2 [2] * tmp2 + vec3 [2] * tmp3;

			///////////////////////////////////////////////////

			TC1 = vec2 [0] + h * F1 + hsq2 * v1;

			TC2 = vec2 [1] + h * F2 + hsq2 * v2;

			TC3 = vec2 [2] + h * F3 + hsq2 * v3;

			////////////////////////////////////////////////

			tmp1 = sqrt ( FC1 * FC1 + FC2 * FC2 + FC3 * FC3 );

			FC1 /= tmp1; FC2 /= tmp1; FC3 /= tmp1;

			/////////////////////////////////////////////////

			tmp1 = FC1 * TC1 + FC2 * TC2 + FC3 * TC3;

			TC1 -= tmp1 * FC1; TC2 -= tmp1 * FC2; TC3 -= tmp1 * FC3;

			/////////////////////////////////////////////////

			tmp1 = sqrt ( TC1 * TC1 + TC2 * TC2 + TC3 * TC3 );

			TC1 /= tmp1; TC2 /= tmp1; TC3 /= tmp1;

			////////////////////////////////////////////////////

			v1 = FC2 * TC3 - FC3 * TC2;

			v2 = FC3 * TC1 - FC1 * TC3;

			v3 = FC1 * TC2 - FC2 * TC1;

			////////////////////////////////////////////////////

			AXES [J2 + 1] = FC1; AXES [J2 + 2] = FC2; AXES [J2 + 3] = FC3;

			AXES [J2 + 4] = TC1; AXES [J2 + 5] = TC2; AXES [J2 + 6] = TC3;

			AXES [J2 + 7] = v1; AXES [J2 + 8] = v2; AXES [J2 + 9] = v3;

			///////////////////////////////////////////////////

			for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
			{
				k = atom_key [i][j];

				/////////////////////////////////////

				x_old [k] = x [k];

				y_old [k] = y [k];

				z_old [k] = z [k];

				//////////////////////////////////////

				x [k] = CMS [J1 + 1] + RX [i][j] * FC1;

				y [k] = CMS [J1 + 2] + RX [i][j] * FC2;

				z [k] = CMS [J1 + 3] + RX [i][j] * FC3;
			}
		}

		//////////////////////////////////////////

		else
		{
			w1 = W [J1 + 1];

			w2 = W [J1 + 2];

			w3 = W [J1 + 3];

			//////////////////////////////////////////

			vec1 [0] = AXES [J2 + 1]; vec1 [1] = AXES [J2 + 2]; vec1 [2] = AXES [J2 + 3];

			vec2 [0] = AXES [J2 + 4]; vec2 [1] = AXES [J2 + 5]; vec2 [2] = AXES [J2 + 6];

			vec3 [0] = AXES [J2 + 7]; vec3 [1] = AXES [J2 + 8]; vec3 [2] = AXES [J2 + 9];

			//Calculate conservative forces in the box coordinate system
			//and conservative torques in the inertia coordinate system

			FC1 = 0; FC2 = 0; FC3 = 0;

			T1 = 0; T2 = 0; T3 = 0;

			for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
			{
				k = atom_key [i][j];

				/////////////////////////////////////

				FC1 += fx [k];

				FC2 += fy [k];

				FC3 += fz [k];

				/////////////////////////////////////

				tmp1 = x [k] - CMS [J1 + 1];

				tmp2 = y [k] - CMS [J1 + 2];

				tmp3 = z [k] - CMS [J1 + 3];

				/////////////////////////////////////

				T1 += ( tmp2 * fz [k] - tmp3 * fy [k] );

				T2 += ( tmp3 * fx [k] - tmp1 * fz [k] );

				T3 += ( tmp1 * fy [k] - tmp2 * fx [k] );
			}

			//////////////////////////////////////

			TC1 = T1 * vec1 [0] + T2 * vec1 [1] + T3 * vec1 [2];

			TC2 = T1 * vec2 [0] + T2 * vec2 [1] + T3 * vec2 [2];

			TC3 = T1 * vec3 [0] + T2 * vec3 [1] + T3 * vec3 [2];

			//Step1

			//Calculate full forces in the box coordinate system

			F1 = w2 * RXYZ [J1 + 3] - w3 * RXYZ [J1 + 2];

			F2 = w3 * RXYZ [J1 + 1] - w1 * RXYZ [J1 + 3];

			F3 = w1 * RXYZ [J1 + 2] - w2 * RXYZ [J1 + 1];

			//////////////////////////////////////

			tmp1 = F1 * vec1 [0] + F2 * vec2 [0] + F3 * vec3 [0];

			tmp2 = F1 * vec1 [1] + F2 * vec2 [1] + F3 * vec3 [1];

			tmp3 = F1 * vec1 [2] + F2 * vec2 [2] + F3 * vec3 [2];

			//////////////////////////////////////

			F1 = ( FC1 - RVISC [i] * v1 - tmp1 ) / RMASS [i];

			F2 = ( FC2 - RVISC [i] * v2 - tmp2 ) / RMASS [i];

			F3 = ( FC3 - RVISC [i] * v3 - tmp3 ) / RMASS [i];

			//Calculate full torques in the inertia coordinate system

			tmp1 = v1 * vec1 [0] + v2 * vec1 [1] + v3 * vec1 [2];

			tmp2 = v1 * vec2 [0] + v2 * vec2 [1] + v3 * vec2 [2];

			tmp3 = v1 * vec3 [0] + v2 * vec3 [1] + v3 * vec3 [2];

			//////////////////////////////////////

			T1 = TC1 - ( RXYZ [J1 + 2] * tmp3 - RXYZ [J1 + 3] * tmp2 );

			T2 = TC2 - ( RXYZ [J1 + 3] * tmp1 - RXYZ [J1 + 1] * tmp3 );

			T3 = TC3 - ( RXYZ [J1 + 1] * tmp2 - RXYZ [J1 + 2] * tmp1 );

			//////////////////////////////////////

			T1 = ( T1 - AV [i] * w1 + HV [i] * w2 + GV [i] * w3 + (IR2 [i] - IR3 [i]) * w2 * w3 ) / IR1 [i];

			T2 = ( T2 + HV [i] * w1 - BV [i] * w2 + FV [i] * w3 + (IR3 [i] - IR1 [i]) * w1 * w3 ) / IR2 [i];

			T3 = ( T3 + GV [i] * w1 + FV [i] * w2 - CV [i] * w3 + (IR1 [i] - IR2 [i]) * w1 * w2 ) / IR3 [i];

			//Step2///////////////////////////////

			//Find new velocities in the box coordinate system

			VCMS [J1 + 1] += h * F1;

			VCMS [J1 + 2] += h * F2;

			VCMS [J1 + 3] += h * F3;

			//Find new angular velocities in the inertia coordinate system

			W [J1 + 1] += h * T1;

			W [J1 + 2] += h * T2;

			W [J1 + 3] += h * T3;

			//Step3///////////////////////////////

			//Calcultate new positions in the box coordinate system

			CMS [J1 + 1] += h * v1 + hsq2 * F1;

			CMS [J1 + 2] += h * v2 + hsq2 * F2;

			CMS [J1 + 3] += h * v3 + hsq2 * F3;

			//Calcultate new vectors of the inertia coordinate system

			F1 = w3 * vec2 [0] - w2 * vec3 [0];

			F2 = w3 * vec2 [1] - w2 * vec3 [1];

			F3 = w3 * vec2 [2] - w2 * vec3 [2];

			/////////////////////////////////////////

			tmp1 = - (w2 * w2 + w3 * w3);

			tmp2 = w1 * w2 + T3;

			tmp3 = w1 * w3 - T2;

			//////////////////////////////////////////

			v1 = vec1 [0] * tmp1 + vec2 [0] * tmp2 + vec3 [0] * tmp3;

			v2 = vec1 [1] * tmp1 + vec2 [1] * tmp2 + vec3 [1] * tmp3;

			v3 = vec1 [2] * tmp1 + vec2 [2] * tmp2 + vec3 [2] * tmp3;

			////////////////////////////////////////////////////

			FC1 = vec1 [0] + h * F1 + hsq2 * v1;

			FC2 = vec1 [1] + h * F2 + hsq2 * v2;

			FC3 = vec1 [2] + h * F3 + hsq2 * v3;

			//////////////////////////////////////////
			//////////////////////////////////////////

			F1 = w1 * vec3 [0] - w3 * vec1 [0];

			F2 = w1 * vec3 [1] - w3 * vec1 [1];

			F3 = w1 * vec3 [2] - w3 * vec1 [2];

			//////////////////////////////////////////

			tmp1 = w1 * w2 - T3;

			tmp2 = - (w1 * w1 + w3 * w3);

			tmp3 = w2 * w3 + T1;

			//////////////////////////////////////////

			v1 = vec1 [0] * tmp1 + vec2 [0] * tmp2 + vec3 [0] * tmp3;

			v2 = vec1 [1] * tmp1 + vec2 [1] * tmp2 + vec3 [1] * tmp3;

			v3 = vec1 [2] * tmp1 + vec2 [2] * tmp2 + vec3 [2] * tmp3;

			///////////////////////////////////////////////////

			TC1 = vec2 [0] + h * F1 + hsq2 * v1;

			TC2 = vec2 [1] + h * F2 + hsq2 * v2;

			TC3 = vec2 [2] + h * F3 + hsq2 * v3;

			////////////////////////////////////////////////

			tmp1 = sqrt ( FC1 * FC1 + FC2 * FC2 + FC3 * FC3 );

			FC1 /= tmp1; FC2 /= tmp1; FC3 /= tmp1;

			/////////////////////////////////////////////////

			tmp1 = FC1 * TC1 + FC2 * TC2 + FC3 * TC3;

			TC1 -= tmp1 * FC1; TC2 -= tmp1 * FC2; TC3 -= tmp1 * FC3;

			/////////////////////////////////////////////////

			tmp1 = sqrt ( TC1 * TC1 + TC2 * TC2 + TC3 * TC3 );

			TC1 /= tmp1; TC2 /= tmp1; TC3 /= tmp1;

			////////////////////////////////////////////////////

			v1 = FC2 * TC3 - FC3 * TC2;

			v2 = FC3 * TC1 - FC1 * TC3;

			v3 = FC1 * TC2 - FC2 * TC1;

			////////////////////////////////////////////////////

			AXES [J2 + 1] = FC1; AXES [J2 + 2] = FC2; AXES [J2 + 3] = FC3;

			AXES [J2 + 4] = TC1; AXES [J2 + 5] = TC2; AXES [J2 + 6] = TC3;

			AXES [J2 + 7] = v1; AXES [J2 + 8] = v2; AXES [J2 + 9] = v3;

			////////////////////////////////////////////////////

			for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
			{
				k = atom_key [i][j];

				////////////////////////////////////////////////////

				x_old [k] = x [k];

				y_old [k] = y [k];

				z_old [k] = z [k];

				////////////////////////////////////////////////////

				x [k] = CMS [J1 + 1] + RX [i][j] * FC1 + RY [i][j] * TC1 + RZ [i][j] * v1;

				y [k] = CMS [J1 + 2] + RX [i][j] * FC2 + RY [i][j] * TC2 + RZ [i][j] * v2;

				z [k] = CMS [J1 + 3] + RX [i][j] * FC3 + RY [i][j] * TC3 + RZ [i][j] * v3;
			}
		}
	}

	/////////////////////////////////////////////////////////////////

	//Periodic boundary for crowders

	for (i = NBP_atom + 1; i < NB_atom + 1; i++)
	{
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
}
