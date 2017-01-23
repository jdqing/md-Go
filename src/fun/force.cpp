#include "force.h"

int force_init()
{
  int In1, In2, In12, InInt;

  for (In1 = 0; In1 < ns; In1++)          // QM: ns=6.
	{
		for (In2 = In1; In2 < ns; In2++)
		{
			In12 = ns * In1 - In1 * (In1 + 1) / 2 + In2;

			PD [ In12 ] = ( D [ In1 ] + D [ In2 ] ) / 2;

			PD_sq [ In12 ] = PD [ In12 ] * PD [ In12 ];
		}

		K1 [ In1 ] = exp ( - VISC [ In1 ] / MASS [ In1 ] * h );

		K2 [ In1 ] = ( 1.0 - K1 [ In1 ] ) / VISC [ In1 ];

		SIGMA_FORCE [ In1 ] = sqrt ( 2.0 * VISC [ In1 ] * T / h );

    // printf("%lf\n",SIGMA_FORCE[ In1 ]);
	}

}

void CC_force_single_pair ( int i1, int i2, int Jn12, double scale_factor,
                            double ( * force) ( double x2, int Jn12 ) )
{
	double r2, xx [3];

	/////////////////////////////////////////////////////////////////

	xx [0] = x [i1] - x [i2];

	xx [1] = y [i1] - y [i2];

	xx [2] = z [i1] - z [i2];

	half_shift (xx);

	/////////////////////////////////////////////////////////////////

	r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

	r2 = scale_factor * force ( r2, Jn12 );

	/////////////////////////////////////////////////////////////////

	xx [0] = xx [0] * r2;

	xx [1] = xx [1] * r2;

	xx [2] = xx [2] * r2;

	/////////////////////////////////////////////////////////////////

	fx [i1] += xx [0]; fx [i2] -= xx [0];

	fy [i1] += xx [1]; fy [i2] -= xx [1];

	fz [i1] += xx [2]; fz [i2] -= xx [2];
}

double LJ_FORCE ( double x2, int Jn12 )
{
	if ( x2 > D2_LJ_CUTOFF [ Jn12 ] ) return (0);

	else
	{
		double r2, r4, r8, r14;

        r2 = D2_LJ [ Jn12 ] / x2;

		r4 = r2 * r2; r8 = r4 * r4; r14 = r8 * r4 * r2;

		r2 = F_LJ [ Jn12 ] * (r14 - r8);

		return ( r2 );
	}
}

double LJ_FORCE_n ( double x2, int Jn12 )
{
	if ( x2 > D2_LJ_CUTOFF_n [ Jn12 ] ) return (0);

	else
	{
		double r2, r4, r8, r14;

		r2 = D2_LJ_n [ Jn12 ] / x2;

		r4 = r2 * r2; r8 = r4 * r4; r14 = r8 * r4 * r2;

		r2 = r14 - r8;

		return ( r2 );
	}
}

void LJ_interactions_lists ( void )
{
	int i, j, k, In12, i1, i2, j1;

	///////////////////////////////////////////////////////

	for (In12 = 0; In12 < nps; In12++)
	{
		for (k = 1; k < list_mass [ In12 ][ 0 ] + 1; k += 2)
		{
			i = list_content [ In12 ][ 0 ][ k ];

			j = list_content [ In12 ][ 0 ][ k + 1 ];

			if ((Flags2Atom[i][j]=='b')||(Flags2Atom[i][j]=='a')||(Flags2Atom[i][j]=='c'))
                continue;

			///////////////////////////////////////////////////////

			i1 = maxi_key [i];

			i2 = maxi_key [j];

			if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

			else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

			/////////////////////////////////////////////////////////////////

			CC_force_single_pair ( i, j, j1, 1.0, LJ_FORCE );
		}
	}
}

void bond_force ( int i1, int i2, double ( * force) ( double r ) )
{
	double v1 [3], r1, tmp1;

	/////////////////////////////////////////////////////////////////

	v1 [0] = x [i2] - x [i1];

	v1 [1] = y [i2] - y [i1];

	v1 [2] = z [i2] - z [i1];

	/////////////////////////////////////////////////////////////////

	r1 = sqrt ( v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2] );

	tmp1 = force ( r1 ) / r1;

	/////////////////////////////////////////////////////////////////

	r1 = tmp1 * v1 [0]; fx [i1] -= r1; fx [i2] += r1;

	r1 = tmp1 * v1 [1]; fy [i1] -= r1; fy [i2] += r1;

	r1 = tmp1 * v1 [2]; fz [i1] -= r1; fz [i2] += r1;

	/////////////////////////////////////////////////////////////////

}

void valence_force ( int i1, int i2, int i3, double ( * force) ( double phi ) )
{
	double  v1 [3], v2 [3], v1_norm, v2_norm, v12_norm, kosinus, sinus, theta,
          tmp1, tmp2, tmp3, C1, C2;

	int k, l, I [4];

	/////////////////////////////////////////////////////////////////

	I [1] = i1; I [2] = i2; I [3] = i3;

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

	v1_norm = 1.0 / v1_norm;

	/////////////////////////////////////////////////////////////////

	v2_norm = v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2];

	v2_norm = 1.0 / v2_norm;

	/////////////////////////////////////////////////////////////////

	v12_norm = sqrt ( v1_norm * v2_norm );

	/////////////////////////////////////////////////////////////////

	kosinus = ( v1 [0] * v2 [0] + v1 [1] * v2 [1] + v1 [2] * v2 [2] ) * v12_norm;

	sinus =  sqrt ( 1.0 - kosinus * kosinus );

	theta = acos ( kosinus );

	/////////////////////////////////////////////////////////////////

	v1_norm = kosinus * v1_norm;

	v2_norm = kosinus * v2_norm;

	sinus = - force ( theta ) / sinus;

	/////////////////////////////////////////////////////////////////

	for (k = 0; k <= 2; k++)
	{
		for (l = 1; l <= 3; l++)
		{
			tmp1 = Kdelta (l, 1);

			tmp2 = Kdelta (l, 2);

			tmp3 = Kdelta (l, 3);

			/////////////////////////////////////////////////////////////////

			C1 = tmp1 - tmp2;

			C2 = tmp3 - tmp2;

			/////////////////////////////////////////////////////////////////

			tmp1 = C1 * v2 [k] + C2 * v1 [k];

			tmp2 = C1 * v1 [k];

			tmp3 = C2 * v2 [k];

			/////////////////////////////////////////////////////////////////

			C1 = tmp1 * v12_norm - tmp2 * v1_norm - tmp3 * v2_norm;

			/////////////////////////////////////////////////////////////////

			tmp1 = sinus * C1;

			if (k == 0) fx [ I [l] ] += tmp1;

			if (k == 1) fy [ I [l] ] += tmp1;

			if (k == 2) fz [ I [l] ] += tmp1;
		}
	}

	/////////////////////////////////////////////////////////////////

}

void dihedral_force ( int i1, int i2, int i3, int i4, double ( * force) ( double phi ) )
{
	double v1 [3], v2 [3], v3 [3], m [3], n [3], m_norm, n_norm, mn_norm, kosinus, sinus, psi, tmp1, tmp2, tmp3, tmp4;

	double A11, A22, A33, A12, A23, A13, C1, C2, D1, D2, D3;

	int k, l, I [5];

	/////////////////////////////////////////////////////////////////

	I [1] = i1; I [2] = i2; I [3] = i3; I [4] = i4;

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

	m_norm = m [0] * m [0] + m [1] * m [1] + m [2] * m [2];

	m_norm = 1.0 / m_norm;

	/////////////////////////////////////////////////////////////////

	n_norm = n [0] * n [0] + n [1] * n [1] + n [2] * n [2];

	n_norm = 1.0 / n_norm;

	/////////////////////////////////////////////////////////////////

	mn_norm = sqrt ( m_norm * n_norm );

	/////////////////////////////////////////////////////////////////

	tmp1 = sqrt ( v2 [0] * v2 [0] + v2 [1] * v2 [1] + v2 [2] * v2 [2] );

	kosinus = ( m [0] * n [0] + m [1] * n [1] + m [2] * n [2] ) * mn_norm;

	sinus =  ( v1 [0] * n [0] + v1 [1] * n [1] + v1 [2] * n [2] ) * tmp1 * mn_norm;

	psi = atan2 ( sinus, kosinus );

	/////////////////////////////////////////////////////////////////

	m_norm = kosinus * m_norm;

	n_norm = kosinus * n_norm;

	sinus = - force ( psi ) / sinus;

	/////////////////////////////////////////////////////////////////

	for (k = 0; k <= 2; k++)
	{
		A11 = 0; A22 = 0; A33 = 0; A12 = 0; A23 = 0; A13 = 0;

		for (l = 0; l <= 2; l++)
		{
			tmp1 = 1 - Kdelta (l, k);

			A11 += tmp1 * v1 [l] * v1 [l];

			A22 += tmp1 * v2 [l] * v2 [l];

			A33 += tmp1 * v3 [l] * v3 [l];

			A12 += tmp1 * v1 [l] * v2 [l];

			A23 += tmp1 * v2 [l] * v3 [l];

			A13 += tmp1 * v1 [l] * v3 [l];
		}

		/////////////////////////////////////////////////////////////////

		for (l = 1; l <= 4; l++)
		{
			tmp1 = Kdelta (l, 1);

			tmp2 = Kdelta (l, 2);

			tmp3 = Kdelta (l, 3);

			tmp4 = Kdelta (l, 4);

			/////////////////////////////////////////////////////////////////

			C1 = A22 * (tmp3 - tmp4) + A23 * (tmp3 - tmp2);

			C2 = A12 * (tmp3 - tmp2) + A22 * (tmp1 - tmp2);

			D1 = A12 * (tmp4 - tmp3) + A23 * (tmp2 - tmp1) + 2 * A13 * (tmp2 - tmp3);

			D2 = A11 * (tmp2 - tmp3) + A12 * (tmp2 - tmp1);

			D3 = A33 * (tmp2 - tmp3) + A23 * (tmp4 - tmp3);

			/////////////////////////////////////////////////////////////////

			tmp1 = C1 * v1 [k] + D1 * v2 [k] + C2 * v3 [k];

			tmp2 = C2 * v1 [k] + D2 * v2 [k];

			tmp3 = C1 * v3 [k] + D3 * v2 [k];

			/////////////////////////////////////////////////////////////////

			C1 = tmp1 * mn_norm + tmp2 * m_norm + tmp3 * n_norm;

			/////////////////////////////////////////////////////////////////

			tmp1 = sinus * C1;

			/////////////////////////////////////////////////////////////////

			if (k == 0) fx [ I [l] ] += tmp1;

			if (k == 1) fy [ I [l] ] += tmp1;

			if (k == 2) fz [ I [l] ] += tmp1;
		}
	}

	/////////////////////////////////////////////////////////////////

}

void polymer_constraints ( void )
{
	int i, k02, k00, k01, k10, k20, k11;

	/////////////////////////////////////////////////////////////////

	for (i = 2; i < N_amino + 1; i++)
	{
		k02 = atom_key [i - 1][2];  //Ca

		k00 = atom_key [i - 1][1];  //C

		k01 = atom_key [i - 1][0];  //O

		k10 = atom_key [i][3];      //N

		k20 = atom_key [i][2];      //Ca

		/////////////////////////////////////////////////////////////////

		if ( amino_key [i] == 'P' )   //proline
		{
			k11 = atom_key [i][6];

			/////////////////////////////////////////////////////////////////

			bond_force ( k00, k01, C_O_force );

			bond_force ( k00, k10, C_N_force );

			/////////////////////////////////////////////////////////////////

			valence_force ( k02, k00, k01, CT_C_O_force );

			valence_force ( k02, k00, k10, CT_C_N_force );

			valence_force ( k01, k00, k10, N_C_O_force );

			/////////////////////////////////////////////////////////////////

			valence_force ( k00, k10, k20, C_N_CT_force );

			valence_force ( k00, k10, k11, C_N_CT_force );

			/////////////////////////////////////////////////////////////////

			dihedral_force ( k02, k00, k10, k11, X2_C_N_X_force );

			dihedral_force ( k02, k00, k10, k20, X2_C_N_X_force );
		}

		/////////////////////////////////////////////////////////////////

		else
		{
			bond_force ( k00, k01, C_O_force );

			bond_force ( k00, k10, C_N_force );

			/////////////////////////////////////////////////////////////////

			valence_force ( k02, k00, k01, CT_C_O_force );

			valence_force ( k02, k00, k10, CT_C_N_force );

			valence_force ( k01, k00, k10, N_C_O_force );

			/////////////////////////////////////////////////////////////////

			valence_force ( k00, k10, k20, C_N_CT_force );

			/////////////////////////////////////////////////////////////////

			dihedral_force ( k02, k00, k10, k20, X2_N_C_X2_force );
		}
	}

	/////////////////////////////////////////////////////////////////

	k02 = atom_key [ N_amino ][2];

	k00 = atom_key [ N_amino ][1];

	k01 = atom_key [ N_amino ][0];

	/////////////////////////////////////////////////////////////////

	bond_force ( k00, k01, C_O_force );

	valence_force ( k02, k00, k01, CT_C_O_force );
}

void LJ_interactions_native_lists ( char * file_contacts )
{
	int i, j, k, i1, i2, j1;

	double tmp, r2, xx [3];

	//double Q;

	///////////////////////////////////////////////////////

	//Q = remaining_native_contacts ();

	//Q = k_NATIVE * ( Q0_NATIVE - Q ) / N0_NATIVE;

	///////////////////////////////////////////////////////

	for (k = 1; k < native_list_mass + 1; k += 2)
	{
		i = native_list_content [k];

		j = native_list_content [k + 1];

		///////////////////////////////////////////////////////

		xx [0] = x [i] - x [j];

		xx [1] = y [i] - y [j];

		xx [2] = z [i] - z [j];

		/////////////////////////////////////////////////////////////////

		r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

		/////////////////////////////////////////////////////////////////


    j1 = (k + 1) / 2;

		tmp = LJ_SCALE_n * F_LJ_n [j1] * LJ_FORCE_n ( r2, j1 );

		xx [0] = xx [0] * tmp;

		xx [1] = xx [1] * tmp;

		xx [2] = xx [2] * tmp;

		/////////////////////////////////////////////////////////////////

		fx [i] += xx [0]; fx [j] -= xx [0];

		fy [i] += xx [1]; fy [j] -= xx [1];

		fz [i] += xx [2]; fz [j] -= xx [2];
	}
}

void deterministic_forces ( char * file_contacts )
{
	int i;

	////////////////////////////////////////////

	for (i = 1; i < NT_atom + 1; i++) { fx [i] = 0; fy [i] = 0; fz [i] = 0; }

	////////////////////////////////////////////

	LJ_interactions_lists ();

	polymer_constraints ();

	//crowder_constraints ();

	LJ_interactions_native_lists ( file_contacts );

  // MTM_force();add by jdq

	//rg_harmonic_constraint ( 1, NTP_atom );
}

void full_forces ( void )
{
	int i, In1;

	double * maxwell_force;

	maxwell_force = (double *) calloc (3 * NT_atom + 1, sizeof(double));

	maxwell_1 ( maxwell_force, 3 * NT_atom );

	for (i = 1; i < NT_atom + 1; i++)
	{
		In1 = part_key [i];

		fx [i] = fx [i] + maxwell_force [i] * SIGMA_FORCE [In1];

		fy [i] = fy [i] + maxwell_force [i + NT_atom] * SIGMA_FORCE [In1];

		fz [i] = fz [i] + maxwell_force [i + 2 * NT_atom] * SIGMA_FORCE [In1];
	}

	free ( maxwell_force );
}
