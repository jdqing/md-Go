#include "energy.h"

double WCA_ENERGY ( double x2, int Jn12 )//( double x2, int In12, int Jn12 )
{
	if ( x2 > D2_LJ_CUTOFF [ Jn12 ] ) return (0);

	else
	{
		double r2, r4, r6;

		r2 = D2_LJ [ Jn12 ] / x2;

		r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;

		r4 = E_LJ [ Jn12 ] * ( r2 - 2.0 * r6 + 1.0 );

		return ( r4 );
	}
}

double LJ_ENERGY ( double x2, int Jn12 )//( double x2, int In12, int Jn12 )
{
	if ( x2 > D2_LJ_CUTOFF [ Jn12 ] ) return (0);

	else
	{
		double r2, r4, r6;

		r2 = D2_LJ [ Jn12 ] / x2;

		r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;

		r4 = E_LJ [ Jn12 ] * ( r2 - 2.0 * r6 );

		return ( r4 );
	}
}

double LJ_ENERGY_n ( double x2, int Jn12 )
{
	if ( x2 > D2_LJ_CUTOFF_n [ Jn12 ] ) return (0);

	else
	{
		double r2, r4, r6;

		r2 = D2_LJ_n [ Jn12 ] / x2;

		r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;

		r4 = E_LJ_n [ Jn12 ] * ( r2 - 2.0 * r6 );

		return ( r4 );
	}
}

double bond_energy ( int i1, int i2, double ( * energy) ( double r ) )
{
	double v1 [3], r1, tmp1;

	/////////////////////////////////////////////////////////////////

	v1 [0] = x [i2] - x [i1];

	v1 [1] = y [i2] - y [i1];

	v1 [2] = z [i2] - z [i1];

	/////////////////////////////////////////////////////////////////

	r1 = sqrt ( v1 [0] * v1 [0] + v1 [1] * v1 [1] + v1 [2] * v1 [2] );

	tmp1 = energy ( r1 );

	/////////////////////////////////////////////////////////////////

	return ( tmp1 );
}


double valence_energy ( int i1, int i2, int i3, double ( * energy) ( double phi ) )
{
	double v1 [3], v2 [3], v1_norm, v2_norm, v12_norm, kosinus, theta, tmp1;

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

	tmp1 = energy ( theta );

	/////////////////////////////////////////////////////////////////

	return ( tmp1 );
}



double dihedral_energy ( int i1, int i2, int i3, int i4, double ( * energy) ( double phi ) )
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

	sinus =  ( v1 [0] * n [0] + v1 [1] * n [1] + v1 [2] * n [2] ) * tmp1;

	psi = atan2 ( sinus, kosinus );

	tmp1 = energy ( psi );

	/////////////////////////////////////////////////////////////////

	return ( tmp1 );
}

double E_total_bvd()
{
	int i, k02, k00, k01, k10, k20, k11;

	double Etotal=0;

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

			Etotal = Etotal + bond_energy ( k00, k01, C_O_energy );

			Etotal = Etotal + bond_energy ( k00, k10, C_N_energy );

			/////////////////////////////////////////////////////////////////

			Etotal = Etotal + valence_energy ( k02, k00, k01, CT_C_O_energy );

			Etotal = Etotal + valence_energy ( k02, k00, k10, CT_C_N_energy );

			Etotal = Etotal + valence_energy ( k01, k00, k10, N_C_O_energy );

			/////////////////////////////////////////////////////////////////

			Etotal = Etotal + valence_energy ( k00, k10, k20, C_N_CT_energy );

			Etotal = Etotal + valence_energy ( k00, k10, k11, C_N_CT_energy );

			/////////////////////////////////////////////////////////////////

			Etotal = Etotal + dihedral_energy ( k02, k00, k10, k11, X2_C_N_X_energy );

			Etotal = Etotal + dihedral_energy ( k02, k00, k10, k20, X2_C_N_X_energy );
		}

		/////////////////////////////////////////////////////////////////

		else
		{
			Etotal = Etotal + bond_energy ( k00, k01, C_O_energy );

			Etotal = Etotal + bond_energy ( k00, k10, C_N_energy );

			/////////////////////////////////////////////////////////////////

			Etotal = Etotal + valence_energy ( k02, k00, k01, CT_C_O_energy );

			Etotal = Etotal + valence_energy ( k02, k00, k10, CT_C_N_energy );

			Etotal = Etotal + valence_energy ( k01, k00, k10, N_C_O_energy );

			/////////////////////////////////////////////////////////////////

			Etotal = Etotal + valence_energy ( k00, k10, k20, C_N_CT_energy );

			/////////////////////////////////////////////////////////////////

			Etotal = Etotal + dihedral_energy ( k02, k00, k10, k20, X2_N_C_X2_energy );
		}
	}

	/////////////////////////////////////////////////////////////////

	k02 = atom_key [ N_amino ][2];

	k00 = atom_key [ N_amino ][1];

	k01 = atom_key [ N_amino ][0];

	/////////////////////////////////////////////////////////////////

	Etotal = Etotal + bond_energy ( k00, k01, C_O_energy );

	Etotal = Etotal + valence_energy ( k02, k00, k01, CT_C_O_energy );

	return Etotal;

}
