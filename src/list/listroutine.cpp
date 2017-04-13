#include "listroutine.h"

int set_dr1()
{
  dr1 [0][0] = 1.5 * 0.75; //N-N

  dr1 [1][0] = 1.5 * 0.75; //N-C

  dr1 [2][0] = 1.5 * 0.75; //N-O

  dr1 [3][0] = 1.5 * 1.00; //N-H

  dr1 [4][0] = 1.5 * 0.75; //N-S

  dr1 [5][0] = 1.5 * 0.50; //N-a260                   what's mean of a260 ?

  dr1 [6][0] = 1.5 * 0.75; //C-C

  dr1 [7][0] = 1.5 * 0.75; //C-O

  dr1 [8][0] = 1.5 * 1.00; //C-H

  dr1 [9][0] = 1.5 * 0.75; //C-S

  dr1 [10][0] = 1.5 * 0.50; //C-a260

  dr1 [11][0] = 1.5 * 0.75; //O-O

  dr1 [12][0] = 1.5 * 1.00; //O-H

  dr1 [13][0] = 1.5 * 0.75; //O-S

  dr1 [14][0] = 1.5 * 0.50; //O-a260

  dr1 [15][0] = 1.5 * 1.25; //H-H

  dr1 [16][0] = 1.5 * 1.00; //H-S

  dr1 [17][0] = 1.5 * 1.00; //H-a260

  dr1 [18][0] = 1.5 * 0.75; //S-S

  dr1 [19][0] = 1.5 * 0.50; //S-a260

  dr1 [20][0] = 1.5 * 0.50; //a260-a260

  return 1;
}

int list_crowder_types ()
{
	int i, j;

	/////////////////////////////////////////////////

	crowder_mass = (int *) calloc ( n_crwd, sizeof(int) );

	crowder_content = (int **) calloc ( n_crwd, sizeof(int *) );

	/////////////////////////////////////////////////////////////////

	for (i = 0; i < n_crwd; i++)
	{
		switch (i)
		{

		case 0://water

			m_crwd [i] = 3;

			V_crwd [i] = 4.0 * pi * pow ( RLJ [258], 3 ) / 3.0;

			/////////////////////////////////////////////////
			//rigid water

			RIGID_SET_crwd [i] = 0;

			RIGID_END_crwd [i] = m_crwd [i] - 1;

			/////////////////////////////////////////////////
			//flexible water

			//RIGID_SET_crwd [i] = 0;

			//RIGID_END_crwd [i] = -1;

			/////////////////////////////////////////////////

			part_key_crwd [i] = NULL;

			part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			part_key_crwd [i][0] = 2;

			part_key_crwd [i][1] = 3;

			part_key_crwd [i][2] = 3;

			/////////////////////////////////////////////////

			maxi_key_crwd [i] = NULL;

			maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			maxi_key_crwd [i][0] = 258;

			maxi_key_crwd [i][1] = 259;

			maxi_key_crwd [i][2] = 259;

			/////////////////////////////////////////////////

			RX_crwd [i] = NULL;

			RY_crwd [i] = NULL;

			RZ_crwd [i] = NULL;

			RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RX_crwd [i][0] = 0;

			RY_crwd [i][0] = 0;

			RZ_crwd [i][0] = 0;

			RX_crwd [i][1] = r_OW_HW * cos ( 0.5 * (pi - theta_HW_OW_HW) );

			RY_crwd [i][1] = r_OW_HW * sin ( 0.5 * (pi - theta_HW_OW_HW) );

			RZ_crwd [i][1] = 0;

			RX_crwd [i][2] = - r_OW_HW * cos ( 0.5 * (pi - theta_HW_OW_HW) );

			RY_crwd [i][2] = r_OW_HW * sin ( 0.5 * (pi - theta_HW_OW_HW) );

			RZ_crwd [i][2] = 0;

			/////////////////////////////////////////////////

			break;

			/////////////////////////////////////////////////

		case 1://sphere 260

			m_crwd [i] = 1;

			V_crwd [i] = m_crwd [i] * 4.0 * pi * pow ( RLJ [260], 3 ) / 3.0;

			/////////////////////////////////////////////////

			RIGID_SET_crwd [i] = 0;

			RIGID_END_crwd [i] = -1;

			/////////////////////////////////////////////////

			part_key_crwd [i] = NULL;

			part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			part_key_crwd [i][0] = 5;

			/////////////////////////////////////////////////

			maxi_key_crwd [i] = NULL;

			maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			maxi_key_crwd [i][0] = 260;

			/////////////////////////////////////////////////

			RX_crwd [i] = NULL;

			RY_crwd [i] = NULL;

			RZ_crwd [i] = NULL;

			RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RX_crwd [i][0] = 0;

			RY_crwd [i][0] = 0;

			RZ_crwd [i][0] = 0;

			/////////////////////////////////////////////////

			break;

			/////////////////////////////////////////////////

		case 2://rigid cylinder 260

			m_crwd [i] = 4;

			V_crwd [i] = m_crwd [i] * 4.0 * pi * pow ( RLJ [260], 3 ) / 3.0;

			/////////////////////////////////////////////////

			RIGID_SET_crwd [i] = 0;

			RIGID_END_crwd [i] = m_crwd [i] - 1;

			/////////////////////////////////////////////////

			part_key_crwd [i] = NULL;

			part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			for (j = 0; j < m_crwd [i]; j++)

				part_key_crwd [i][j] = 5;

			/////////////////////////////////////////////////

			maxi_key_crwd [i] = NULL;

			maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			for (j = 0; j < m_crwd [i]; j++)

				maxi_key_crwd [i][j] = 260;

			/////////////////////////////////////////////////

			RX_crwd [i] = NULL;

			RY_crwd [i] = NULL;

			RZ_crwd [i] = NULL;

			RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			for (j = 0; j < m_crwd [i]; j++)
			{
				RX_crwd [i][j] = j * r_a260_a260;

				RY_crwd [i][j] = 0;

				RZ_crwd [i][j] = 0;
			}

			/////////////////////////////////////////////////

			break;

			/////////////////////////////////////////////////

		case 3://peptide 260

			m_crwd [i] = 4;

			V_crwd [i] = m_crwd [i] * 4.0 * pi * pow ( RLJ [260], 3 ) / 3.0;

			/////////////////////////////////////////////////
			//rigid peptide

			//RIGID_SET_crwd [i] = 0;

			//RIGID_END_crwd [i] = m_crwd [i] - 1;

			/////////////////////////////////////////////////
			//flexible peptide

			RIGID_SET_crwd [i] = 0;

			RIGID_END_crwd [i] = -1;

			/////////////////////////////////////////////////

			part_key_crwd [i] = NULL;

			part_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			for (j = 0; j < m_crwd [i]; j++)

				part_key_crwd [i][j] = 5;

			/////////////////////////////////////////////////

			maxi_key_crwd [i] = NULL;

			maxi_key_crwd [i] = (int *) calloc ( m_crwd [i], sizeof(int) );

			for (j = 0; j < m_crwd [i]; j++)

				maxi_key_crwd [i][j] = 260;

			/////////////////////////////////////////////////

			RX_crwd [i] = NULL;

			RY_crwd [i] = NULL;

			RZ_crwd [i] = NULL;

			RX_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RY_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			RZ_crwd [i] = (double *) calloc ( m_crwd [i], sizeof(double) );

			for (j = 0; j < m_crwd [i]; j++)
			{
				RX_crwd [i][j] = j * r_a260_a260 * cos ( 0.5 * (pi - theta_a260_a260_a260) );

				RY_crwd [i][j] = 0.5 * ( 1 - pow (-1, j) ) * r_a260_a260 * sin ( 0.5 * (pi - theta_a260_a260_a260) );

				RZ_crwd [i][j] = 0;
			}

			/////////////////////////////////////////////////

			break;

			/////////////////////////////////////////////////

		default:

			break;
		}
	}
  return 1;
}

void add_bond ( int i1, int i2 )
{
	bond_list_mass += 2;

	bond_list_content = (int *) realloc ( bond_list_content, ( bond_list_mass + 1 ) * sizeof(int) );

	bond_list_content [ bond_list_mass - 1 ] = i1;

	bond_list_content [ bond_list_mass ] = i2;

}

void add_valence ( int i1, int i2 )
{
	valence_list_mass += 2;

	valence_list_content = (int *) realloc ( valence_list_content, ( valence_list_mass + 1 ) * sizeof(int) );

	valence_list_content [ valence_list_mass - 1 ] = i1;

	valence_list_content [ valence_list_mass ] = i2;
}

void add_dihedral ( int i1, int i2 )
{
	dihedral_list_mass += 2;

	dihedral_list_content = (int *) realloc ( dihedral_list_content, ( dihedral_list_mass + 1 ) * sizeof(int) );

	dihedral_list_content [ dihedral_list_mass - 1 ] = i1;

	dihedral_list_content [ dihedral_list_mass ] = i2;
}

int populate_bvd_lists ()
{
	int i, k02, k00, k01, k10, k20, k11;

	int k22, k30, k32;

	/////////////////////////////////////////////////////////////////

	for (i = 2; i < N_amino + 1; i++) ////marked by QM.
	{
		k02 = atom_key [i - 1][2];   // Ca

		k00 = atom_key [i - 1][1];   // C

		k01 = atom_key [i - 1][0];   // O

		k10 = atom_key [i][3];       // N

		k20 = atom_key [i][2];       // Ca

		/////////////////////////////////////////////////////////////////

		if ( amino_key [i] == 'P' ) //proline
		{
			k11 = atom_key [i][6];  // Atom C_d

			/////////////////////////////////////////////////////////////////

			add_bond ( k00, k01 );  // bond C=O

			add_bond ( k00, k10 );  // peptide bond C-N

			/////////////////////////////////////////////////////////////////

			add_valence ( k02, k01 );

			add_valence ( k02, k10 );

			add_valence ( k01, k10 );

			/////////////////////////////////////////////////////////////////

			add_valence ( k00, k20 );

			add_valence ( k00, k11 );   //// This is especially for Proline

			/////////////////////////////////////////////////////////////////

			add_dihedral ( k02, k11 );   /// This is especially for Proline.

			add_dihedral ( k02, k20 );

			add_dihedral ( k01, k11 );   /// This is especially for Proline.

			add_dihedral ( k01, k20 );
		}

		/////////////////////////////////////////////////////////////////

		else
		{
			add_bond ( k00, k01 );

			add_bond ( k00, k10 );

			/////////////////////////////////////////////////////////////////

			add_valence ( k02, k01 );

			add_valence ( k02, k10 );

			add_valence ( k01, k10 );

			/////////////////////////////////////////////////////////////////

			add_valence ( k00, k20 );

			/////////////////////////////////////////////////////////////////

			add_dihedral ( k02, k20 );   // Ca-C-N-Ca

			add_dihedral ( k01, k20 );   // O=C-N-Ca
		}
	}

	/////////////////////////////////////////////////////////////////

	k02 = atom_key [ N_amino ][2];  // Ca for last amino acid in protein sequence

	k00 = atom_key [ N_amino ][1];  // O

	k01 = atom_key [ N_amino ][0];  // C

	/////////////////////////////////////////////////////////////////

	add_bond ( k00, k01 );        // C=O for last amino acid

	add_valence ( k02, k01 );     // Ca-C=O

	/////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////

	for (i = 1; i < N_amino + 1; i++)
	{
		if (i > 1) k00 = atom_key [i - 1][1];

		k10 = atom_key [i][3];

		k22 = atom_key [i][4];

		k30 = atom_key [i][1];

		if (i < N_amino) k32 = atom_key [i + 1][3];

		/////////////////////////////////////////////////////////////////

		if ( amino_key [i] == 'G' ) //glycine
		{
			if (i > 1) add_dihedral ( k00, k30 );

			if (i < N_amino) add_dihedral ( k10, k32 );
		}

		/////////////////////////////////////////////////////////////////

		else if ( amino_key [i] == 'P' ) //proline
		{
			if (i > 1)
			{
				add_dihedral ( k00, k30 );

				add_dihedral ( k00, k22 );

				add_dihedral ( k00, atom_key [i][5] );
			}

			if (i < N_amino)
			{
				add_dihedral ( k10, k32 );

				add_dihedral ( k22, k32 );
			}
		}

		/////////////////////////////////////////////////////////////////

		else
		{
			if (i > 1)
			{
				add_dihedral ( k00, k30 );

				add_dihedral ( k00, k22 );
			}

			if (i < N_amino)
			{
				add_dihedral ( k10, k32 );

				add_dihedral ( k22, k32 );
			}
		}
	}
  return 1;
}

void add_to_native_list ( int i, int j, double r2, double E )
{
	int k;

	/////////////////////////////////////////////////////////////////

	native_list_mass += 2;

	native_list_content = (int *) realloc ( native_list_content, ( native_list_mass + 1 ) * sizeof(int) );

	native_list_content [ native_list_mass - 1 ] = i;

	native_list_content [ native_list_mass ] = j;

	/////////////////////////////////////////////////////////////////

	k = native_list_mass / 2;

	D2_LJ_n = (double *) realloc ( D2_LJ_n, (k + 1) * sizeof(double) );

	D2_LJ_CUTOFF_n = (double *) realloc ( D2_LJ_CUTOFF_n, (k + 1) * sizeof(double) );

	E_LJ_n = (double *) realloc ( E_LJ_n, (k + 1) * sizeof(double) );

	F_LJ_n = (double *) realloc ( F_LJ_n, (k + 1) * sizeof(double) );

	/////////////////////////////////////////////////////////////////

	D2_LJ_n [k] = r2;

	D2_LJ_CUTOFF_n [k] = 9.0 * D2_LJ_n [k];

	E_LJ_n [k] = E;

	F_LJ_n [k] = 12.0 * E_LJ_n [k] / D2_LJ_n [k];
}

int populate_native_lists ( char * file_native_list, double cutoff )
{
	FILE * f1;

	int i, j, i1, i2, j1, k1;

	double xx [3], r2;

	// double r21, r22;

	/////////////////////////////////////////////////////////////////

	for (i = 1; i < NTP_atom + 1; i++)
	{
		for (j = i + 1; j < NTP_atom + 1; j++)
		{
			//if ( INDX [i] == INDX [j] && JNDX [i] > 0 && JNDX [j] > 0 ) continue;

			if ( INDX [i] == INDX [j] ) continue;

			///////////////////////////////////////////////////////

			k1 = 0;

			///////////////////////////////////////////////////////
			for (j1 = 1; j1 < bond_list_mass + 1; j1 += 2)       ///// exclude bond case between atom i and i+1
			{
				i1 = bond_list_content [j1];

				i2 = bond_list_content [j1 + 1];

				if ( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) ) { k1 = 1; break; }
			}

			if ( k1 ) continue;

			///////////////////////////////////////////////////////
			for (j1 = 1; j1 < valence_list_mass + 1; j1 += 2) ///// exclude valence case between atom i and i+2
			{
				i1 = valence_list_content [j1];

				i2 = valence_list_content [j1 + 1];

				if ( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) ) { k1 = 1; break; }
			}

			if ( k1 ) continue;

			///////////////////////////////////////////////////////

			for (j1 = 1; j1 < dihedral_list_mass + 1; j1 += 2)   ///// exclude dihedral case between atom i and i+3
			{
				i1 = dihedral_list_content [j1];

				i2 = dihedral_list_content [j1 + 1];

				if ( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) ) { k1 = 1; break; }
			}

			if ( k1 ) continue;

			///////////////////////////////////////////////////////

			i1 = maxi_key [i];

			i2 = maxi_key [j];

			if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

			else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

			/////////////////////////////////////////////////////////////////

			xx [0] = x [i] - x [j];

			xx [1] = y [i] - y [j];

			xx [2] = z [i] - z [j];

			r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];  // squre of the distance between two atoms in each pair

			if ( r2 < cutoff * cutoff * D2_LJ [j1] ) add_to_native_list ( i, j, r2, E_LJ [j1] );

			/////////////////////////////////////////////////////////////////

		}
	}

	/////////////////////////////////////////////////////////////////

	f1 = fopen ( file_native_list, "w" );

	fprintf ( f1, "%d\n", native_list_mass );

	for (i = 1; i < native_list_mass + 1; i += 2)

		fprintf ( f1, "%d %d\n", native_list_content [i], native_list_content [i + 1] );

	i1 = native_list_mass / 2;

	for (i = 1; i < i1 + 1; i++)

		fprintf ( f1, "%le %le %le %le\n", D2_LJ_n [i], D2_LJ_CUTOFF_n [i], E_LJ_n [i], F_LJ_n [i] );

	fclose (f1);

	/////////////////////////////////////////////////////////////////

	N0_NATIVE = i1;

	THETA_NATIVE = (double *) calloc ( i1 + 1, sizeof(double) );

  return 1;
}

int separate_in_cells ( int In1, int In2, int InInt )
{
	FILE * f1;

	int i, m, In12;

	int i1, i2, i3;

	double tmp, lX, lY, lZ;

	/////////////////////////////////////////////////////////////

	if ( McT )
	{
		for (m = 1; m < McT + 1; m++) free ( cell_content [m] );

		free ( cell_content ); free ( cell_mass );

		cell_content = NULL; cell_mass = NULL;

		McX = 0; McY = 0; McZ = 0; McT = 0;
	}

	/////////////////////////////////////////////////////////////

	if ( InInt == 0 )
	{
		In12 = In1 * ns - In1 * (1 + In1) / 2 + In2;

		tmp = 1.0 * PD [ In12 ]; //excluded volume only
	}

	R2_CUTOFF [ In12 ][ InInt ] = tmp * tmp;

	R1_LIST [ In12 ][ InInt ] = tmp + dr1 [ In12 ][ InInt ];

	R2_LIST [ In12 ][ InInt ] = R1_LIST [ In12 ][ InInt ] * R1_LIST [ In12 ][ InInt ];

	/////////////////////////////////////////////////////////////

	tmp = R1_LIST [ In12 ][ InInt ];

	if ( tmp > side_X || tmp > side_Y || tmp > side_Z )
	{
		f1 = fopen ( "error.dat", "a" );

		fprintf ( f1, "Error in <separate_in_cells> for list [%d][%d]: cell size = %le, side_X = %le, side_Y = %le, side_Z = %le\n", In12, InInt, tmp, side_X, side_Y, side_Z );

		fclose ( f1 );

		abort ();
	}

	/////////////////////////////////////////////////////////////

	McX = (int) ( side_X / tmp ); lX = side_X / McX;

	McY = (int) ( side_Y / tmp ); lY = side_Y / McY;

	McZ = (int) ( side_Z / tmp ); lZ = side_Z / McZ;

	McT = McX * McY * McZ;

	/////////////////////////////////////////////////////////////

	cell_content = (int **) calloc ( McT + 1, sizeof(int *) );

	for (m = 0; m < McT + 1; m++) cell_content [m] = NULL;

	cell_mass = (int *) calloc ( McT + 1, sizeof(int) );

	/////////////////////////////////////////////////////////////

	for (i = 1; i < NT_atom + 1; i++)
	{
		In12 = part_key [i];

		if ( In12 == In1 || In12 == In2 )
		{
			tmp = x [i];

			if ( tmp < 0 ) tmp += side_X;

			else if ( tmp > side_X ) tmp -= side_X;

			i1 = (int) ceil (tmp / lX);

			/////////////////////////////////////////////////////////////

			tmp = y [i];

			if ( tmp < 0 ) tmp += side_Y;

			else if ( tmp > side_Y ) tmp -= side_Y;

			i2 = (int) ceil (tmp / lY);

			/////////////////////////////////////////////////////////////

			tmp = z [i];

			if ( tmp < 0 ) tmp += side_Z;

			else if ( tmp > side_Z ) tmp -= side_Z;

			i3 = (int) ceil (tmp / lZ);

			/////////////////////////////////////////////////////////////

			m = (i1 - 1) * McY * McZ + (i2 - 1) * McZ + i3;

			cell_mass [m] += 1;

			cell_content [m] = (int *) realloc ( cell_content [m], (cell_mass [m] + 1) * sizeof(int) );

			cell_content [m] [ cell_mass [m] ] = i;
		}
	}
  return 1;
}

int populate_list ( int In1, int In2, int InInt ) ///QM: In1=0;<ns and In2=In1;<ns
{
	int J1, J2, k;

	int i1, i2;

	int j1, j2;

	int k1, k2;

	int l1, l2;

	int m1, m2, m3;

	int n1, n2, n3;

	int In12, tmp;

	double r2, xx [3];

	int index [13][3] = { { -1, -1, -1 }, { -1, -1, 0 },  { -1, -1, 1 },  { -1, 0, -1 }, { -1, 0, 0 }, { -1, 0, 1 },

	{ -1, 1, -1 }, { -1, 1, 0 }, { -1, 1, 1 }, { 0, -1, -1 }, { 0, -1, 0 }, { 0, -1, 1 }, { 0, 0, -1 } };

	//////////////////////////////////////////

	In12 = In1 * ns - In1 * (1 + In1) / 2 + In2; ////////QM: In12 [0, 20]

	if ( list_mass [ In12 ][ InInt ] )
	{
		free ( list_content [ In12 ][ InInt ] );

		list_content [ In12 ][ InInt ] = NULL;

		list_mass [ In12 ][ InInt ] = 0;
	}

	//////////////////////////////////////////

	separate_in_cells ( In1, In2, InInt );

	//////////////////////////////////////////

	for (m1 = 1; m1 < McX + 1; m1++)
	{
		for (m2 = 1; m2 < McY + 1; m2++)
		{
			for (m3 = 1; m3 < McZ + 1; m3++)
			{
				J1 = (m1 - 1) * McY * McZ + (m2 - 1) * McZ + m3;

				for (k = 0; k < 13; k++)
				{
					n1 = m1 + index [k][0];

					if (n1 > McX) { n1 -= McX; if ( n1 == m1 || n1 == m1 - 1 ) continue; }

					else if (n1 < 1) { n1 += McX; if ( n1 == m1 || n1 == m1 + 1 ) continue; }

					n2 = m2 + index [k][1];

					if (n2 > McY) { n2 -= McY; if ( n2 == m2 || n2 == m2 - 1 ) continue; }

					else if (n2 < 1) { n2 += McY; if ( n2 == m2 || n2 == m2 + 1 ) continue; }

					n3 = m3 + index [k][2];

					if (n3 > McZ) { n3 -= McZ; if ( n3 == m3 || n3 == m3 - 1 ) continue; }

					else if (n3 < 1) { n3 += McZ; if ( n3 == m3 || n3 == m3 + 1 ) continue; }

					J2 = (n1 - 1) * McY * McZ + (n2 - 1) * McZ + n3;

					for (j1 = 1; j1 < cell_mass [J1] + 1; j1++)
					{
						for (j2 = 1; j2 < cell_mass [J2] + 1; j2++)
						{
							i1 = cell_content [J1][j1]; i2 = cell_content [J2][j2];

							if (In1 < In2)
							{
								n1 = part_key [i1]; n2 = part_key [i2];

								if (n1 == n2) continue;

								if (n2 < n1) { tmp = i1; i1 = i2; i2 = tmp; }
							}

							///////////////////////////////////////////////////////

							n1 = INDX [i1];

							n2 = INDX [i2];

							//////////////////////////////////////////// How to consider Rigid Part? July 17, 2014

							if (n1 == n2)
							{
								l1 = JNDX [i1];

								l2 = JNDX [i2];

								k1 = RIGID_SET [n1];

								k2 = RIGID_END [n1];

								if (l1 >= k1 && l1 <= k2 && l2 >= k1 && l2 <= k2) continue;  // two heavy atoms in ONE rigid body.
							}

							///////////////////////////////////////////////////////

							CC_vector ( i1, i2, xx );

							if ( fabs (xx [0]) > R1_LIST [ In12 ][ InInt ] ) continue;

							if ( fabs (xx [1]) > R1_LIST [ In12 ][ InInt ] ) continue;

							if ( fabs (xx [2]) > R1_LIST [ In12 ][ InInt ] ) continue;

							r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

							if ( r2 <= R2_LIST [ In12 ][ InInt ] )
							{
								list_mass [ In12 ][ InInt ] += 2; tmp = list_mass [ In12 ][ InInt ];

								list_content [ In12 ][ InInt ] = (int *) realloc ( list_content [ In12 ][ InInt ], (tmp + 1) * sizeof(int) );

								list_content [ In12 ][ InInt ][ tmp - 1 ] = i1; list_content [ In12 ][ InInt ][ tmp ] = i2;
							}
						}
					}
				}

				/////////////////////////////////////////

				for (j1 = 1; j1 < cell_mass [J1] + 1; j1++)
				{
					for (j2 = j1 + 1; j2 < cell_mass [J1] + 1; j2++)
					{
						i1 = cell_content [J1][j1]; i2 = cell_content [J1][j2];

						if (In1 < In2)
						{
							n1 = part_key [i1]; n2 = part_key [i2];

							if (n1 == n2) continue;

							if (n2 < n1) { tmp = i1; i1 = i2; i2 = tmp; }
						}

						///////////////////////////////////////////////////////

						n1 = INDX [i1];

						n2 = INDX [i2];

						///////////////////////////////////////////////////////

						if (n1 == n2)
						{
							l1 = JNDX [i1];

							l2 = JNDX [i2];

							k1 = RIGID_SET [n1];

							k2 = RIGID_END [n1];

							if (l1 >= k1 && l1 <= k2 && l2 >= k1 && l2 <= k2) continue;
						}

						///////////////////////////////////////////////////////

						CC_vector ( i1, i2, xx );

						r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

						if ( r2 <= R2_LIST [ In12 ][ InInt ] )
						{
							list_mass [ In12 ][ InInt ] += 2; tmp = list_mass [ In12 ][ InInt ];

							list_content [ In12 ][ InInt ] = (int *) realloc ( list_content [ In12 ][ InInt ], (tmp + 1) * sizeof(int) );

							list_content [ In12 ][ InInt ][ tmp - 1 ] = i1; list_content [ In12 ][ InInt ][ tmp ] = i2;
						}
					}
				}
			}
		}
	}

  return 1;
}

void check_shifts ( void )
{
	int i;

	int In1, In2, In12, InInt;

	double x_tmp, y_tmp, z_tmp, r2;

	double shift [ ns ];

	/////////////////////////////////////

	for (In1 = 0; In1 < ns; In1++) shift [ In1 ] = 0;

	for (i = 1; i < NT_atom + 1; i++)
	{
		In1 = part_key [i];

		x_tmp = x [i] - x_old [i];

		y_tmp = y [i] - y_old [i];

		z_tmp = z [i] - z_old [i];

		r2 = x_tmp * x_tmp + y_tmp * y_tmp + z_tmp * z_tmp;

		if ( r2 > shift [ In1 ] ) shift [ In1 ] = r2;
	}

	for (In1 = 0; In1 < ns; In1++) shift [ In1 ] = sqrt ( shift [ In1 ] );

	/////////////////////////////////////

	//for (In1 = 0; In1 < ns; In1++) printf ("shift [%d] = %le ", In1, shift [In1]);

	//printf ("\n\n");

	//getchar ();

	/////////////////////////////////////

	for (In1 = 0; In1 < ns; In1++)
	{
		for (In2 = In1; In2 < ns; In2++)
		{
			In12 = In1 * ns - In1 * (1 + In1) / 2 + In2;

			for (InInt = 0; InInt < npi; InInt++)
			{
				PROGRESS [ In12 ][ InInt ] += ( shift [ In1 ] + shift [ In2 ] );

				if ( PROGRESS [ In12 ][ InInt ] > dr1 [ In12 ][ InInt ] )
				{
					//printf ("%ld, PROGRESS [ %d ][ %d ] = %le: populate list?\n", klok, In12, InInt, PROGRESS [ In12 ][ InInt ] );

					//getchar ();

					populate_list ( In1, In2, InInt );

					PROGRESS [ In12 ][ InInt ] = 0;

					//printf ("list_mass [%d][%d] = %d\n\n", In12, InInt, list_mass [In12][InInt] );
				}
			}
		}
	}
}

int flag_any_2_atoms()
{
    //'b'->bonds;'a'->angles;'d'->dihedral;'c'->native contacts;
    //'s'->same;'o'->others
    int i,j,i1,i2,j1,k1;

    Flags2Atom = (char **)malloc((NTP_atom+1)*sizeof(char *));
    for(i=1;i<NTP_atom+1;i++)
    {
        Flags2Atom[i]=(char *)malloc((NTP_atom+1)*sizeof(char));
    }

    for(i=1;i<NTP_atom+1;i++)
    {
        for(j=1;j<NTP_atom+1;j++)
        {

            if (i==j)
            {
                Flags2Atom[i][j]='s';
                continue;
            }

            k1=0;

            for(j1=1;j1<bond_list_mass+1;j1+=2)
            {
                i1=bond_list_content[j1];
                i2=bond_list_content[j1+1];
                if( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) )
                {
                    Flags2Atom[i][j]='b';
                    k1=1;
                    break;
                }
            }

            if(k1) continue;

            for (j1 = 1; j1 < valence_list_mass + 1; j1 += 2) ///// exclude valence case between atom i and i+2
			{
				i1 = valence_list_content [j1];

				i2 = valence_list_content [j1 + 1];

				if ( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) )
                {
                    Flags2Atom[i][j]='a';
                    k1=1;
                    break;
                }
			}

			if ( k1 ) continue;

			for (j1 = 1; j1 < dihedral_list_mass + 1; j1 += 2)   ///// exclude dihedral case between atom i and i+3
			{
				i1 = dihedral_list_content [j1];

				i2 = dihedral_list_content [j1 + 1];

				if ( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) )
                {
                    Flags2Atom[i][j]='d';
                    k1=1;
                    break;
                }
			}

			if ( k1 ) continue;

			for (j1 = 1; j1 < native_list_mass + 1; j1 += 2)
            {
                i1 = native_list_content [j1];

                i2 = native_list_content [j1 + 1];

                if ( ( i == i1 && j == i2 ) || ( i == i2 && j == i1 ) )
                {
                    Flags2Atom[i][j]='c';
                    k1=1;
                    break;
                }
            }

            if ( k1 ) continue;

			Flags2Atom[i][j]='o';

        }
    }

    return 1;

}
