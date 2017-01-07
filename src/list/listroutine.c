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
