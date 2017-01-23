#include "dealfor.h"

double radius_of_gyration ( int n1, int n2 )
{
	int i;

	double xx [3], r1 = 0;

	double xcm = 0, ycm = 0, zcm = 0;

	int N = 0;

	////////////////////////////////////////////////////

	for (i = n1; i < n2 + 1; i++) { N++; xcm += x [i]; ycm += y [i]; zcm += z [i]; }

	xcm /= N; ycm /= N; zcm /= N;

	////////////////////////////////////////////////////

	for (i = n1; i < n2 + 1; i++)
	{
		xx [0] = x [i] - xcm;

		xx [1] = y [i] - ycm;

		xx [2] = z [i] - zcm;

		r1 += xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];
	}

	r1 = sqrt ( r1 / N );

	return ( r1 );
}

double remaining_native_contacts ()
{
//	FILE * f1;

	int i, j, k, j1;

	double count = 0;

	double xx [3], r2;

	///////////////////////////////////////////////////////

	for (k = 1; k < native_list_mass + 1; k += 2)
	{
		j1 = (k + 1) / 2;

		///////////////////////////////////////////////////////

		i = native_list_content [k];

		j = native_list_content [k + 1];

		///////////////////////////////////////////////////////

		xx [0] = x [i] - x [j];

		xx [1] = y [i] - y [j];

		xx [2] = z [i] - z [j];

		r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

		///////////////////////////////////////////////////////

		if ( r2 < D2_LJ_n [j1] )	count += 1.0;

		else if ( r2 < D2_LJ_CUTOFF_n [j1] )
		{
			xx [0] = D2_LJ_n [j1] / r2;

			xx [1] = xx [0] * xx [0];

			xx [2] = xx [1] * xx [0];

			xx [0] = xx [2] * xx [2];

			count += ( 2.0 * xx [2] - xx [0] );
		}


		else count += 0;
	}

	///////////////////////////////////////////////////////

	return ( count / N0_NATIVE );
}
int deal_for( const char *flag)
{

}
