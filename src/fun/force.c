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
