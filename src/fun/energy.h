#ifndef _ENERGY_H
#define _ENERGY_H
#include "../init/init.h"
#include "potential.h"

double WCA_ENERGY ( double x2, int Jn12 );

double LJ_ENERGY ( double x2, int Jn12 );//( double x2, int In12, int Jn12 )

double LJ_ENERGY_n ( double x2, int Jn12 );

double bond_energy ( int i1, int i2, double ( * energy) ( double r ) );

double valence_energy ( int i1, int i2, int i3, double ( * energy) ( double phi ) );

double dihedral_energy ( int i1, int i2, int i3, int i4, double ( * energy) ( double phi ) );

double E_total_bvd();

#endif
