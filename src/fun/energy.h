#ifndef _ENERGY_H
#define _ENERGY_H
#include "../init/init.h"

double WCA_ENERGY ( double x2, int Jn12 );

double LJ_ENERGY ( double x2, int Jn12 );//( double x2, int In12, int Jn12 )

double LJ_ENERGY_n ( double x2, int Jn12 );

double C_N_energy ( double r );
double C_O_energy ( double r );

double bond_energy ( int i1, int i2, double ( * energy) ( double r ) );

double C_N_CT_energy ( double theta );
double CT_C_N_energy ( double theta );
double CT_C_O_energy ( double theta );
double N_C_O_energy ( double theta );

double valence_energy ( int i1, int i2, int i3, double ( * energy) ( double phi ) );

double X2_N_C_X2_energy ( double phi );
double X2_C_N_X_energy ( double phi );

double dihedral_energy ( int i1, int i2, int i3, int i4, double ( * energy) ( double phi ) );

double E_total_bvd();

#endif
