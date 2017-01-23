#ifndef _FORCE_H
#define _FORCE_H

#include "../init/glbcova.h"
#include "../tools/auxiliarytool.h"
#include "../fun/potential.h"

int force_init();

void CC_force_single_pair ( int i1, int i2, int Jn12, double scale_factor,
                            double ( * force) ( double x2, int Jn12 ) );

double LJ_FORCE ( double x2, int Jn12 );

double LJ_FORCE_n ( double x2, int Jn12 );

void LJ_interactions_lists ( void );

void bond_force ( int i1, int i2, double ( * force) ( double r ) );

void valence_force ( int i1, int i2, int i3, double ( * force) ( double phi ) );

void dihedral_force ( int i1, int i2, int i3, int i4, double ( * force) ( double phi ) );

void polymer_constraints ( void );

void LJ_interactions_native_lists ( char * file_contacts );

void deterministic_forces ( char * file_contacts );

void full_forces ( void );

#endif
