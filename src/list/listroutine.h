#ifndef _LISTROUTINE_H
#define _LISTROUTINE_H

#include "../init/glbcova.h"
#include "../tools/auxiliarytool.h"

int set_dr1();

int list_crowder_types ();

void add_bond ( int i1, int i2 );
void add_valence ( int i1, int i2 );
void add_dihedral ( int i1, int i2 );

int populate_bvd_lists ();

void add_to_native_list ( int i, int j, double r2, double E );

int populate_native_lists ( char * file_native_list, double cutoff );

int populate_list ( int In1, int In2, int InInt );

int separate_in_cells ( int In1, int In2, int InInt );

#endif
