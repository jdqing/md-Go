#include "../init/glbcova.h"

int Kdelta ( int i, int j );
void generate_perpendicular_vector ( double * n3, double * res );
void cross_product ( double * v1, double * v2, double * res );
void maxwell_1 ( double * array, int array_size );
void generate_atom_velocity ( int i, int j );
void move_atom ( int i, int j );
void leap_frog_move_atom ( int i, int j );
void half_shift ( double * xx );
void CC_vector ( int i1, int i2, double * xx );
double bond_length ( int i1, int i2 );
double valence_angle ( int i1, int i2, int i3 );
double dihedral_angle ( int i1, int i2, int i3, int i4 );
