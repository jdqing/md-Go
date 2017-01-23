#ifndef _WRITEOUT_H
#define _WRITEOUT_H

#include "../init/init.h"
#include "../diagnostics/diagnostics.h"
#include "../datadeal/dealfor.h"


void write_coordinates ( char * file_coordinates );

void write_configuration ( char * file_config1, char * file_config2 );

void write_protein_and_crowders ( char * file_PDB, char * frame_root, int frame_num );

int make_records();

#endif
