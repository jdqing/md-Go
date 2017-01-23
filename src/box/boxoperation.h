#ifndef _BOXOEPRATION_H
#define _BOXOEPRATION_H

#include "../init/init.h"
#include "../tools/auxiliarytool.h"

void add_crowder ( int tp );
void remove_crowder ( int i );


int box_init();
int set_box ();

void treat_overlaps ( void );
void adjust_box ( void );
void print_overlaps ( void );
int check_box ( void );
void Bfactor ( char * file_Bfactor );

#endif
