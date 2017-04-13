#ifndef _READIN_H
#define _READIN_H

#include "../init/glbcova.h"
#include "../tools/auxiliarytool.h"

typedef struct SimulationConfig {
  char          PDBID[20];
  double        TemperatureK;
  double        LJScaleN;
  double        CUTOFF;
  unsigned long SEED;
  unsigned long T_STEPS;
  unsigned long STEPS_RE;

#if MTM_yn
  double        CONCENT;
#endif

} SConfig;

extern SConfig *pSimuCfg;


int read_Simu_Config( FILE * fp );

int read_maxi_key( const char * file_maxi_key );

int read_PDB ( char * file_PDB ) ;

int read_original_structure  ( char * file_original_structure );

#endif
