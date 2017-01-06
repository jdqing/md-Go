#ifndef _READIN_H
#define _READIN_H

#include "../init/glbcova.h"

typedef struct SimulationConfig {
  char          PDBID[20];
  double        TemperatureK;
  double        LJScaleN;
  double        CUTOFF;
  unsigned long SEED;
} SConfig;

extern SConfig *pSimuCfg;


int read_Simu_Config(FILE *fp, SConfig *pCfg);

#endif
