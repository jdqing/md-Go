/*
 * This is the init part and ending part of the project.
 */
#ifndef _INIT_H
#define _INIT_H

#include "glbcova.h"

typedef struct SimulationConfig {
  char          PDBID[20];
  double        TemperatureK;
  double        LJScaleN;
  double        CUTOFF;
  unsigned long SEED;
} SConfig;

extern SConfig *pSimuCfg;

int init_Go(const char *cfgfile);

int end_Go();

#endif
