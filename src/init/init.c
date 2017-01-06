/*
 * initialization part and ending part of the program
 */
#include "init.h"

int init_Go(const char *cfgfile)
{
  printf("////Initialization of the whole program\n");
  printf("////////Initialization\n");

  FILE *fp;//read config file begin.
  fp = fopen (cfgfile,"r");
  if (fp!=NULL)
  {
    pSimuCfg = ( SConfig * )calloc( 1, sizeof ( SConfig ) );
    read_Simu_Config( fp, pSimuCfg);
    fclose (fp);
  }
  else{
    printf("\n////////Can not find file: %s\n",cfgfile );
    return 0;
  }
  printf(".");//read config file completed.

  T = 0.596 * pSimuCfg->TemperatureK / 300.0;

  srand(pSimuCfg->SEED);

  box_init();//set box parameters:BORDER_MIN, BORDER_MAX.
  printf(".");//box parameters set completed.

  force_init();//sigma force pre factors . Ref. Langevin Dynamics.
  printf(".");//force factors preparation completed.

  set_dr1();
  printf(".");



  printf(".\n");
  printf("////////Initialization Complete\n");
}

int end_Go()
{

}
