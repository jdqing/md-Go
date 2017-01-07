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
    read_Simu_Config(fp);
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

  read_maxi_key("MAXIKEY.dat");//read data from force field
  printf(".");

  list_crowder_types();
  printf(".");

  char str_filename[50];
  sprintf(str_filename,"Config1_%s_T%d.dat",pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  fp = fopen (str_filename,"r");//Check the config1 file exist or not. This
                                //decide run from the beginning or continue run.
  if (fp!=NULL)
  {
    // printf("////////%s exist.\n",str_filename);
    fclose (fp);
  }
  else
  {
    // printf("////////%s does not exist.\n",str_filename);
    sprintf(str_filename,"%s.pdb",pSimuCfg->PDBID);
    read_PDB (str_filename);
    printf(".");

  }

  printf(".\n");
  printf("////////Initialization Complete\n");
}

int end_Go()
{

}
