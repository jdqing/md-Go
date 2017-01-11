/*
 * initialization part and ending part of the program
 */
#include "init.h"

int init_Go(const char *cfgfile)
{
  char str_filename[50];

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
  displayProgress(10);//read config file completed.

  T = 0.596 * pSimuCfg->TemperatureK / 300.0;

  srand(pSimuCfg->SEED);

  box_init();//set box parameters:BORDER_MIN, BORDER_MAX.
  displayProgress(20);//box parameters set completed.

  force_init();//sigma force pre factors . Ref. Langevin Dynamics.
  displayProgress(30);//force factors preparation completed.

  set_dr1();
  displayProgress(40);

  if(!read_maxi_key("MAXIKEY.dat")) return 0;//read data from force field
  displayProgress(50);

  list_crowder_types();
  displayProgress(60);


  sprintf(str_filename,"Config1_%s_T%d.dat",pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  fp = fopen (str_filename,"r");//Check the config1 file exist or not. This
                                //decide run from the beginning or continue run.
  if (fp!=NULL)
  {
    //Have not written this continue part that well ,for the moment.
    printf("////////%s exist.\n",str_filename);
    fclose (fp);
    return 0;
  }
  else
  {
    // printf("////////%s does not exist.\n",str_filename);
    sprintf(str_filename,"%s.pdb",pSimuCfg->PDBID);
    if(!read_PDB (str_filename)) return 0;
    displayProgress(70);

    populate_bvd_lists();//populate bvd = bond valence dihedral.
    displayProgress(80);

    sprintf(str_filename, "%s_noh.dat", pSimuCfg->PDBID);
    populate_native_lists(str_filename, pSimuCfg->CUTOFF);
    displayProgress(90);
  }

  flag_any_2_atoms();
  displayProgress(99);

  displayProgress(100);
  printf("////////Initialization Complete\n");
  return 1;
}

int end_Go()
{

}
