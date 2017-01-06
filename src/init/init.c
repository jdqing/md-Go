/*
 * initialization part and ending part of the program
 */
#include "init.h"

SConfig *pSimuCfg;

int get_num (FILE *fp, char *flag, double *num){

}
int get_str (FILE *fp, char *flag, char   *str){

}
int get_bool(FILE *fp, char *flag, bool    *bl){

}
int read_Simu_Config(FILE *fp, SConfig *pCfg){

}

int init_Go(const char *cfgfile)
{
  printf("////Initialization of the whole program\n");
  printf("////////Initialization");

  FILE *fp;
  fp = fopen (cfgfile,"r");
  if (fp!=NULL)
  {
    pSimuCfg = ( SConfig * )calloc( sizeof ( SConfig ) );
    read_Simu_Config( fp, pSimuCfg);
    fclose (fp);
  }
  else{
    printf("\n////////Can not find file: %s\n",cfgfile );
    return 0;
  }

  printf(".\n");
  printf("////////Initialization Complete\n");
  return 1;
}

int end_Go()
{

}
