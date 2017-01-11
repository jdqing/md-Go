#include "runGo.h"
int init_run_Go()
{
  int In1, In2, In12, InInt;
  char str_filename[50];

  printf("////Initialization of running a coordinate\n");

  printf("////////Initialization\n");

  sprintf(str_filename, "Snap0_%s.pdb", pSimuCfg->PDBID);
  if(!read_original_structure ( str_filename )) return 0;

  set_box();
  displayProgress(40);

  for (In1 = 0; In1 < ns; In1++)
  {
    for (In2 = In1; In2 < ns; In2++)
    {
      In12 = In1 * ns - In1 * (1 + In1) / 2 + In2;

      for (InInt = 0; InInt < NPI; InInt++)
      {
        list_mass [ In12 ][ InInt ] = 0;

        list_content [ In12 ][ InInt ] = NULL;

        populate_list ( In1, In2, InInt );

        PROGRESS [ In12 ][ InInt ] = 0;

        // printf ("list_mass [%d][%d] = %d\n\n", In12, InInt, list_mass [In12][InInt] );
      }
    }
  }
  displayProgress(80);

  MTM_init();//molecular transfer model , ref.


  displayProgress(100);
  printf("////////Initialization complete\n");

  return 1;

}
int run_Go()
{
  if(!init_run_Go()) return 0;

  return 1;
}
int init_con_Go()
{

}
int con_Go()
{
  init_con_Go();
}
