#include "runGo.h"
int init_run_Go()
{
  int In1, In2, In12, InInt;
  char str_filename[50];

  printf("////Initialization of running\n");

  printf("////////Initialization\n");

  sprintf(str_filename, "Snap0_%s.pdb", pSimuCfg->PDBID);
  if(!read_original_structure ( str_filename )) return 0;
  displayProgress(30);

  set_box();
  displayProgress(60);

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
  displayProgress(90);

  // MTM_init();//molecular transfer model , ref.


  displayProgress(100);
  printf("////////Initialization complete\n");

  return 1;

}
int run_Go()
{
  if(!init_run_Go()) return 0;

  FILE *fp;

  char str_filename[50];

  long int klok_old = 0;

  const long int  total_steps = 50000000001;
  const int       steps_record= 10000;

  for (klok = klok_old; klok < total_steps; klok++ )
  {
    if(klok % steps_record == 0)
    {
      velocity_rescale ( sc );
      maxwell_tail = check_maxwell ();
      cm_maxwell_tail = check_cm_maxwell ();
      ang_maxwell_tail = check_angular_maxwell ();

      sprintf(str_filename, "Diagnost_%s_T%d.dat", pSimuCfg->PDBID,
              int(pSimuCfg->TemperatureK));
      fp = fopen (str_filename,"a");
      if (fp!=NULL)
      {
        fprintf ( fp, "%ld %le %le %le %le %le %le %le\n", klok, sc [3], sc [0],
              sc [1], sc [2], maxwell_tail, cm_maxwell_tail, ang_maxwell_tail );
        fclose (fp);
      }

    }

  }

  return 1;
}
int init_con_Go()
{

}
int con_Go()
{
  init_con_Go();
}
