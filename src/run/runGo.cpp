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

#if MTM_yn
  MTM_init(NTP_atom,pSimuCfg->CONCENT);//molecular transfer model , ref.
#endif

  displayProgress(100);
  printf("////////Initialization complete\n");

  return 1;

}
int run_Go()
{
  if(!init_run_Go()) return 0;

  printf("////MD run\n" );
  printf("////////MD run begin\n" );

  FILE *fp;

  int In1;

  long int klok_old = 0;

  // const long int  total_steps = 500000001;
  // const int       steps_record= 10000;
  long       total_steps = pSimuCfg->T_STEPS;
  long       steps_record= pSimuCfg->STEPS_RE;

  for (klok = klok_old; klok <= total_steps; klok++ )
  {
    if(klok % steps_record == 0)
    {
      make_records();
    }

    char str_filename[50];
    sprintf(str_filename, "Contacts_%s_T%d.dat", pSimuCfg->PDBID,
            int(pSimuCfg->TemperatureK));
    deterministic_forces ( str_filename );

    full_forces ();

#if MTM_yn
    MTM_ene(NTP_atom);
#endif

		move_rigid_units ();

		In1 = check_box ();
		if ( !In1 ) check_shifts ();

  }
  printf("////////MD run complete\n" );

  return 1;
}
int init_con_Go()
{

}
int con_Go()
{
  init_con_Go();
}
