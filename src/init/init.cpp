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

  LJ_SCALE_n = pSimuCfg->LJScaleN;
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

  int In1, In2, In12, InInt;

  free(pSimuCfg);

	free ( NS_atom );

	free ( RIGID_SET );

	free ( RIGID_END );

	for (In1 = 1; In1 < NB_atom + 1; In1++) free ( atom_key [ In1 ] );

	free ( atom_key );

	free ( amino_key );

	free ( crowder_key );

	free ( part_key );

	free ( maxi_key );

	free ( residue_mass );

	for (In1 = 0; In1 < 20; In1++) free ( residue_content [In1] );

	free ( residue_content );

	////////////////////////////////////////

	for (In1 = 0; In1 < n_crwd; In1++)
	{
		free ( RX_crwd [ In1 ] );

		free ( RY_crwd [ In1 ] );

		free ( RZ_crwd [ In1 ] );

		free ( part_key_crwd [ In1 ] );

		free ( maxi_key_crwd [ In1 ] );
	}

	free ( crowder_mass );

	for (In1 = 0; In1 < n_crwd; In1++) free ( crowder_content [In1] );

	free ( crowder_content );

	////////////////////////////////////////

	free ( CMS );

	free ( VCMS );

	////////////////////////////////////////

	free ( AXES );

	free ( W );

	////////////////////////////////////////

	free ( x ); free ( y ); free ( z );

	free ( x_old ); free ( y_old ); free ( z_old );

	free ( vx ); free ( vy ); free ( vz );

	free ( fx ); free ( fy ); free ( fz );

	////////////////////////////////////////
	free (x0);
  free (yy0);
  free (z0);

  free (x1);
  free (yy1);
  free (z1);

	////////////////////////////////////////

	free ( INDX ); free ( JNDX );

	free ( RMASS );

	free ( RVISC );

	free ( RXYZ );

	for (In1 = 1; In1 < NB_atom + 1; In1++) { free ( RX [ In1 ] ); free ( RY [ In1 ] ); free ( RZ [ In1 ] ); }

	free ( RX ); free ( RY ); free ( RZ );

	free ( IR1 ); free ( IR2 ); free ( IR3 );

	free ( AV ); free ( BV ); free ( CV );

	free ( FV ); free ( GV ); free ( HV );

	////////////////////////////////////////

	for (In12 = 0; In12 < nps; In12++)

		for (InInt = 0; InInt < npi; InInt++)

			free ( list_content [ In12 ][ InInt ] );

	////////////////////////////////////////

	for (In1 = 1; In1 < McT + 1; In1++) free ( cell_content [ In1 ] );

	free (cell_content); free (cell_mass);

	////////////////////////////////////////

	free ( bond_list_content );

	free ( valence_list_content );

	free ( dihedral_list_content );

	////////////////////////////////////////

	free ( native_list_content );

	free ( D2_LJ_n );

	free ( D2_LJ_CUTOFF_n );

	free ( E_LJ_n );

	free ( F_LJ_n );

	////////////////////////////////////////

	free ( THETA_NATIVE );

}
