#include "writeout.h"

void write_coordinates ( char * file_coordinates )
{
	FILE * f1;

	/////////////////////////////////////////////////

	f1 = fopen ( file_coordinates, "a" );

	/////////////////////////////////////////////////

	fwrite ( x, sizeof(double), NT_atom + 1, f1 );

	fwrite ( y, sizeof(double), NT_atom + 1, f1 );

	fwrite ( z, sizeof(double), NT_atom + 1, f1 );

	/////////////////////////////////////////////////

	fclose (f1);
}

void write_configuration ( char * file_config1, char * file_config2 )
{
	FILE * f1;

	int i;

	/////////////////////////////////////////////////

	f1 = fopen ( file_config1, "w" );

	fseek ( f1, 0L, SEEK_SET );

	fprintf ( f1, "%d %d %d %d %d %d %d %d %d %le %le %le %le\n", N_amino, NBP_atom, NTP_atom, NB_atom, NT_atom, TER1, TER2, TER3, TER4, side_X, side_Y, side_Z, side_XYZ );

	fclose (f1);

	/////////////////////////////////////////////////

	f1 = fopen ( file_config2, "w" );

	fseek ( f1, 0L, SEEK_SET );

	/////////////////////////////////////////////////

	fwrite ( residue_mass, sizeof(int), 20, f1 );

	/////////////////////////////////////////////////

	for (i = 0; i < 20; i++)
	{
		if ( residue_mass [i] )

			fwrite ( residue_content [i], sizeof(int), residue_mass [i] + 1, f1 );
	}

	/////////////////////////////////////////////////

	fwrite ( crowder_mass, sizeof(int), n_crwd, f1 );

	/////////////////////////////////////////////////

	for (i = 0; i < n_crwd; i++)
	{
		if ( crowder_mass [i] )

			fwrite ( crowder_content [i], sizeof(int), crowder_mass [i] + 1, f1 );
	}

	/////////////////////////////////////////////////

	fwrite ( NS_atom, sizeof(int), NB_atom + 1, f1 );

	/////////////////////////////////////////////////////////////////

	fwrite ( RIGID_SET, sizeof(int), NB_atom + 1, f1 );

	/////////////////////////////////////////////////////////////////

	fwrite ( RIGID_END, sizeof(int), NB_atom + 1, f1 );

	/////////////////////////////////////////////////////////////////

	for (i = 1; i < NB_atom + 1; i++) fwrite ( atom_key [i], sizeof(int), NS_atom [i] + 1, f1 );

	/////////////////////////////////////////////////////////////////

	fwrite ( amino_key, sizeof(char), N_amino + 1, f1 );

	/////////////////////////////////////////////////////////////////

	fwrite ( crowder_key, sizeof(int), NB_atom + 1, f1 );

	/////////////////////////////////////////////////////////////////

	fwrite ( maxi_key, sizeof(int), NT_atom + 1, f1 );

	/////////////////////////////////////////////////

	fwrite ( part_key, sizeof(int), NT_atom + 1, f1 );

	/////////////////////////////////////////////////

	fwrite ( IR1, sizeof(double), NB_atom + 1, f1 );

	fwrite ( IR2, sizeof(double), NB_atom + 1, f1 );

	fwrite ( IR3, sizeof(double), NB_atom + 1, f1 );

	/////////////////////////////////////////////////

	fwrite ( CMS, sizeof(double), 3 * NB_atom + 1, f1 );

	fwrite ( VCMS, sizeof(double), 3 * NB_atom + 1, f1 );

	/////////////////////////////////////////////////

	fwrite ( AXES, sizeof(double), 9 * NB_atom + 1, f1 );

	fwrite ( W, sizeof(double), 3 * NB_atom + 1, f1 );

	/////////////////////////////////////////////////

	fwrite ( x, sizeof(double), NT_atom + 1, f1 );

	fwrite ( y, sizeof(double), NT_atom + 1, f1 );

	fwrite ( z, sizeof(double), NT_atom + 1, f1 );

	fwrite ( x_old, sizeof(double), NT_atom + 1, f1 );

	fwrite ( y_old, sizeof(double), NT_atom + 1, f1 );

	fwrite ( z_old, sizeof(double), NT_atom + 1, f1 );

	fwrite ( vx, sizeof(double), NT_atom + 1, f1 );

	fwrite ( vy, sizeof(double), NT_atom + 1, f1 );

	fwrite ( vz, sizeof(double), NT_atom + 1, f1 );

	/////////////////////////////////////////////////

	fclose (f1);
}

void write_protein_and_crowders ( char * file_PDB, char * frame_root, int frame_num )
{
	FILE * f1;

	long position, finish;

	char str [10], file_frame [200], syscall [200];

	int i, j, k, l, entry;

	/////////////////////////////////////////////////////////////////

	strcpy ( file_frame, frame_root );

	sprintf ( file_frame + strlen ( file_frame ), "%04d.pdb", frame_num ); //QM: change 03d to 04d

	sprintf ( syscall, "cp %s %s", file_PDB, file_frame );  //////////////// whats mean of "cp" ????????????????????????????????????????? QM
    //sprintf ( syscall, "copy %s %s", file_PDB, file_frame );  //////////////// use new command xcopy or copy by QM ????????????????????????????????????????? QM
	system ( syscall );

	/////////////////////////////////////////////////////////////////
	//4 11 14 20 22 26 38 46 54 60 66

	f1 = fopen ( file_frame, "r+" );

	fseek ( f1, 0L, SEEK_END );

	finish = ftell ( f1 );

	/////////////////////////////////////////////////////////////////

	fseek ( f1, 0L, SEEK_SET );

	k = 1;

	entry = 0;

	do
	{
		entry ++;

		fscanf ( f1, "%s", str );

		switch ( entry )
		{

		case 6:

			position = ftell ( f1 );

			fseek ( f1, position, SEEK_SET );

			fprintf ( f1, "%12.3f%8.3f%8.3f", x [k], y [k], z [k] );

			//fseek ( f1, ftell ( f1 ), SEEK_SET ); //////////////////////////////// added by QM
			break;

		//case 8:// last column - 3
		case 9:

			entry = 0;

			k++;

			break;

		default: break;

		}

		position = ftell ( f1 );

		//printf("%ld,%ld,%ld,%ld \n",frame_num,k,position,finish);//////////////////////////////////////????????????????????????//
	}

	while ( position < finish );

	/////////////////////////////////////////////////////////////////

	fprintf ( f1, "\n" );

	//fprintf ( f1, "\nTER%8d%15s\n", k, "GLY D 356" );

	/////////////////////////////////////////////////////////////////

	entry = k;

	for (i = NBP_atom + 1; i < NB_atom + 1; i++)     /// QM: write x y z coordinates with pdb style for crowders
	{
		for (j = 0; j < NS_atom [i] + 1; j++)
		{
			entry ++;

			k = atom_key [i][j];

			l = part_key [k];

			/////////////////////////////////////////////////////////////////

			if ( l == 0 ) strcpy ( str, "N" );

			if ( l == 1 ) strcpy ( str, "C" );

			if ( l == 2 ) strcpy ( str, "O" );

			if ( l == 3 ) strcpy ( str, "H" );

			if ( l == 4 ) strcpy ( str, "S" );

			if ( l == 5 ) strcpy ( str, "K" );

			/////////////////////////////////////////////////////////////////

			fprintf ( f1, "HETATM%5d%3s%6s%2c%4d%12.3f%8.3f%8.3f\n", entry, str, str, 'E', i - NBP_atom, x [k], y [k], z [k] );
		}
	}

	/////////////////////////////////////////////////////////////////

	fclose (f1);
}
int make_records()
{
  double sc [4];
	double maxwell_tail, cm_maxwell_tail, ang_maxwell_tail;
  char str_filename[50],str_filename2[50];
  FILE *fp;

	static clock_t cl1;
	if(klok == 0)
	{
		cl1 = clock();
	}

	sprintf(str_filename, "Time_%s_T%d.dat", pSimuCfg->PDBID,
					int(pSimuCfg->TemperatureK));
	fp = fopen (str_filename,"a");
	fprintf(fp,"%ld %ld\n",klok,clock()-cl1);
	fclose(fp);
	cl1 = clock();

  // record diagnostics
  velocity_rescale ( sc );
  maxwell_tail = check_maxwell ();
  cm_maxwell_tail = check_cm_maxwell ();
  ang_maxwell_tail = check_angular_maxwell ();
  sprintf(str_filename, "Diagnost_%s_T%d.dat", pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  fp = fopen (str_filename,"a");
  fprintf ( fp, "%ld %le %le %le %le %le %le %le\n", klok, sc [3], sc [0],
        sc [1], sc [2], maxwell_tail, cm_maxwell_tail, ang_maxwell_tail );
  fclose (fp);

  // snapshot
  // sprintf(str_filename, "%s.pdb", pSimuCfg->PDBID);
  // sprintf(str_filename2, "Snap_%s_T%d_", pSimuCfg->PDBID,
  //         int(pSimuCfg->TemperatureK));
  // if ( (klok % 500000) == 0 )
  // write_protein_and_crowders ( str_filename, str_filename2, klok / 500000 );

  // record energy
  // sprintf(str_filename, "Energy_%s_T%d.dat", pSimuCfg->PDBID,
  //         int(pSimuCfg->TemperatureK));
  // fp = fopen (str_filename,"a");
  // fprintf ( fp, "%ld, %le, %le \n", klok,total_energy_LJ(),
  //           total_energy_LJ2() );
  // fclose (fp);

  sprintf(str_filename, "Dimensions_%s_T%d.dat", pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  fp = fopen (str_filename,"a");
  fprintf ( fp, "%ld %le %le\n", klok, radius_of_gyration ( 1, NTP_atom ),
            remaining_native_contacts () );
  fclose (fp);

  sprintf(str_filename,"Config1_%s_T%d.dat",pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  sprintf(str_filename2,"Config2_%s_T%d.dat",pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  write_configuration(str_filename, str_filename2);

  sprintf(str_filename, "Coordinates_%s_T%d.dat", pSimuCfg->PDBID,
          int(pSimuCfg->TemperatureK));
  write_coordinates(str_filename);
}
