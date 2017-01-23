#include "readin.h"

SConfig *pSimuCfg = NULL;

bool get_words_after_keywords(char *get_words, const char *keywords, FILE *fp){
  long position, finish;
  char str[100];
  fseek( fp, 0L, SEEK_END);
  finish = ftell(fp);
  fseek( fp, 0L, SEEK_SET);
  do{
    fscanf(fp, "%s", str);
    if(strcmp(str, keywords)==0){
      fscanf(fp, "%s", get_words);
      return 1;
    }
    position = ftell(fp);
  }while( position < finish);

  return 0;
}
int get_long(FILE *fp, const char *flag, unsigned long *num){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    *num = (unsigned long)strtoul( str , NULL, 0);
    return 1;
  }
  return 0;
}
int get_dou (FILE *fp, const char *flag, double        *num){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    *num = atof( str );
    return 1;
  }
  return 0;
}
int get_str (FILE *fp, const char *flag, char       *strget){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    strcpy(strget,str);
    return 1;
  }
  return 0;
}
int get_bool(FILE *fp, const char *flag, bool           *bl){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    if(strcmp(str,"YES")==0||strcmp(str,"yes")==0||strcmp(str,"Yes")==0||
    strcmp(str,"Y")==0||strcmp(str,"y")==0){
      *bl = 1;
      return 1;
    }
    else if(strcmp(str,"NO")==0||strcmp(str,"no")==0||strcmp(str,"N")==0||
    strcmp(str,"n")==0){
      *bl = 0;
      return 1;
    }
  }
  return 0;
}
int read_Simu_Config(FILE *fp)
{
  pSimuCfg = ( SConfig * )calloc( 1, sizeof ( SConfig ) );

  get_str (fp, "PDBID",        pSimuCfg->PDBID);
  get_dou (fp, "TemperatureK",&pSimuCfg->TemperatureK);
  get_dou (fp, "LJScaleN",    &pSimuCfg->LJScaleN);
  get_dou (fp, "CUTOFF",      &pSimuCfg->CUTOFF);
  get_long(fp, "SEED",        &pSimuCfg->SEED);

  // printf("/////////////%s\n", pSimuCfg->PDBID);
  // printf("////////////%lf\n", pSimuCfg->TemperatureK);
  // printf("////////////%lf\n", pSimuCfg->LJScaleN);
  // printf("////////////%lf\n", pSimuCfg->CUTOFF);
  // printf("////////////%ld\n", pSimuCfg->SEED);

  return 1;
}
int read_maxi_key ( const char * file_maxi_key )
{
	FILE * f1;
	int Jn1, Jn2, Jn12;

	//////////////////////////////////////////////////////////////////////////////
	f1 = fopen ( file_maxi_key, "r" );

  if(f1==NULL) {printf("////////Can not find file MAXIKEY.dat\n");return 0;}

	fseek ( f1, 0L, SEEK_SET );

	for (Jn1 = 1; Jn1 < na + 1; Jn1++)
  fscanf ( f1, "%d %le %le %le\n", &Jn2, RLJ + Jn1, ELJ + Jn1, QC + Jn1 );

	fclose ( f1 );
	//////////////////////////////////////////////////////////////////////////////
	for (Jn1 = 1; Jn1 < na + 1; Jn1++)
	{
		for (Jn2 = Jn1; Jn2 < na + 1; Jn2++)
		{
			Jn12 = (Jn1 - 1) * na - (Jn1 - 1) * Jn1 / 2 + Jn2;

			D_LJ [ Jn12 ] = RLJ [ Jn1 ] + RLJ [ Jn2 ];

      if ( fabs ( D_LJ [ Jn12 ] ) < 1.0e-15 ) D_LJ [ Jn12 ] = 0;

			D2_LJ [ Jn12 ] = D_LJ [ Jn12 ] * D_LJ [ Jn12 ];

			D2_LJ_CUTOFF [ Jn12 ] = 1.0 * D2_LJ [ Jn12 ];//excluded volume interactions

			E_LJ [ Jn12 ] = sqrt ( ELJ [ Jn1 ] * ELJ [ Jn2 ] );

			if ( D2_LJ [ Jn12 ] ) F_LJ [ Jn12 ] = 12.0 * E_LJ [ Jn12 ] / D2_LJ [ Jn12 ];               //QM; value of 12 means?

			else F_LJ [ Jn12 ] = 0;
		}
	}
}

int read_PDB ( char * file_PDB )         // March 10, 2014
{
	FILE * f1;

	long position, finish;

	char str [10];

	char chain = 'A';

	int i, j, k, l, J1, J2, entry;

	double tmp1, tmp2, tmp3;

	double x_tmp, y_tmp, z_tmp;

	double A, B, C, F, G, H;

	double a0, a1, a2;

	double ep1 [3], ep2 [3], ep3 [3];

	/////////////////////////////////////////////////////////////////

	residue_mass = (int *) calloc ( 20, sizeof(int) );

	residue_content = (int **) calloc ( 20, sizeof(int *) );

	/////////////////////////////////////////////////////////////////

	f1 = fopen ( file_PDB, "r" );       // QM: file_PDB is the pdb file used in this program.

  if(f1==NULL) {printf("////////Can not find file %s\n",file_PDB);return 0;}

  fseek ( f1, 0L, SEEK_END );

	finish = ftell ( f1 );              // QM: finish = the value of the last pointer address in the given pdb file

	/////////////////////////////////////////////////////////////////

	fseek ( f1, 0L, SEEK_SET );

	k = 0;

	entry = 0;

	do
	{
		entry ++;

		fscanf ( f1, "%s", str );

		switch ( entry )
		{

		case 3:                  // QM:  entry==3, means read the " Name of Heavy Atom" in pdb file.

			if ( strcmp ( str, "N" ) == 0 )
			{
				k = 1;

				N_amino += 1;

				amino_key = (char *) realloc ( amino_key, (N_amino + 1) * sizeof(char) );
			}

			break;

		case 4:                 // QM:  entry==4, means read the " Name of Amino Acid.

			if (k)              // QM:  k==1, only when "case 3" has been implemented. After "case 4", K==0;
			{
				if ( strcmp ( str, "ALA" ) == 0 )
				{
					amino_key [ N_amino ] = 'A';   // QM: One key alphabet represents the amino acid

					/////////////////////////////////////////////////////////////////

					residue_mass [0] += 1;          // QM: statistic the number of each type amino acid

					residue_content [0] = (int *) realloc ( residue_content [0], (residue_mass [0] + 1) * sizeof(int) );

					residue_content [0][ residue_mass [0] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "ARG" ) == 0 )
				{
					amino_key [ N_amino ] = 'R';

					/////////////////////////////////////////////////////////////////

					residue_mass [1] += 1;

					residue_content [1] = (int *) realloc ( residue_content [1], (residue_mass [1] + 1) * sizeof(int) );

					residue_content [1][ residue_mass [1] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "ASN" ) == 0 )
				{
					amino_key [ N_amino ] = 'N';

					/////////////////////////////////////////////////////////////////

					residue_mass [2] += 1;

					residue_content [2] = (int *) realloc ( residue_content [2], (residue_mass [2] + 1) * sizeof(int) );

					residue_content [2][ residue_mass [2] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "ASP" ) == 0 )
				{
					amino_key [ N_amino ] = 'D';

					/////////////////////////////////////////////////////////////////

					residue_mass [3] += 1;

					residue_content [3] = (int *) realloc ( residue_content [3], (residue_mass [3] + 1) * sizeof(int) );

					residue_content [3][ residue_mass [3] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "CYS" ) == 0 )
				{
					amino_key [ N_amino ] = 'C';

					/////////////////////////////////////////////////////////////////

					residue_mass [4] += 1;

					residue_content [4] = (int *) realloc ( residue_content [4], (residue_mass [4] + 1) * sizeof(int) );

					residue_content [4][ residue_mass [4] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "GLU" ) == 0 )
				{
					amino_key [ N_amino ] = 'E';

					/////////////////////////////////////////////////////////////////

					residue_mass [5] += 1;

					residue_content [5] = (int *) realloc ( residue_content [5], (residue_mass [5] + 1) * sizeof(int) );

					residue_content [5][ residue_mass [5] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "GLN" ) == 0 )
				{
					amino_key [ N_amino ] = 'Q';

					/////////////////////////////////////////////////////////////////

					residue_mass [6] += 1;

					residue_content [6] = (int *) realloc ( residue_content [6], (residue_mass [6] + 1) * sizeof(int) );

					residue_content [6][ residue_mass [6] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "GLY" ) == 0 )
				{
					amino_key [ N_amino ] = 'G';

					/////////////////////////////////////////////////////////////////

					residue_mass [7] += 1;

					residue_content [7] = (int *) realloc ( residue_content [7], (residue_mass [7] + 1) * sizeof(int) );

					residue_content [7][ residue_mass [7] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "HIS" ) == 0 )
				{
					amino_key [ N_amino ] = 'H';

					/////////////////////////////////////////////////////////////////

					residue_mass [8] += 1;

					residue_content [8] = (int *) realloc ( residue_content [8], (residue_mass [8] + 1) * sizeof(int) );

					residue_content [8][ residue_mass [8] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "ILE" ) == 0 )
				{
					amino_key [ N_amino ] = 'I';

					/////////////////////////////////////////////////////////////////

					residue_mass [9] += 1;

					residue_content [9] = (int *) realloc ( residue_content [9], (residue_mass [9] + 1) * sizeof(int) );

					residue_content [9][ residue_mass [9] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "LEU" ) == 0 )
				{
					amino_key [ N_amino ] = 'L';

					/////////////////////////////////////////////////////////////////

					residue_mass [10] += 1;

					residue_content [10] = (int *) realloc ( residue_content [10], (residue_mass [10] + 1) * sizeof(int) );

					residue_content [10][ residue_mass [10] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "LYS" ) == 0 )
				{
					amino_key [ N_amino ] = 'K';

					/////////////////////////////////////////////////////////////////

					residue_mass [11] += 1;

					residue_content [11] = (int *) realloc ( residue_content [11], (residue_mass [11] + 1) * sizeof(int) );

					residue_content [11][ residue_mass [11] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "MET" ) == 0 )
				{
					amino_key [ N_amino ] = 'M';

					/////////////////////////////////////////////////////////////////

					residue_mass [12] += 1;

					residue_content [12] = (int *) realloc ( residue_content [12], (residue_mass [12] + 1) * sizeof(int) );

					residue_content [12][ residue_mass [12] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "PHE" ) == 0 )
				{
					amino_key [ N_amino ] = 'F';

					/////////////////////////////////////////////////////////////////

					residue_mass [13] += 1;

					residue_content [13] = (int *) realloc ( residue_content [13], (residue_mass [13] + 1) * sizeof(int) );

					residue_content [13][ residue_mass [13] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "PRO" ) == 0 )
				{
					amino_key [ N_amino ] = 'P';

					/////////////////////////////////////////////////////////////////

					residue_mass [14] += 1;

					residue_content [14] = (int *) realloc ( residue_content [14], (residue_mass [14] + 1) * sizeof(int) );

					residue_content [14][ residue_mass [14] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "SER" ) == 0 )
				{
					amino_key [ N_amino ] = 'S';

					/////////////////////////////////////////////////////////////////

					residue_mass [15] += 1;

					residue_content [15] = (int *) realloc ( residue_content [15], (residue_mass [15] + 1) * sizeof(int) );

					residue_content [15][ residue_mass [15] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "THR" ) == 0 )
				{
					amino_key [ N_amino ] = 'T';

					/////////////////////////////////////////////////////////////////

					residue_mass [16] += 1;

					residue_content [16] = (int *) realloc ( residue_content [16], (residue_mass [16] + 1) * sizeof(int) );

					residue_content [16][ residue_mass [16] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "TRP" ) == 0 )
				{
					amino_key [ N_amino ] = 'W';

					/////////////////////////////////////////////////////////////////

					residue_mass [17] += 1;

					residue_content [17] = (int *) realloc ( residue_content [17], (residue_mass [17] + 1) * sizeof(int) );

					residue_content [17][ residue_mass [17] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "TYR" ) == 0 )
				{
					amino_key [ N_amino ] = 'Y';

					/////////////////////////////////////////////////////////////////

					residue_mass [18] += 1;

					residue_content [18] = (int *) realloc ( residue_content [18], (residue_mass [18] + 1) * sizeof(int) );

					residue_content [18][ residue_mass [18] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				if ( strcmp ( str, "VAL" ) == 0 )
				{
					amino_key [ N_amino ] = 'V';

					/////////////////////////////////////////////////////////////////

					residue_mass [19] += 1;

					residue_content [19] = (int *) realloc ( residue_content [19], (residue_mass [19] + 1) * sizeof(int) );

					residue_content [19][ residue_mass [19] ] = N_amino;
				}

				/////////////////////////////////////////////////////////////////

				k = 0;
			}

			break;

		case 5:                  // QM: entry==5, read the chain name, A or B or C ...

			if ( strcmp ( str, "B" ) == 0 && chain == 'A' ) { TER1 = N_amino - 1; chain = 'B'; }

			if ( strcmp ( str, "C" ) == 0 && chain == 'B' ) { TER2 = N_amino - 1; chain = 'C'; }

			if ( strcmp ( str, "D" ) == 0 && chain == 'C' ) { TER3 = N_amino - 1; chain = 'D'; }

			break;

		case 12:   // entry==12, means the last arrow. Then entry==0, read the next line.

			entry = 0;

			break;

		default:

			break;
		}

		position = ftell ( f1 );  // QM: the position of the pointer "ָ��"
	}

	while ( position < finish );    // QM: finish == "the last position of the pointer"  ???

	/////////////////////////////////////////////////////////////////

	TER4 = N_amino;

	if ( chain == 'C' ) TER3 = TER4;

	if ( chain == 'B' ) { TER2 = TER4; TER3 = TER4; }

	if ( chain == 'A' ) { TER1 = TER4; TER2 = TER4; TER3 = TER4; }




	/////////////////////////////////////////////////////////////////

	NB_atom = N_amino;                // why the number of protein backbone is equal to that of the number of residues?????????????????????

	crowder_key = (int *) calloc ( NB_atom + 1, sizeof(int) );

	NS_atom = (int *) calloc ( NB_atom + 1, sizeof(int) );

	RIGID_SET = (int *) calloc ( NB_atom + 1, sizeof(int) );

	RIGID_END = (int *) calloc ( NB_atom + 1, sizeof(int) );

	atom_key = (int **) calloc ( NB_atom + 1, sizeof(int *) );

	/////////////////////////////////////////////////////////////////

	fseek ( f1, 0L, SEEK_SET );   ///// QM: Need not reopen the file, this statement can realize similar function.

	k = 0;

	entry = 0;

	do
	{
		entry ++;

		fscanf ( f1, "%s", str );

		switch ( entry )
		{

		case 3:

			if ( strcmp ( str, "N" ) == 0 )
			{
				k += 1;           // when read "N", k ++; That means "k" equals to the amino acid sequence in the protein chain

				/////////////////////////////////////////////////////////////////

				NT_atom += 1;      // QM: "NT_atom", means the total number of atoms in the protein

				/////////////////////////////////////////////////////////////////

				x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );              // QM: Why "Nt_atom +1" ? ???

				y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

				z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );
			}

			/////////////////////////////////////////////////////////////////

			else if ( strcmp ( str, "CA" ) == 0 )
			{
				NT_atom += 1;

				/////////////////////////////////////////////////////////////////

				x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );

				y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

				z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );
			}

			/////////////////////////////////////////////////////////////////

			else if ( strcmp ( str, "C" ) == 0 )
			{
				NT_atom += 1;

				/////////////////////////////////////////////////////////////////

				x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );

				y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

				z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );
			}

			/////////////////////////////////////////////////////////////////

			else if ( strcmp ( str, "O" ) == 0 )
			{
				NT_atom += 1;

				/////////////////////////////////////////////////////////////////

				x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );

				y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

				z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );

				/////////////////////////////////////////////////////////////////
                /// QM: maxi_key records all the types of the heavy atoms in protein

				maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );

				part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );    /// part_key

				part_key [ NT_atom - 3 ] = 0;

				part_key [ NT_atom - 2 ] = 1;

				part_key [ NT_atom - 1 ] = 1;

				part_key [ NT_atom ] = 2;

				/////////////////////////////////////////////////////////////////

				atom_key [k] = (int *) calloc ( 4, sizeof(int) );         /////////////////////QM:  atom_key[k][i], i=0,1,2,3,4....

				atom_key [k][0] = NT_atom;

				atom_key [k][1] = NT_atom - 1;

				atom_key [k][2] = NT_atom - 2;

				atom_key [k][3] = NT_atom - 3;

				/////////////////////////////////////////////////////////////////

				NS_atom [k] = 3;              // Original value (0-O, 1-C, 2-Ca, 3-N ????) before calculating the heavy atoms in each side chain of AA.

				switch ( amino_key [k] )      // QM: animo acid sequence in the ojection protein chain.
				{

				case 'A':

					NS_atom [k] += 1;         // calculate the number of heavy atoms in amino acid side chains.  but why originally add "3" ????? not "4" ????

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 1;

					maxi_key [ NT_atom - 2 ] = 2;

					maxi_key [ NT_atom - 1 ] = 3;

					maxi_key [ NT_atom ] = 4;

					break;

				case 'R':

					NS_atom [k] += 7;           // QM:  "7" is the number of heavy atoms in the side chain of amino acid Arg.

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 6;

					maxi_key [ NT_atom - 2 ] = 7;

					maxi_key [ NT_atom - 1 ] = 8;

					maxi_key [ NT_atom ] = 9;

					break;

				case 'N':

					NS_atom [k] += 4;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 16;

					maxi_key [ NT_atom - 2 ] = 17;

					maxi_key [ NT_atom - 1 ] = 18;

					maxi_key [ NT_atom ] = 19;

					break;

				case 'D':

					NS_atom [k] += 4;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 24;

					maxi_key [ NT_atom - 2 ] = 25;

					maxi_key [ NT_atom - 1 ] = 26;

					maxi_key [ NT_atom ] = 27;

					break;

				case 'C':

					NS_atom [k] += 2;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 31;

					maxi_key [ NT_atom - 2 ] = 32;

					maxi_key [ NT_atom - 1 ] = 33;

					maxi_key [ NT_atom ] = 34;

					break;

				case 'E':

					NS_atom [k] += 5;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 37;

					maxi_key [ NT_atom - 2 ] = 38;

					maxi_key [ NT_atom - 1 ] = 39;

					maxi_key [ NT_atom ] = 40;

					break;

				case 'Q':

					NS_atom [k] += 5;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 45;

					maxi_key [ NT_atom - 2 ] = 46;

					maxi_key [ NT_atom - 1 ] = 47;

					maxi_key [ NT_atom ] = 48;

					break;

				case 'G':

					NS_atom [k] += 0;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 54;

					maxi_key [ NT_atom - 2 ] = 55;

					maxi_key [ NT_atom - 1 ] = 56;

					maxi_key [ NT_atom ] = 57;

					break;

				case 'H':

					NS_atom [k] += 6;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 58;

					maxi_key [ NT_atom - 2 ] = 59;

					maxi_key [ NT_atom - 1 ] = 60;

					maxi_key [ NT_atom ] = 61;

					break;

				case 'I':

					NS_atom [k] += 4;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 68;

					maxi_key [ NT_atom - 2 ] = 69;

					maxi_key [ NT_atom - 1 ] = 70;

					maxi_key [ NT_atom ] = 71;

					break;

				case 'L':

					NS_atom [k] += 4;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 76;

					maxi_key [ NT_atom - 2 ] = 77;

					maxi_key [ NT_atom - 1 ] = 78;

					maxi_key [ NT_atom ] = 79;

					break;

				case 'K':

					NS_atom [k] += 5;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 83;

					maxi_key [ NT_atom - 2 ] = 84;

					maxi_key [ NT_atom - 1 ] = 85;

					maxi_key [ NT_atom ] = 86;

					break;

				case 'M':

					NS_atom [k] += 4;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 92;

					maxi_key [ NT_atom - 2 ] = 93;

					maxi_key [ NT_atom - 1 ] = 94;

					maxi_key [ NT_atom ] = 95;

					break;

				case 'F':

					NS_atom [k] += 7;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 100;

					maxi_key [ NT_atom - 2 ] = 101;

					maxi_key [ NT_atom - 1 ] = 102;

					maxi_key [ NT_atom ] = 103;

					break;

				case 'P':

					NS_atom [k] += 3;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 109;

					maxi_key [ NT_atom - 2 ] = 110;

					maxi_key [ NT_atom - 1 ] = 111;

					maxi_key [ NT_atom ] = 112;

					break;

				case 'S':

					NS_atom [k] += 2;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 116;

					maxi_key [ NT_atom - 2 ] = 117;

					maxi_key [ NT_atom - 1 ] = 118;

					maxi_key [ NT_atom ] = 119;

					break;

				case 'T':

					NS_atom [k] += 3;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 122;

					maxi_key [ NT_atom - 2 ] = 123;

					maxi_key [ NT_atom - 1 ] = 124;

					maxi_key [ NT_atom ] = 125;

					break;

				case 'W':

					NS_atom [k] += 10;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 129;

					maxi_key [ NT_atom - 2 ] = 130;

					maxi_key [ NT_atom - 1 ] = 131;

					maxi_key [ NT_atom ] = 132;

					break;

				case 'Y':

					NS_atom [k] += 8;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 143;

					maxi_key [ NT_atom - 2 ] = 144;

					maxi_key [ NT_atom - 1 ] = 145;

					maxi_key [ NT_atom ] = 146;

					break;

				case 'V':

					NS_atom [k] += 3;

					/////////////////////////////////////////////////////////////////

					maxi_key [ NT_atom - 3 ] = 153;

					maxi_key [ NT_atom - 2 ] = 154;

					maxi_key [ NT_atom - 1 ] = 155;

					maxi_key [ NT_atom ] = 156;

					break;

				default: break;

				}

				/////////////////////////////////////////////////////////////////    QM: notice their values.    ///////////////////////////////////////////////////////////

				RIGID_SET [k] = 1;

				RIGID_END [k] = NS_atom [k];

				/////////////////////////////////////////////////////////////////

				atom_key [k] = (int *) realloc ( atom_key [k], ( NS_atom [k] + 1 ) * sizeof(int) );
			}

			/////////////////////////////////////////////////////////////////

			else
			{
				NT_atom += 1;

				/////////////////////////////////////////////////////////////////

				x = (double *) realloc ( x, (NT_atom + 1) * sizeof(double) );

				y = (double *) realloc ( y, (NT_atom + 1) * sizeof(double) );

				z = (double *) realloc ( z, (NT_atom + 1) * sizeof(double) );

				/////////////////////////////////////////////////////////////////

				maxi_key = (int *) realloc ( maxi_key, (NT_atom + 1) * sizeof(int) );

				part_key = (int *) realloc ( part_key, (NT_atom + 1) * sizeof(int) );

				/////////////////////////////////////////////////////////////////

				switch ( amino_key [k] )                                             // QM: amino acid sequence in object protein chain
				{

				case 'A':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 5; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'R':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 10; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 11; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD" ) == 0 ) {
						maxi_key [ NT_atom ] = 12; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "NE" ) == 0 ) {
						maxi_key [ NT_atom ] = 13; part_key [ NT_atom ] = 0; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "CZ" ) == 0 ) {
						maxi_key [ NT_atom ] = 14; part_key [ NT_atom ] = 1; atom_key [k][8] = NT_atom; }

					else if ( strcmp ( str, "NH1" ) == 0 ) {
						maxi_key [ NT_atom ] = 15; part_key [ NT_atom ] = 0; atom_key [k][9] = NT_atom; }

					else if ( strcmp ( str, "NH2" ) == 0 ) {
						maxi_key [ NT_atom ] = 15; part_key [ NT_atom ] = 0; atom_key [k][10] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'N':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 20; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 21; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "OD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 22; part_key [ NT_atom ] = 2; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "ND2" ) == 0 ) {
						maxi_key [ NT_atom ] = 23; part_key [ NT_atom ] = 0; atom_key [k][7] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'D':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 28; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 29; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "OD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 30; part_key [ NT_atom ] = 2; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "OD2" ) == 0 ) {
						maxi_key [ NT_atom ] = 30; part_key [ NT_atom ] = 2; atom_key [k][7] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'C':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 35; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "SG" ) == 0 ) {
						maxi_key [ NT_atom ] = 36; part_key [ NT_atom ] = 4; atom_key [k][5] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'E':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 41; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 42; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD" ) == 0 ) {
						maxi_key [ NT_atom ] = 43; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "OE1" ) == 0 ) {
						maxi_key [ NT_atom ] = 44; part_key [ NT_atom ] = 2; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "OE2" ) == 0 ) {
						maxi_key [ NT_atom ] = 44; part_key [ NT_atom ] = 2; atom_key [k][8] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'Q':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 49; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 50; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD" ) == 0 ) {
						maxi_key [ NT_atom ] = 51; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "OE1" ) == 0 ) {
						maxi_key [ NT_atom ] = 52; part_key [ NT_atom ] = 2; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "NE2" ) == 0 ) {
						maxi_key [ NT_atom ] = 53; part_key [ NT_atom ] = 0; atom_key [k][8] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'G': break;

				/////////////////////////////////////////////////////////////////

				case 'H':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 62; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 63; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "ND1" ) == 0 ) {
						maxi_key [ NT_atom ] = 64; part_key [ NT_atom ] = 0; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CD2" ) == 0 ) {
						maxi_key [ NT_atom ] = 65; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "CE1" ) == 0 ) {
						maxi_key [ NT_atom ] = 66; part_key [ NT_atom ] = 1; atom_key [k][8] = NT_atom; }

					else if ( strcmp ( str, "NE2" ) == 0 ) {
						maxi_key [ NT_atom ] = 67; part_key [ NT_atom ] = 0; atom_key [k][9] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'I':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 72; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG1" ) == 0 ) {
						maxi_key [ NT_atom ] = 73; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CG2" ) == 0 ) {
						maxi_key [ NT_atom ] = 74; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 75; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'L':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 80; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 81; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 82; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CD2" ) == 0 ) {
						maxi_key [ NT_atom ] = 82; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'K':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 87; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 88; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD" ) == 0 ) {
						maxi_key [ NT_atom ] = 89; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CE" ) == 0 ) {
						maxi_key [ NT_atom ] = 90; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "NZ" ) == 0 ) {
						maxi_key [ NT_atom ] = 91; part_key [ NT_atom ] = 0; atom_key [k][8] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'M':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 96; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 97; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "SD" ) == 0 ) {
						maxi_key [ NT_atom ] = 98; part_key [ NT_atom ] = 4; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CE" ) == 0 ) {
						maxi_key [ NT_atom ] = 99; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'F':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 104; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 105; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 106; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CD2" ) == 0 ) {
						maxi_key [ NT_atom ] = 106; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "CE1" ) == 0 ) {
						maxi_key [ NT_atom ] = 107; part_key [ NT_atom ] = 1; atom_key [k][8] = NT_atom; }

					else if ( strcmp ( str, "CE2" ) == 0 ) {
						maxi_key [ NT_atom ] = 107; part_key [ NT_atom ] = 1; atom_key [k][9] = NT_atom; }

					else if ( strcmp ( str, "CZ" ) == 0 ) {
						maxi_key [ NT_atom ] = 108; part_key [ NT_atom ] = 1; atom_key [k][10] = NT_atom; }

					break;

				/////////////////////////////////////////////////////////////////

				case 'P':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 113; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 114; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD" ) == 0 ) {
						maxi_key [ NT_atom ] = 115; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'S':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 120; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "OG" ) == 0 ) {
						maxi_key [ NT_atom ] = 121; part_key [ NT_atom ] = 2; atom_key [k][5] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'T':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 126; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "OG1" ) == 0 ) {
						maxi_key [ NT_atom ] = 127; part_key [ NT_atom ] = 2; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CG2" ) == 0 ) {
						maxi_key [ NT_atom ] = 128; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'W':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 133; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 134; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 135; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CD2" ) == 0 ) {
						maxi_key [ NT_atom ] = 136; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "NE1" ) == 0 ) {
						maxi_key [ NT_atom ] = 137; part_key [ NT_atom ] = 0; atom_key [k][8] = NT_atom; }

					else if ( strcmp ( str, "CE2" ) == 0 ) {
						maxi_key [ NT_atom ] = 138; part_key [ NT_atom ] = 1; atom_key [k][9] = NT_atom; }

					else if ( strcmp ( str, "CE3" ) == 0 ) {
						maxi_key [ NT_atom ] = 139; part_key [ NT_atom ] = 1; atom_key [k][10] = NT_atom; }

					else if ( strcmp ( str, "CZ2" ) == 0 ) {
						maxi_key [ NT_atom ] = 140; part_key [ NT_atom ] = 1; atom_key [k][11] = NT_atom; }

					else if ( strcmp ( str, "CZ3" ) == 0 ) {
						maxi_key [ NT_atom ] = 141; part_key [ NT_atom ] = 1; atom_key [k][12] = NT_atom; }

					else if ( strcmp ( str, "CH2" ) == 0 ) {
						maxi_key [ NT_atom ] = 142; part_key [ NT_atom ] = 1; atom_key [k][13] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'Y':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 147; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG" ) == 0 ) {
						maxi_key [ NT_atom ] = 148; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CD1" ) == 0 ) {
						maxi_key [ NT_atom ] = 149; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else if ( strcmp ( str, "CD2" ) == 0 ) {
						maxi_key [ NT_atom ] = 149; part_key [ NT_atom ] = 1; atom_key [k][7] = NT_atom; }

					else if ( strcmp ( str, "CE1" ) == 0 ) {
						maxi_key [ NT_atom ] = 150; part_key [ NT_atom ] = 1; atom_key [k][8] = NT_atom; }

					else if ( strcmp ( str, "CE2" ) == 0 ) {
						maxi_key [ NT_atom ] = 150; part_key [ NT_atom ] = 1; atom_key [k][9] = NT_atom; }

					else if ( strcmp ( str, "CZ" ) == 0 ) {
						maxi_key [ NT_atom ] = 151; part_key [ NT_atom ] = 1; atom_key [k][10] = NT_atom; }

					else if ( strcmp ( str, "OH" ) == 0 ) {
						maxi_key [ NT_atom ] = 152; part_key [ NT_atom ] = 2; atom_key [k][11] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				/////////////////////////////////////////////////////////////////

				case 'V':

					if ( strcmp ( str, "CB" ) == 0 ) {
						maxi_key [ NT_atom ] = 157; part_key [ NT_atom ] = 1; atom_key [k][4] = NT_atom; }

					else if ( strcmp ( str, "CG1" ) == 0 ) {
						maxi_key [ NT_atom ] = 158; part_key [ NT_atom ] = 1; atom_key [k][5] = NT_atom; }

					else if ( strcmp ( str, "CG2" ) == 0 ) {
						maxi_key [ NT_atom ] = 158; part_key [ NT_atom ] = 1; atom_key [k][6] = NT_atom; }

					else printf ( "Unknown atom %s in res %d\n", str, k );

					break;

				default: break;

				}
			}

			break;

		case 7: x [ NT_atom ] = atof ( str ); break;

		case 8: y [ NT_atom ] = atof ( str ); break;

		case 9: z [ NT_atom ] = atof ( str ); break;

		case 12: entry = 0; break;          //when entry==0, read the next line.

		default: break;

		}

		position = ftell ( f1 );
	}

	while ( position < finish );           //"finish" is the value of the end of the file.

	fclose (f1);

	/////////////////////////////////////////////////////////////////

	IR1 = (double *) calloc ( NB_atom + 1, sizeof(double) );

	IR2 = (double *) calloc ( NB_atom + 1, sizeof(double) );

	IR3 = (double *) calloc ( NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	CMS = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	VCMS = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	AXES = (double *) calloc ( 9 * NB_atom + 1, sizeof(double) );

	W = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	x_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

	y_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

	z_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

	vx = (double *) calloc ( NT_atom + 1, sizeof(double) );

	vy = (double *) calloc ( NT_atom + 1, sizeof(double) );

	vz = (double *) calloc ( NT_atom + 1, sizeof(double) );

	fx = (double *) calloc ( NT_atom + 1, sizeof(double) );

	fy = (double *) calloc ( NT_atom + 1, sizeof(double) );

	fz = (double *) calloc ( NT_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	INDX = (int *) calloc ( NT_atom + 1, sizeof(int) );

	JNDX = (int *) calloc ( NT_atom + 1, sizeof(int) );

	for (i = 1; i < NB_atom + 1; i++)           // QM: amino acid sequence in protein chain
	{
		for (j = 0; j < NS_atom [i] + 1; j++)   // QM: heavy atom sequence in each amino acid in protein chain
		{
			k = atom_key [i][j];                // QM: k is the heavy atom position sequence in atom_key rule.

			INDX [k] = i;           // QM: given the atom sequence, deduce the amino acid sequence "i" and heavy atom sequence "j" in the "residue i"

			JNDX [k] = j;
		}
	}

	/////////////////////////////////////////////////////////////////

	RMASS = (double *) calloc ( NB_atom + 1, sizeof(double) );

	RVISC = (double *) calloc ( NB_atom + 1, sizeof(double) );

	RXYZ = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	RX = (double **) calloc ( NB_atom + 1, sizeof(double *) );

	RY = (double **) calloc ( NB_atom + 1, sizeof(double *) );

	RZ = (double **) calloc ( NB_atom + 1, sizeof(double *) );

	for (i = 1; i < NB_atom + 1; i++)
	{
		RX [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

		RY [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

		RZ [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );
	}

	AV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	BV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	CV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	FV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	GV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	HV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	for (i = 1; i < NB_atom + 1; i++)
	{
		for (j = 0; j < RIGID_SET [i]; j++) ///// QM: Notice the value of "RIGID_SET[i]". now still not clear!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

		  generate_atom_velocity ( i, j );

		for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++)

	      generate_atom_velocity ( i, j );

		if ( RIGID_SET [i] > RIGID_END [i] ) continue;

		/////////////////////////////////////////////////////////////////

		J1 = 3 * (i - 1);

		J2 = 9 * (i - 1);

		/////////////////////////////////////////////////////////////////

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)    ///////////// RIGID_SET[i]==1 ????
		{
			k = atom_key [i][j];

			l = part_key [k];

			/////////////////////////////////////////////////////////////////

			RMASS [i] += MASS [l];

			RVISC [i] += VISC [l];

			/////////////////////////////////////////////////////////////////

			CMS [J1 + 1] += MASS [l] * x [k];

			CMS [J1 + 2] += MASS [l] * y [k];

			CMS [J1 + 3] += MASS [l] * z [k];
		}

		CMS [J1 + 1] /= RMASS [i];            //// QM: CMS[J1+1] is the mass center of each residue in protein chain?????

		CMS [J1 + 2] /= RMASS [i];

		CMS [J1 + 3] /= RMASS [i];

		/////////////////////////////////

		A = 0; B = 0; C = 0; F = 0; G = 0; H = 0;

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			l = part_key [k];

			x_tmp = x [k] - CMS [J1 + 1];

			y_tmp = y [k] - CMS [J1 + 2];

			z_tmp = z [k] - CMS [J1 + 3];

			A += MASS [l] * (y_tmp * y_tmp + z_tmp * z_tmp);  //rotation inertia

			B += MASS [l] * (x_tmp * x_tmp + z_tmp * z_tmp);

			C += MASS [l] * (x_tmp * x_tmp + y_tmp * y_tmp);

			F += MASS [l] * y_tmp * z_tmp;

			G += MASS [l] * x_tmp * z_tmp;

			H += MASS [l] * x_tmp * y_tmp;
		}

		/////////////////////////////////

		a2 = - A - B - C;

		a1 = A * B + A * C + B * C - H * H - G * G - F * F;

		a0 = A * F * F + B * G * G + C * H * H + 2.0 * H * G * F - A * B * C;

		/////////////////////////////////

		tmp1 = - a2 / 3.0;

		x_tmp = ( 9.0 * a1 * a2 - 27 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;

		y_tmp = tmp1 * tmp1 - a1 / 3.0;

		z_tmp = x_tmp * pow ( y_tmp, - 1.5 );

		tmp2 = 2.0 * sqrt ( y_tmp );

		tmp3 = acos ( z_tmp ) / 3.0;

		a0 = pi / 3.0;

		/////////////////////////////////

		x_tmp = tmp1 + tmp2 * cos ( tmp3 );

		y_tmp = tmp1 - tmp2 * cos ( tmp3 - a0 );

		z_tmp = tmp1 - tmp2 * cos ( tmp3 + a0 );

		/////////////////////////////////

		if ( x_tmp < y_tmp ) { tmp1 = x_tmp; tmp2 = y_tmp; }

		else { tmp1 = y_tmp; tmp2 = x_tmp; }

		if ( z_tmp < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = z_tmp; }

		else if ( z_tmp < tmp2 ) { tmp3 = tmp2; tmp2 = z_tmp; }

		else tmp3 = z_tmp;

		IR1 [i] = tmp1;

		IR2 [i] = tmp2;

		IR3 [i] = tmp3;

		/////////////////////////////////

		x_tmp = (A - tmp1) * F + G * H;

		y_tmp = (B - tmp1) * G + F * H;

		z_tmp = (C - tmp1) * H + G * F;

		a0 = z_tmp / x_tmp;

		a1 = z_tmp / y_tmp;

		a2 = sqrt ( a0 * a0 + a1 * a1 + 1.0 );

		ep1 [0] = a0 / a2;

		ep1 [1] = a1 / a2;

		ep1 [2] = 1.0 / a2;

		/////////////////////////////////

		x_tmp = (A - tmp2) * F + G * H;

		y_tmp = (B - tmp2) * G + F * H;

		z_tmp = (C - tmp2) * H + G * F;

		a0 = x_tmp / y_tmp;

		a1 = x_tmp / z_tmp;

		a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

		ep2 [0] = 1.0 / a2;

		ep2 [1] = a0 / a2;

		ep2 [2] = a1 / a2;

		/////////////////////////////////

		cross_product ( ep1, ep2, ep3 );

		/////////////////////////////////

		AXES [J2 + 1] = ep1 [0]; AXES [J2 + 2] = ep1 [1]; AXES [J2 + 3] = ep1 [2];

		AXES [J2 + 4] = ep2 [0]; AXES [J2 + 5] = ep2 [1]; AXES [J2 + 6] = ep2 [2];

		AXES [J2 + 7] = ep3 [0]; AXES [J2 + 8] = ep3 [1]; AXES [J2 + 9] = ep3 [2];

		/////////////////////////////////

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			l = part_key [k];

			/////////////////////////////////

			x_tmp = x [k] - CMS [J1 + 1];

			y_tmp = y [k] - CMS [J1 + 2];

			z_tmp = z [k] - CMS [J1 + 3];

			/////////////////////////////////

			RX [i][j] = x_tmp * ep1 [0] + y_tmp * ep1 [1] + z_tmp * ep1 [2];

			RY [i][j] = x_tmp * ep2 [0] + y_tmp * ep2 [1] + z_tmp * ep2 [2];

			RZ [i][j] = x_tmp * ep3 [0] + y_tmp * ep3 [1] + z_tmp * ep3 [2];

			/////////////////////////////////

			if ( fabs ( RX [i][j] ) < 1.0e-10 ) RX [i][j] = 0;

			if ( fabs ( RY [i][j] ) < 1.0e-10 ) RY [i][j] = 0;

			if ( fabs ( RZ [i][j] ) < 1.0e-10 ) RZ [i][j] = 0;

			/////////////////////////////////

			RXYZ [J1 + 1] += VISC [l] * RX [i][j];

			RXYZ [J1 + 2] += VISC [l] * RY [i][j];

			RXYZ [J1 + 3] += VISC [l] * RZ [i][j];

			/////////////////////////////////

			AV [i] += VISC [l] * ( RY [i][j] * RY [i][j] + RZ [i][j] * RZ [i][j] );

			BV [i] += VISC [l] * ( RX [i][j] * RX [i][j] + RZ [i][j] * RZ [i][j] );

			CV [i] += VISC [l] * ( RX [i][j] * RX [i][j] + RY [i][j] * RY [i][j] );

			FV [i] += VISC [l] * RY [i][j] * RZ [i][j];

			GV [i] += VISC [l] * RX [i][j] * RZ [i][j];

			HV [i] += VISC [l] * RX [i][j] * RY [i][j];
		}

		/////////////////////////////////

		if ( fabs ( RXYZ [J1 + 1] ) < 1.0e-10 ) RXYZ [J1 + 1] = 0;

		if ( fabs ( RXYZ [J1 + 2] ) < 1.0e-10 ) RXYZ [J1 + 2] = 0;

		if ( fabs ( RXYZ [J1 + 3] ) < 1.0e-10 ) RXYZ [J1 + 3] = 0;

		/////////////////////////////////

genw1:
		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw1;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR1 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 1] = tmp1 * cos ( tmp2 );

genw2:
		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw2;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR2 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 2] = tmp1 * cos ( tmp2 );

genw3:
		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw3;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR3 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 3] = tmp1 * cos ( tmp2 );

genvcms1:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms1;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 1] = tmp1 * cos ( tmp2 );

genvcms2:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms2;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 2] = tmp1 * cos ( tmp2 );

genvcms3:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms3;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 3] = tmp1 * cos ( tmp2 );
	}

	/////////////////////////////////////////////////////////////////

	NBP_atom = NB_atom;

	NTP_atom = NT_atom;

  return 1;
}

int read_original_structure  ( char * file_original_structure )
{

	FILE * f1;

	long position, finish;

	char str [10];

	//char chain = 'A';

	int i, i_0, j, k, l, J1, J2, entry;

	double tmp1, tmp2, tmp3;

	double x_tmp, y_tmp, z_tmp;

	double A, B, C, F, G, H;

	double a0, a1, a2;

	double ep1 [3], ep2 [3], ep3 [3];


	/////////////////////////////////////////////////////////////////

	f1 = fopen ( file_original_structure, "r" );// QM: file_original_structure
            //is the origina structure with pdb style in the folding trajectory.

  if(f1==NULL){printf("////////Can not find file %s\n",file_original_structure);
                return 0;}

	fseek ( f1, 0L, SEEK_END );

	finish = ftell ( f1 );// QM: finish = the value of the last pointer address
                        //in the given original structure file

	//////////////////////////////////////////////////////////////////////////////

	fseek ( f1, 0L, SEEK_SET );   ///// QM: Need not reopen the file,
                                //this statement can realize similar function.

	k = 0;

	entry = 0;

	i_0 = 1;

	do
	{
		entry ++;

		fscanf ( f1, "%s", str );

		switch ( entry )
		{

		case 7: x [ i_0 ] = atof ( str ); break;   // read the X, Y, Z coordinates

		case 8: y [ i_0 ] = atof ( str ); break;

		case 9: z [ i_0 ] = atof ( str ); break;

		case 12: entry = 0; i_0 ++; break;          //when entry==0, read the next line.

		default: break;

		}

		position = ftell ( f1 );
	}

	while ( position < finish );           //"finish" is the value of the end of the file.

	fclose (f1);

	/////////////////////////////////////////////////////////////////

	IR1 = (double *) calloc ( NB_atom + 1, sizeof(double) );

	IR2 = (double *) calloc ( NB_atom + 1, sizeof(double) );

	IR3 = (double *) calloc ( NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	CMS = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	VCMS = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	AXES = (double *) calloc ( 9 * NB_atom + 1, sizeof(double) );

	W = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	x_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

	y_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

	z_old = (double *) calloc ( NT_atom + 1, sizeof(double) );

	vx = (double *) calloc ( NT_atom + 1, sizeof(double) );

	vy = (double *) calloc ( NT_atom + 1, sizeof(double) );

	vz = (double *) calloc ( NT_atom + 1, sizeof(double) );

	fx = (double *) calloc ( NT_atom + 1, sizeof(double) );

	fy = (double *) calloc ( NT_atom + 1, sizeof(double) );

	fz = (double *) calloc ( NT_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	INDX = (int *) calloc ( NT_atom + 1, sizeof(int) );

	JNDX = (int *) calloc ( NT_atom + 1, sizeof(int) );

	for (i = 1; i < NB_atom + 1; i++)           // QM: amino acid sequence in protein chain
	{
		for (j = 0; j < NS_atom [i] + 1; j++)   // QM: heavy atom sequence in each amino acid in protein chain
		{
			k = atom_key [i][j];                // QM: k is the heavy atom position sequence in atom_key rule.

			INDX [k] = i;           // QM: given the atom sequence, deduce the amino acid sequence "i" and heavy atom sequence "j" in the "residue i"

			JNDX [k] = j;
		}
	}

	/////////////////////////////////////////////////////////////////

	RMASS = (double *) calloc ( NB_atom + 1, sizeof(double) );

	RVISC = (double *) calloc ( NB_atom + 1, sizeof(double) );

	RXYZ = (double *) calloc ( 3 * NB_atom + 1, sizeof(double) );

	RX = (double **) calloc ( NB_atom + 1, sizeof(double *) );

	RY = (double **) calloc ( NB_atom + 1, sizeof(double *) );

	RZ = (double **) calloc ( NB_atom + 1, sizeof(double *) );

	for (i = 1; i < NB_atom + 1; i++)
	{
		RX [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

		RY [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );

		RZ [i] = (double *) calloc ( NS_atom [i] + 1, sizeof(double) );
	}

	AV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	BV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	CV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	FV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	GV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	HV = (double *) calloc ( NB_atom + 1, sizeof(double) );

	/////////////////////////////////////////////////////////////////

	for (i = 1; i < NB_atom + 1; i++)
	{
		for (j = 0; j < RIGID_SET [i]; j++) ///// QM: Notice the value of "RIGID_SET[i]". now still not clear!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

		  generate_atom_velocity ( i, j );

		for (j = RIGID_END [i] + 1; j < NS_atom [i] + 1; j++)

	      generate_atom_velocity ( i, j );

		if ( RIGID_SET [i] > RIGID_END [i] ) continue;

		/////////////////////////////////////////////////////////////////

		J1 = 3 * (i - 1);

		J2 = 9 * (i - 1);

		/////////////////////////////////////////////////////////////////

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)    ///////////// RIGID_SET[i]==1 ????
		{
			k = atom_key [i][j];

			l = part_key [k];

			/////////////////////////////////////////////////////////////////

			RMASS [i] += MASS [l];

			RVISC [i] += VISC [l];

			/////////////////////////////////////////////////////////////////

			CMS [J1 + 1] += MASS [l] * x [k];

			CMS [J1 + 2] += MASS [l] * y [k];

			CMS [J1 + 3] += MASS [l] * z [k];
		}

		CMS [J1 + 1] /= RMASS [i];            //// QM: CMS[J1+1] is the mass center of each residue in protein chain?????

		CMS [J1 + 2] /= RMASS [i];

		CMS [J1 + 3] /= RMASS [i];

		/////////////////////////////////

		A = 0; B = 0; C = 0; F = 0; G = 0; H = 0;

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			l = part_key [k];

			x_tmp = x [k] - CMS [J1 + 1];

			y_tmp = y [k] - CMS [J1 + 2];

			z_tmp = z [k] - CMS [J1 + 3];

			A += MASS [l] * (y_tmp * y_tmp + z_tmp * z_tmp);  //rotation inertia

			B += MASS [l] * (x_tmp * x_tmp + z_tmp * z_tmp);

			C += MASS [l] * (x_tmp * x_tmp + y_tmp * y_tmp);

			F += MASS [l] * y_tmp * z_tmp;

			G += MASS [l] * x_tmp * z_tmp;

			H += MASS [l] * x_tmp * y_tmp;
		}

		/////////////////////////////////

		a2 = - A - B - C;

		a1 = A * B + A * C + B * C - H * H - G * G - F * F;

		a0 = A * F * F + B * G * G + C * H * H + 2.0 * H * G * F - A * B * C;

		/////////////////////////////////

		tmp1 = - a2 / 3.0;

		x_tmp = ( 9.0 * a1 * a2 - 27 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;

		y_tmp = tmp1 * tmp1 - a1 / 3.0;

		z_tmp = x_tmp * pow ( y_tmp, - 1.5 );

		tmp2 = 2.0 * sqrt ( y_tmp );

		tmp3 = acos ( z_tmp ) / 3.0;

		a0 = pi / 3.0;

		/////////////////////////////////

		x_tmp = tmp1 + tmp2 * cos ( tmp3 );

		y_tmp = tmp1 - tmp2 * cos ( tmp3 - a0 );

		z_tmp = tmp1 - tmp2 * cos ( tmp3 + a0 );

		/////////////////////////////////

		if ( x_tmp < y_tmp ) { tmp1 = x_tmp; tmp2 = y_tmp; }

		else { tmp1 = y_tmp; tmp2 = x_tmp; }

		if ( z_tmp < tmp1 ) { tmp3 = tmp2; tmp2 = tmp1; tmp1 = z_tmp; }

		else if ( z_tmp < tmp2 ) { tmp3 = tmp2; tmp2 = z_tmp; }

		else tmp3 = z_tmp;

		IR1 [i] = tmp1;

		IR2 [i] = tmp2;

		IR3 [i] = tmp3;

		/////////////////////////////////

		x_tmp = (A - tmp1) * F + G * H;

		y_tmp = (B - tmp1) * G + F * H;

		z_tmp = (C - tmp1) * H + G * F;

		a0 = z_tmp / x_tmp;

		a1 = z_tmp / y_tmp;

		a2 = sqrt ( a0 * a0 + a1 * a1 + 1.0 );

		ep1 [0] = a0 / a2;

		ep1 [1] = a1 / a2;

		ep1 [2] = 1.0 / a2;

		/////////////////////////////////

		x_tmp = (A - tmp2) * F + G * H;

		y_tmp = (B - tmp2) * G + F * H;

		z_tmp = (C - tmp2) * H + G * F;

		a0 = x_tmp / y_tmp;

		a1 = x_tmp / z_tmp;

		a2 = sqrt ( 1.0 + a0 * a0 + a1 * a1 );

		ep2 [0] = 1.0 / a2;

		ep2 [1] = a0 / a2;

		ep2 [2] = a1 / a2;

		/////////////////////////////////

		cross_product ( ep1, ep2, ep3 );

		/////////////////////////////////

		AXES [J2 + 1] = ep1 [0]; AXES [J2 + 2] = ep1 [1]; AXES [J2 + 3] = ep1 [2];

		AXES [J2 + 4] = ep2 [0]; AXES [J2 + 5] = ep2 [1]; AXES [J2 + 6] = ep2 [2];

		AXES [J2 + 7] = ep3 [0]; AXES [J2 + 8] = ep3 [1]; AXES [J2 + 9] = ep3 [2];

		/////////////////////////////////

		for (j = RIGID_SET [i]; j < RIGID_END [i] + 1; j++)
		{
			k = atom_key [i][j];

			l = part_key [k];

			/////////////////////////////////

			x_tmp = x [k] - CMS [J1 + 1];

			y_tmp = y [k] - CMS [J1 + 2];

			z_tmp = z [k] - CMS [J1 + 3];

			/////////////////////////////////

			RX [i][j] = x_tmp * ep1 [0] + y_tmp * ep1 [1] + z_tmp * ep1 [2];

			RY [i][j] = x_tmp * ep2 [0] + y_tmp * ep2 [1] + z_tmp * ep2 [2];

			RZ [i][j] = x_tmp * ep3 [0] + y_tmp * ep3 [1] + z_tmp * ep3 [2];

			/////////////////////////////////

			if ( fabs ( RX [i][j] ) < 1.0e-10 ) RX [i][j] = 0;

			if ( fabs ( RY [i][j] ) < 1.0e-10 ) RY [i][j] = 0;

			if ( fabs ( RZ [i][j] ) < 1.0e-10 ) RZ [i][j] = 0;

			/////////////////////////////////

			RXYZ [J1 + 1] += VISC [l] * RX [i][j];

			RXYZ [J1 + 2] += VISC [l] * RY [i][j];

			RXYZ [J1 + 3] += VISC [l] * RZ [i][j];

			/////////////////////////////////

			AV [i] += VISC [l] * ( RY [i][j] * RY [i][j] + RZ [i][j] * RZ [i][j] );

			BV [i] += VISC [l] * ( RX [i][j] * RX [i][j] + RZ [i][j] * RZ [i][j] );

			CV [i] += VISC [l] * ( RX [i][j] * RX [i][j] + RY [i][j] * RY [i][j] );

			FV [i] += VISC [l] * RY [i][j] * RZ [i][j];

			GV [i] += VISC [l] * RX [i][j] * RZ [i][j];

			HV [i] += VISC [l] * RX [i][j] * RY [i][j];
		}

		/////////////////////////////////

		if ( fabs ( RXYZ [J1 + 1] ) < 1.0e-10 ) RXYZ [J1 + 1] = 0;

		if ( fabs ( RXYZ [J1 + 2] ) < 1.0e-10 ) RXYZ [J1 + 2] = 0;

		if ( fabs ( RXYZ [J1 + 3] ) < 1.0e-10 ) RXYZ [J1 + 3] = 0;

		/////////////////////////////////

genw1:
		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw1;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR1 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 1] = tmp1 * cos ( tmp2 );

genw2:
		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw2;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR2 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 2] = tmp1 * cos ( tmp2 );

genw3:
		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genw3;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / IR3 [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		W [J1 + 3] = tmp1 * cos ( tmp2 );

genvcms1:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms1;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 1] = tmp1 * cos ( tmp2 );

genvcms2:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms2;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 2] = tmp1 * cos ( tmp2 );

genvcms3:

		a1 = (double) rand () / RAND_MAX; if ( a1 == 0 ) goto genvcms3;

		a2 = (double) rand () / RAND_MAX;

		tmp1 = sqrt ( T / RMASS [i] ) * sqrt ( -2.0 * log (a1) );

		tmp2 = 2.0 * pi * a2;

		VCMS [J1 + 3] = tmp1 * cos ( tmp2 );
	}

	/////////////////////////////////////////////////////////////////

	NBP_atom = NB_atom;

	NTP_atom = NT_atom;

  return 1;
}
