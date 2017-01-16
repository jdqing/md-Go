#include "estimateGo.h"

int k1,k2;
int *k2_1;//this: k2 = k2_1[k1]; Just the map of k1->k2
int gen_k2_k1(const char *kfile){
  k2_1 = (int *)calloc(NTP_atom+1,sizeof(int));

  FILE *fp;
  fp = fopen (kfile,"r");
  if (fp!=NULL)
  {
    for (size_t i = 1; i < NTP_atom + 1; i++) {
      fscanf(fp,"%d %d",&k1,k2_1+i);
      // printf("%d %d\n",k1,k2_1[k1] );
    }
    fclose (fp);
    return 1;
  }
  else
  {
    printf("////////Can not find file: %s\n",kfile );
    return 0;
  }

}

double cal_ELJ_FF(const char *fileCoor){
  int i,j;

  FILE *fp;
  fp = fopen (fileCoor,"r");
  if (fp!=NULL)
  {
    char s1[10],s2[10],s3[10],s4[10];
    double aa,bb;
    int cc,countn=0;

    long position=0L,finish;
    fseek(fp,0L,SEEK_END);
    finish = ftell(fp);
    fseek(fp,0L,SEEK_SET);
    while(position<finish){
      countn++;

      for( i = 1; i <= NTP_atom; i++)
      {
        j = k2_1[i];//k1->k2
        fscanf(fp,"%s %s %s %s %d %lf %lf %lf %lf %lf",
                   s1,s2,s3,s4,&cc,x+j,y+j,z+j,&aa,&bb);
      }
      double E_total=0,x2,xx[3],r2,r4,r6;
      int i1,i2,j1;

      for(i=1;i<=NTP_atom;i++)
      {
        for(j=i+1;j<=NTP_atom;j++)
        {
          if( INDX[i] == INDX[j] ) continue;
          if( Flags2Atom[i][j]=='b' || Flags2Atom[i][j]=='a') continue;

          i1 = maxi_key[i];
          i2 = maxi_key[j];
          if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;
          else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

          xx[0] = x[i] - x[j];
          xx[1] = y[i] - y[j];
          xx[2] = z[i] - z[j];
          x2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

          if (x2 > D2_LJ_CUTOFF[j1])
          {
            r2 = D2_LJ [ j1 ] / x2;
            r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;
            r4 = E_LJ [ j1 ] * ( r2 - 2.0 * r6 );
            E_total+=r4;
          }
        }
      }
      FILE *fw;
      fw = fopen ("ELJ_FF_500.dat","a");
      if (fw!=NULL)
      {
        fprintf(fw,"%d %lf\n",countn,E_total);
        fclose (fw);
      }

      position = ftell(fp);
    }
    fclose (fp);
  }
  else
  {
    printf("////////Can not find file: %s\n",fileCoor );
  }
}

double cal_EFNLJ_FF(const char *fileCoor){
  int i,j;

  FILE *fp;
  fp = fopen (fileCoor,"r");
  if (fp!=NULL)
  {
    char s1[10],s2[10],s3[10],s4[10];
    double aa,bb;
    int cc,countn=0;
    double Eav_bvd=0,Eav_pc=0,Eav_xy1=0,Eav_xy2=0;

    long position=0L,finish;
    fseek(fp,0L,SEEK_END);
    finish = ftell(fp);
    fseek(fp,0L,SEEK_SET);
    while(position<finish){
      countn++;

      for( i = 1; i <= NTP_atom; i++)
      {
        j = k2_1[i];//k1->k2
        fscanf(fp,"%s %s %s %s %d %lf %lf %lf %lf %lf",
                   s1,s2,s3,s4,&cc,x+j,y+j,z+j,&aa,&bb);
      }

      set_box();
      int In1, In2, In12, InInt;

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

      double E_total_b=0,E_total_fnlj=0,E_total_xy1=0,E_total_xy2=0,x2,xx[3],r2,r4,r6;
      int i1,i2,j1,k;

      //E_bvd
      E_total_b = E_total_bvd();

      // E_nonnative_contact LJ, x2<D2_LJ_CUTOFF[j1]
      for(In12=0; In12 < nps; In12++)
      {
        for(k=1;k<list_mass[In12][0]+1;k+=2)
        {
          i = list_content[In12][0][k];
          j = list_content[In12][0][k+1];

          if( Flags2Atom[i][j]=='b' || Flags2Atom[i][j]=='a') continue;

          i1 = maxi_key[i];
          i2 = maxi_key[j];
          if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;
          else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

          xx[0] = x[i] - x[j];
          xx[1] = y[i] - y[j];
          xx[2] = z[i] - z[j];
          half_shift (xx);
          x2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

          if (x2 <= D2_LJ_CUTOFF[j1])
          {
            r2 = D2_LJ [ j1 ] / x2;
            r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;
            r4 = E_LJ [ j1 ] * ( r2 - 2.0 * r6 );
            E_total_fnlj+=r4;
          }
        }
      }


      for(i=1;i<=NTP_atom;i++)
      {
        for(j=i+1;j<=NTP_atom;j++)
        {
          if( INDX[i] == INDX[j] ) continue;
          if( Flags2Atom[i][j]=='b' || Flags2Atom[i][j]=='a') continue;

          i1 = maxi_key[i];
          i2 = maxi_key[j];
          if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;
          else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

          xx[0] = x[i] - x[j];
          xx[1] = y[i] - y[j];
          xx[2] = z[i] - z[j];
          x2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

          if (x2 > D2_LJ_CUTOFF[j1])
          {
            r2 = D2_LJ [ j1 ] / x2;
            r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;
            r4 = E_LJ [ j1 ] * ( r2 - 2.0 * r6 );
            E_total_xy1+=r4;
          }
        }
      }


      for(i=1;i<=NTP_atom;i++)
      {
        for(j=i+1;j<=NTP_atom;j++)
        {
          //if( INDX[i] == INDX[j] ) continue;
          if( Flags2Atom[i][j]=='b' || Flags2Atom[i][j]=='a') continue;

          i1 = maxi_key[i];
          i2 = maxi_key[j];
          if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;
          else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

          xx[0] = x[i] - x[j];
          xx[1] = y[i] - y[j];
          xx[2] = z[i] - z[j];
          x2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

          if (x2 > D2_LJ_CUTOFF[j1])
          {
            r2 = D2_LJ [ j1 ] / x2;
            r4 = r2 * r2; r6 = r4 * r2; r2 = r6 * r6;
            r4 = E_LJ [ j1 ] * ( r2 - 2.0 * r6 );
            E_total_xy2+=r4;
          }
        }
      }

      Eav_bvd=(Eav_bvd+E_total_b);
      Eav_pc =(Eav_pc +E_total_fnlj);
      Eav_xy1=(Eav_xy1+E_total_xy1);
      Eav_xy2=(Eav_xy2+E_total_xy2);

      // FILE *fw;
      // fw = fopen ("E_FF_2f21_300.dat","a");
      // if (fw!=NULL)
      // {
      //   fprintf(fw,"%d %lf %lf %lf %lf %lf %lf\n",countn,E_total_b,E_total_fnlj,
      //   E_total_xy1,E_total_xy2,E_total_b+E_total_fnlj+E_total_xy1,
      //   E_total_b+E_total_fnlj+E_total_xy2);
      //   fclose (fw);
      // }
      // FILE *fw;
      // fw = fopen ("Rc.dat","a");
      // if (fw!=NULL)
      // {
      //   double Q=remaining_native_contacts ();
      //   fprintf(fw,"%d %lf %lf\n",countn,radius_of_gyration ( 1, NTP_atom ), Q);
      //   fclose (fw);
      // }

      position = ftell(fp);
    }

    Eav_bvd=Eav_bvd/countn;
    Eav_pc =Eav_pc /countn;
    Eav_xy1=Eav_xy1/countn;
    Eav_xy2=Eav_xy2/countn;

    FILE *fw;
    fw = fopen ("E_FF_1pin.dat","a");
    if (fw!=NULL)
    {
      fprintf(fw,"E_FF_1PIN:\n");
      fprintf(fw,"%lf %lf %lf %lf %lf %lf %lf\n",Eav_bvd,Eav_pc,Eav_bvd+Eav_pc,
                Eav_xy1,Eav_xy2,Eav_bvd+Eav_pc+Eav_xy1,Eav_bvd+Eav_pc+Eav_xy2);
      fclose (fw);
    }
    fclose (fp);
  }
  else
  {
    printf("////////Can not find file: %s\n",fileCoor );
  }

}

double cal_ELJ_SBM( const char *fileCoor, double cutoff){

  int i,j, i1, i2, j1, k1, k;
  double xx[3],r2,x2;

  char str_filename[50];

  sprintf(str_filename,"%s.pdb",pSimuCfg->PDBID);
  read_original_structure ( str_filename );

  double Eav_na=0,Eav_non=0;

  native_list_mass = 0;

  for (i = 1; i < NTP_atom + 1; i++)
	{
		for (j = i + 1; j < NTP_atom + 1; j++)
		{
      if ( INDX [i] == INDX [j] ) continue;

      if( Flags2Atom[i][j]=='b' || Flags2Atom[i][j]=='a'|| Flags2Atom[i][j]=='d') continue;

      i1 = maxi_key [i];

			i2 = maxi_key [j];

			if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

			else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

      xx [0] = x [i] - x [j];

			xx [1] = y [i] - y [j];

			xx [2] = z [i] - z [j];

			r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];  // squre of the distance between two atoms in each pair

			if ( r2 < cutoff * cutoff * D2_LJ [j1] ) add_to_native_list ( i, j, r2, E_LJ [j1] );

    }
  }
  printf("//native_list_mass:%d\n", native_list_mass);

  FILE *fp;
  fp = fopen (fileCoor,"r");
  if (fp!=NULL)
  {
    char s1[10],s2[10],s3[10],s4[10];
    double aa,bb;
    int cc,countn=0;


    long position=0L,finish;
    fseek(fp,0L,SEEK_END);
    finish = ftell(fp);
    fseek(fp,0L,SEEK_SET);
    while(position<finish){
      countn++;

      double E_total_na = 0, E_total_non = 0;

      for( i = 1; i <= NTP_atom; i++)
      {
        j = k2_1[i];//k1->k2
        fscanf(fp,"%s %s %s %s %d %lf %lf %lf %lf %lf",
                   s1,s2,s3,s4,&cc,x+j,y+j,z+j,&aa,&bb);
      }

      set_box();
      int In1, In2, In12, InInt;

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

      for(In12=0; In12 < nps; In12++)
      {
        for(k=1;k<list_mass[In12][0]+1;k+=2)
        {
          i = list_content[In12][0][k];
          j = list_content[In12][0][k+1];

          if( Flags2Atom[i][j]=='b' || Flags2Atom[i][j]=='a') continue;

          i1 = maxi_key[i];
          i2 = maxi_key[j];
          if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;
          else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

          xx[0] = x[i] - x[j];
          xx[1] = y[i] - y[j];
          xx[2] = z[i] - z[j];
          // half_shift (xx);
          x2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

          E_total_non += WCA_ENERGY ( x2, j1 );


        }
      }



      for(k=1;k<native_list_mass + 1; k += 2)
      {
        i = native_list_content[k];
        j = native_list_content[k+1];

        xx [0] = x [i] - x [j];

  			xx [1] = y [i] - y [j];

  			xx [2] = z [i] - z [j];

        x2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

        i1 = maxi_key [i];

  			i2 = maxi_key [j];

  			if (i1 < i2) j1 = (i1 - 1) * na - (i1 - 1) * i1 / 2 + i2;

  			else j1 = (i2 - 1) * na - (i2 - 1) * i2 / 2 + i1;

        E_total_na-= WCA_ENERGY ( x2, j1 );



        j1 = (k+1)/2;
        if (x2 <= 0.81*D2_LJ_n [j1]) continue;


        E_total_na += LJ_ENERGY_n(x2,j1);
        // printf("%lf  %lf  %lf  %lf  %lf  %lf\n",x2,D2_LJ_n [j1],x2/D2_LJ_n [j1], LJ_ENERGY_n(x2,j1),D2_LJ_CUTOFF [j1],LJ_ENERGY(x2,j1));

      }

      Eav_na = Eav_na + E_total_na;
      Eav_non= Eav_non+ E_total_non;

      position = ftell(fp);
    }
    fclose (fp);

    Eav_na = Eav_na/countn;
    Eav_non = Eav_non/countn;

    FILE *fw;
    fw = fopen ("E_SBM_1pin.dat","a");
    if (fw!=NULL)
    {
      // fprintf(fw,"Eav_bvd,Eav_pc,Eav_xy1,Eav_xy2,Eav_all1,Eav_all2\n");
      fprintf(fw,"%lf %lf %lf\n",cutoff,Eav_na,Eav_non);
      fclose (fw);
    }
  }
  else
  {
    printf("////////Can not find file: %s\n",fileCoor );
  }



}

int estimate_Go()
{
  printf("////Estimation of parameters\n" );
  printf("////////Estimation begin\n" );

  //
  // gen_k2_k1 ("1pin_k2to1.dat");
  //
  // for(size_t i = 20; i <=40; i++)
  // {
  //   double cut = i * 0.05;
  //   cal_ELJ_SBM("1pinCoor300.pdb",cut);
  // }
  // for(size_t i = 20; i <=40; i++)
  // {
  //   double cut = i * 0.05;
  //   cal_ELJ_SBM("1pinCoor500.pdb",cut);
  // }



  //cal_ELJ_FF("1pinCoor500.pdb");

  // cal_EFNLJ_FF("1pinCoor300.pdb");//use
  //
  // cal_EFNLJ_FF("1pinCoor500.pdb");

  // cal_E_SBM("1pinCoor500.pdb",1.35);



  //get_E_FF();

  printf("////////Estimation complete\n");
  return 1;

}
