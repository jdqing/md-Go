/*
 * This is a tool for using the F90 Code in the C code project.
 */
#include "gbsa_f2c.h"

// common from F90
double surften = 0.005;// used for generalized Born
double *tran_ene_guan = NULL;
int nrespa = 1;
double rgbmax = 25.0;
double gb_fs_max = 0.0;

double *foldx,*foldy,*foldz;
double *oldx,*oldy,*oldz;
int *jj;
double *r2x;

double  * gbsafrc;
double  * vdwrad;
double  * lcpo_p1;
double  * lcpo_p2;
double  * lcpo_p3;
double  * lcpo_p4;

int * ineighbor;
int * j_vals;
int * k_vals;
const int na_tmp = 260;
double VDWR[na_tmp+1],LCPOP1[na_tmp+1],LCPOP2[na_tmp+1],LCPOP3[na_tmp+1],LCPOP4[na_tmp+1];

void MTM_readmaxi(const char * file_maxi_key) {

  FILE *f1;
  int Jn1,Jn2;
  char str1[10],str2[10],str3[10],str4[10];
  f1 = fopen ( file_maxi_key, "r" );

  if(f1==NULL) {printf("////////Can not find file lcpo_maix.dat\n");}

	fseek ( f1, 0L, SEEK_SET );

	for (Jn1 = 1; Jn1 < na + 1; Jn1++){
  fscanf ( f1, "%d %s %s %s %lf %s %lf %lf %lf %lf\n", &Jn2,str1,str2,str3,
          VDWR+Jn1,str4,LCPOP1+Jn1,LCPOP2+Jn1,LCPOP3+Jn1,LCPOP4+Jn1  );
  // printf ( "%d,%s,%s,%s,%lf,%s,%lf,%lf,%lf,%lf\n", Jn2,str1,str2,str3,
  //                 VDWR[Jn1],str4,LCPOP1[Jn1],LCPOP2[Jn1],LCPOP3[Jn1],LCPOP4[Jn1]  );
  }



	fclose ( f1 );
}

void MTM_init(int atm_cnt,const double Concent)
{

  // Local variables:
  int i,nmk;
  // num_ints and num_reals are used to return allocation counts. Don't zero.

  gbsafrc = (double *) calloc( 3*( atm_cnt + 1 ) , sizeof ( double ) );// *

  lcpo_p1 = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  lcpo_p2 = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  lcpo_p3 = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  lcpo_p4 = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  vdwrad  = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  j_vals  = (int *) calloc( atm_cnt + 1 , sizeof ( int ) );
  k_vals  = (int *) calloc( atm_cnt + 1 , sizeof ( int ) );
  ineighbor = (int *) calloc( (atm_cnt + 1) * 40 , sizeof ( int ) );// *
  foldx = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  foldy = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  foldz = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  oldx = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  oldy = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );
  oldz = (double *) calloc( atm_cnt + 1 , sizeof ( double ) );


  jj = (int *)calloc(atm_cnt+1,sizeof(int));
  r2x = (double *)calloc(atm_cnt+1,sizeof(double));


  MTM_readmaxi("lcpo_maxi.dat");

  for (i = 1; i <= atm_cnt; i++) {

    nmk = maxi_key[i];
    vdwrad[i] = VDWR[nmk] + 1.4;
    lcpo_p1[i]= LCPOP1[nmk];
    lcpo_p2[i]= LCPOP2[nmk];
    lcpo_p3[i]= LCPOP3[nmk];
    lcpo_p4[i]= LCPOP4[nmk];
  }

  double BB_mk = -39.21, BB_bk =-31.86;
  double SC_mk [20] ={   0.0,  -7.20, -41.77,  -75.99,  -68.10,  -85.86, -124.57, -50.86, -18.75, -18.75,
		                -102.03, -56.9, -123.19, -196.25, -102.03, -56.90, -65.0,   -67.95,  42.34, -65.0 };  // The last value is Cys, not sure!!!!!

	double SC_bk [20] ={ 0.0, -2.28,  -23.41, -41.42,   -37.93,  -42.76, -61.12, -27.76, -31.25, -31.25,
		              -65.73, -57.07, -78.71, -138.75,  -65.73,  -57.07, -85.0,   0.0,    0.0,    -85.0 };    // The last value is Cys, not sure!!!!!

	double BB_GkG [20] ={ 84.9817, 62.5336, 53.8005, 50.2648, 50.2648, 50.2648, 48.3400, 56.8578, 60.9227, 56.2249,
                          55.6039, 52.0904, 47.3327, 43.7780, 56.6455, 53.0339, 51.3509, 48.3400, 46.1820, 57.7220
	                    };

	double SC_GkG [20] ={  0.000,  108.259, 147.128, 164.683, 164.683, 164.683, 174.605, 132.713, 114.940, 135.721,
	                      138.647, 155.429, 179.916, 198.715, 133.722, 150.786, 159.160, 174.605, 185.982, 128.640
	                     };
	char amino_name [20] ={'G', 'A', 'V', 'L', 'I', 'M', 'F', 'P', 'S', 'T', 'N', 'Q', 'Y', 'W', 'D', 'E', 'H', 'K', 'R', 'C' };

  int i_aa,j_aa,j;

	char AKey;

  tran_ene_guan = (double *) calloc ( atm_cnt+1, sizeof(double) );

  for (i=1; i< atm_cnt+1; i++)
  {
      i_aa = INDX [i];
      j_aa = JNDX [i];
      AKey = amino_key[i_aa];

      for(j=0; j<20; j++)
      {
          if(AKey == amino_name[j])
          {
              if(AKey == 'G')
              {
                  tran_ene_guan[i] = (BB_mk * Concent + BB_bk )/BB_GkG[j]/1000.0;

                  break;
              }
              else
              {
                  if(j_aa<=3)
                      tran_ene_guan[i] = (BB_mk * Concent + BB_bk )/BB_GkG[j]/1000.0;
                  else
                      tran_ene_guan[i] = (SC_mk [j] * Concent + SC_bk[j] )/SC_GkG[j]/1000.0;

                  break;
              }

          }
      }

  }



}

void MTM_ene( int atm_cnt)
{
  double esurf=0.0;
  double frespa;
  double rgbmaxpsmax2;
  double xij, yij, zij;
  double xi, yi, zi;
  double r2;
  double dij;
  int max_i;
  int neibr_cnt;
  int icount;
  int i, j, k;



  frespa = nrespa;
  neibr_cnt = 0; // we have no neighbors yet for LCPO gbsa
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)*(rgbmax + gb_fs_max);
  max_i = atm_cnt;
  // if (natbel > 0) max_i = natbel;

  for (i = 1; i<= max_i; i++) {
    xi = x[i];
    yi = y[i];
    zi = z[i];

    icount = 0;
    for (j = 1; j<= atm_cnt; j++) {

      if ( i == j ) continue;

      xij = xi - x[j];
      yij = yi - y[j];
      zij = zi - z[j];
      r2 = xij * xij + yij * yij + zij * zij;
      if (r2 > rgbmaxpsmax2) continue;

      // pairlist contains only atoms within rgbmax + safety margin

      icount = icount + 1;
      jj[icount] = j;
      r2x[icount] = r2;

    }

    for (k = 1; k<= icount; k++) {
      j = jj[k];
      dij = sqrt(r2x[k]);

      if (vdwrad[i] + vdwrad[j] > dij) {
        // Only count LCPO for non-Hydrogens
        // if (vdwrad[i] > 2.5 && vdwrad[j] > 2.5) {
          neibr_cnt = neibr_cnt + 1;
          ineighbor[neibr_cnt] = j;
          // }
      }
    }

    neibr_cnt = neibr_cnt + 1;
    ineighbor[neibr_cnt] = 0;

  }
  clock_t c11;
  c11 = clock();
  MTM_lcpo(neibr_cnt, max_i, frespa);
// printf("%ld\n",clock()-c11);
}

void MTM_lcpo( int max_count, int max_i, double frespa)
{

  double esurf=0.0;
  int nei_count;  //neighbor count
  int i, j, k;    //atom indices
  int num_j_vals;
  int num_k_vals;
  int icount;
  int count2;
  int count2_fin;
  int j_loop, k_loop;

  double surf_i;
  double xi, yi, zi;
  double xj, yj, zj;
  double xk, yk, zk;
  double xji, yji, zji;
  double xjk, yjk, zjk;
  double aij, ajk;
  double rij, one_rij;
  double rjk, one_rjk;
  double p3p4aij;
  double vdw2dif;
  double lastxj;
  double lastyj;
  double lastzj;
  double ai;
  double totsasa;

//Running sums (i.e. sum a_ij)
  double sumaij;
  double sumajk;
  double sumajk2;
  double sumaijajk;
  double sumdaijddijdxi;
  double sumdaijddijdyi;
  double sumdaijddijdzi;
  double sumdaijddijdxiajk;
  double sumdaijddijdyiajk;
  double sumdaijddijdziajk;
  double sumdajkddjkdxj;
  double sumdajkddjkdyj;
  double sumdajkddjkdzj;
//Derivatives of pairwise overlaps (i.e., d a_ij / d d_ij)
  double daijddij;
  double daijddijdxj;
  double daijddijdyj;
  double daijddijdzj;
  double daidxi;
  double daidyi;
  double daidzi;
  double daidxj;
  double daidyj;
  double daidzj;
  double dajkddjk;
  double dajkddjkdxj;
  double dajkddjkdyj;
  double dajkddjkdzj;

  nei_count = 1;
  totsasa = 0;

  for( i = 1 ; i <= max_i; i++){
    foldx[i] = fx[i];
    foldy[i] = fy[i];
    foldz[i] = fz[i];
    oldx[i] = x[i];
    oldy[i] = y[i];
    oldz[i] = z[i];
  }


  for( i = 1 ; i <= max_i; i++){

    if( ineighbor[nei_count] == 0) nei_count++;
    else{
// printf("%d    %lf    %lf    %lf   ",i,fx[i],fy[i],fz[i]);
// if(klok % 1000 == 0)
// printf("%d    ",i);
// printf("%d    %lf    %lf    %lf   ",i,foldx[i],foldy[i],foldz[i]);
      surf_i = 4.0 * pi * vdwrad[i] * vdwrad[i];
      sumaij            = 0.0;
      sumajk            = 0.0;
      sumaijajk         = 0.0;
      sumdaijddijdxi    = 0.0;
      sumdaijddijdyi    = 0.0;
      sumdaijddijdzi    = 0.0;
      sumdaijddijdxiajk = 0.0;
      sumdaijddijdyiajk = 0.0;
      sumdaijddijdziajk = 0.0;

      icount = nei_count;
      num_j_vals = 0;

      do {
        j = ineighbor[nei_count];
        num_j_vals = num_j_vals + 1;
        j_vals[num_j_vals] = j;
        nei_count = nei_count + 1;

        if (ineighbor[nei_count] != 0) continue;

        nei_count = nei_count + 1;
        break;
      } while(1);

      xi = x[i];
      yi = y[i];
      zi = z[i];

      for (j_loop = 1; j_loop <= num_j_vals; j_loop++) {

        j = j_vals[j_loop];

        xj = x[j];
        yj = y[j];
        zj = z[j];

        xji = xj - xi;
        yji = yj - yi;
        zji = zj - zi;
        rij = sqrt(xji * xji + yji * yji + zji * zji);

        one_rij = 1 / rij;

        // Equation 3 of Weiser, et. al.; J Comput Chem, v20 p217
        aij = 2.0 * pi * vdwrad[i] * ( vdwrad[i] - rij * 0.5 -
              (vdwrad[i] * vdwrad[i] - vdwrad[j] * vdwrad[j]) *
              one_rij * 0.50);

        // First derivatives: Appendix A of Weiser, et. al.
        daijddij = pi * vdwrad[i] * (one_rij *one_rij *
                   (vdwrad[i] * vdwrad[i] - vdwrad[j] * vdwrad[j]) - 1);
        daijddijdxj = daijddij * xji * one_rij;
        daijddijdyj = daijddij * yji * one_rij;
        daijddijdzj = daijddij * zji * one_rij;

        // Add it up
        sumaij = sumaij + aij;

        count2 = icount;

        sumajk2 = 0.0;
        sumdajkddjkdxj = 0.0;
        sumdajkddjkdyj = 0.0;
        sumdajkddjkdzj = 0.0;

        p3p4aij = - (lcpo_p3[i] + lcpo_p4[i] * aij) * frespa;

        for (k_loop = count2; k_loop <= max_count; k_loop++) {
          if(ineighbor[k_loop] == 0){count2_fin = k_loop - 1;break;}
        }

        num_k_vals = 0;
        for (k_loop = count2; k_loop <= count2_fin; k_loop++) {
          k = ineighbor[k_loop];
          if(k != j){
            num_k_vals = num_k_vals + 1;
            k_vals[num_k_vals] = k;
          }
        }

        count2 = icount;
        for (k_loop = 1; k_loop <= num_k_vals; k_loop++) {

          k = k_vals[k_loop];
          xk = x[k];
          yk = y[k];
          zk = z[k];

          xjk = xj - xk;
          yjk = yj - yk;
          zjk = zj - zk;
          rjk = sqrt(xjk * xjk + yjk * yjk + zjk * zjk);

          one_rjk = 1.0 / rjk;

          if (vdwrad[j] + vdwrad[k] > rjk) {

            vdw2dif = vdwrad[j] * vdwrad[j] - vdwrad[k] * vdwrad[k];

            ajk = pi * vdwrad[j] *
                 (2.0 * vdwrad[j] - rjk - vdw2dif * one_rjk);

            sumajk = sumajk + ajk;
            sumajk2 = sumajk2 + ajk;

            // First derivatives

              dajkddjk = pi * vdwrad[j] * one_rjk *
                    (one_rjk * one_rjk * vdw2dif - 1.0);

              dajkddjkdxj = dajkddjk * xjk;
              dajkddjkdyj = dajkddjk * yjk;
              dajkddjkdzj = dajkddjk * zjk;

              sumdajkddjkdxj = sumdajkddjkdxj + dajkddjkdxj;
              sumdajkddjkdyj = sumdajkddjkdyj + dajkddjkdyj;
              sumdajkddjkdzj = sumdajkddjkdzj + dajkddjkdzj;

              double tmp1,tmp2,tmp3;
              tmp1=dajkddjkdxj * p3p4aij;
              tmp2=dajkddjkdyj * p3p4aij;
              tmp3=dajkddjkdzj * p3p4aij;
              // if(fabs(tmp1)>1)tmp1=0;
              // if(fabs(tmp2)>1)tmp2=0;
              // if(fabs(tmp3)>1)tmp3=0;


              fx[k] = fx[k] - tran_ene_guan [k] *tmp1;
              fy[k] = fy[k] - tran_ene_guan [k] *tmp2;
              fz[k] = fz[k] - tran_ene_guan [k] *tmp3;
          }
        }

        sumaijajk = sumaijajk+aij*sumajk2;

        // First derivatives

        sumdaijddijdxi = sumdaijddijdxi - daijddijdxj;
        sumdaijddijdyi = sumdaijddijdyi - daijddijdyj;
        sumdaijddijdzi = sumdaijddijdzi - daijddijdzj;

        sumdaijddijdxiajk = sumdaijddijdxiajk - daijddijdxj * sumajk2;
        sumdaijddijdyiajk = sumdaijddijdyiajk - daijddijdyj * sumajk2;
        sumdaijddijdziajk = sumdaijddijdziajk - daijddijdzj * sumajk2;

        lastxj = daijddijdxj * sumajk2 + aij * sumdajkddjkdxj;
        lastyj = daijddijdyj * sumajk2 + aij * sumdajkddjkdyj;
        lastzj = daijddijdzj * sumajk2 + aij * sumdajkddjkdzj;

        daidxj =  (lcpo_p2[i] * daijddijdxj + lcpo_p3[i] *
                            sumdajkddjkdxj + lcpo_p4[i] * lastxj);
        daidyj =  (lcpo_p2[i] * daijddijdyj + lcpo_p3[i] *
                            sumdajkddjkdyj + lcpo_p4[i] * lastyj);
        daidzj =  (lcpo_p2[i] * daijddijdzj + lcpo_p3[i] *
                            sumdajkddjkdzj + lcpo_p4[i] * lastzj);

        // if(fabs(daidxj)>1)daidxj=0;
        // if(fabs(daidyj)>1)daidyj=0;
        // if(fabs(daidzj)>1)daidzj=0;

        fx[j] = fx[j] - tran_ene_guan [j] * daidxj * frespa;
        fy[j] = fy[j] - tran_ene_guan [j] * daidyj * frespa;
        fz[j] = fz[j] - tran_ene_guan [j] * daidzj * frespa;
      }

      ai = lcpo_p1[i] * surf_i + lcpo_p2[i] * sumaij + lcpo_p3[i] * sumajk +
           lcpo_p4[i] * sumaijajk;

      daidxi = (lcpo_p2[i] * sumdaijddijdxi + lcpo_p4[i] *
                          sumdaijddijdxiajk);
      daidyi = (lcpo_p2[i] * sumdaijddijdyi + lcpo_p4[i] *
                          sumdaijddijdyiajk);
      daidzi = (lcpo_p2[i] * sumdaijddijdzi + lcpo_p4[i] *
                          sumdaijddijdziajk);

      // if(fabs(daidxi)>2) daidxi=0;
      // if(fabs(daidyi)>2) daidyi=0;
      // if(fabs(daidzi)>2) daidzi=0;

      fx[i] = fx[i] - tran_ene_guan [i] * daidxi*frespa;
      fy[i] = fy[i] - tran_ene_guan [i] * daidyi*frespa;
      fz[i] = fz[i] - tran_ene_guan [i] * daidzi*frespa;
// if(klok % 1000 == 0)
// printf("%d    %lf   %lf    %lf    %lf\n", i,tran_ene_guan [i],daidxi,daidyi,daidzi);
// printf("%lf    %lf    %lf\n",fx[i]-foldx[i],fy[i]-foldy[i],fz[i]-foldz[i]);
// printf("%lf    %lf    %lf\n",x[i]-oldx[i],y[i]-oldy[i],z[i]-oldz[i]);
      // totsasa = totsasa + tran_ene_guan [i] * ai;

    }
  }
  // esurf = surften * totsasa;

}
