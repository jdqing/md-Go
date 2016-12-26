/*This is a tool for caculte the sasa of a protein using LCPO method.
 *
 */

double *gbsafrc;
double *vdwrad;
double *lcpo_p1;
double *lcpo_p2;
double *lcpo_p3;
double *lcpo_p4;

int *ineighbor;
int *j_vals;
int *k_vals;

int sasa_lcpo(double *crd, double *frc, int max_count,
              int max_i, double frespa, double *esurf)
{
  int    nei_count;
  int    i, j, k;
  int    num_j_vals;
  int    num_k_vals;
  int    icount;
  int    count2;
  int    count2_fin;
  int    j_loop, k_loop;

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

  //Loop over i
  for(i = 1; i <= max_i; i++)
  {
    if(ineighbor[nei_count] == 0) nei_count = nei_count + 1;
    else
    {
      surf_i = 4.0 * PI * vdwrad[i] * vdwrad[i];
      sumaij            = 0;
      sumajk            = 0;
      sumaijajk         = 0;
      sumdaijddijdxi    = 0;
      sumdaijddijdyi    = 0;
      sumdaijddijdzi    = 0;
      sumdaijddijdxiajk = 0;
      sumdaijddijdyiajk = 0;
      sumdaijddijdziajk = 0;

      icount = nei_count;
      num_j_vals = 0;

      do{
        j = ineighbor[nei_count];
        num_j_vals = num_j_vals + 1;
        j_vals[num_j_vals] = j;
        nei_count = nei_count + 1;

        if(ineighbor[nei_count] != 0) continue;

        nei_count = nei_count + 1;
        break;
      }while(1);

      xi = crd[3 * i - 2];
      yi = crd[3 * i - 1];
      zi = crd[3 * i    ];

      //Loop over j
      for(j_loop = 1; j_loop <= num_reals; j_loop++)
      {
        j = j_vals[j_loop];

        xj = crd[3 * j - 2];
        yj = crd[3 * j - 1];
        zj = crd[3 * j    ];

        xji = xj - xi;
        yji = yj - yi;
        zji = zj - zi;
        rij = sqrt(xji * xji + yji * yji + zji * zji);
        one_rij = 1 / rij;

        //Equation 3 of Weiser, et. al.; J Comput Chem, v20 p217

        aij = 2.0 * PI * vdwrad[i] * ( vdwrad[i] - rij * 0.5 -
              (vdwrad[i] * vdwrad[i] - vdwrad[j] * vdwrad[j]) * one_rij * 0.5);

        //printf();

        //First derivatives: Appendix A of Weiser, et. al.

        daijddij = PI * vdwrad[i] * (one_rij * one_rij *
                   (vdwrad[i] * vdwrad[i] - vdwrad[j] * vdwrad[j] - 1);
        daijddijdxj = daijddij * xji * one_rij;
        daijddijdyj = daijddij * yji * one_rij;
        daijddijdzj = daijddij * zji * one_rij;

        //Add it up

        sumaij = sumaij + aij;

        count2 = icount;

        sumajk2 = 0.0;
        sumdajkddjkdxj = 0.0;
        sumdajkddjkdyj = 0.0;
        sumdajkddjkdzj = 0.0;

        p3p4aij = -surften * (lcpo_p3[i] + lcpo_p4[i] * aij) * frespa;

        //Loop over k
        for(k_loop = count2; k_loop <= max_count; k_loop++)
        {
          if(ineighbor[k_loop] == 0)
          {
            count2_fin = k_loop - 1;
            break;
          }
        }

        num_k_vals = 0;

        for(k_loop = count2; k_loop <= count2_fin; k_loop++)
        {
          k = ineighbor[k_loop];

          if(k == j)
          {
            num_k_vals = num_k_vals + 1;
            k_vals[num_k_vals] = k;
          }
        }

        count2 = icount;

        for(k_loop = 1;k_loop <= num_k_vals; k_loop++)
        {
          k  = k_vals[k_loop];
          xk = crd[3 * k - 2];
          yk = crd[3 * k - 1];
          zk = crd[3 * k    ];

          xjk = xj - xk;
          yjk = yj - yk;
          zjk = zj - zk;
          rjk = sqrt(xjk * xjk + yjk * yjk + zjk * zjk);
          one_rjk = 1.0 / rjk;


        }




      }

    }

  }

}
