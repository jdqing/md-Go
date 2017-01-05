/*
 * LCPO
 *
 * This is a tool for caculte the sasa of a protein using LCPO method.
 *This code was used to try converting the gbsa.F90 to C Code.
 *But I want to try link the F90 file to the C project first.
 *or you can say mixed-language compiling.
 *
 * Not a complete code!!!
 */

double **gbsafrc;
double *vdwrad;
double *lcpo_p1;
double *lcpo_p2;
double *lcpo_p3;
double *lcpo_p4;

int **ineighbor;
int *j_vals;
int *k_vals;

int sasa_setup(int atm_cnt, int *num_ints, int *num_reals)
{
  int alloc_failed;
  int *bond_partner;
  int *total_bond_partner;

  int i;
  int atm_1;
  int atm_2;
  int num_bonds;
  int total_bonds;
  int atomicnumber;
  char isymbl[4]; //holder for atom type

  //num_ints and num_reals are used to return allocation counts. Don't zero.

  //ROMEE Add this before calling  if (gbsa == 1) then

  //First we calloc the main arrays

  gbsafrc = (double **)calloc((atm_cnt+1) * sizeof(double *));
  for(i=1;i<=atm_cnt;i++)
  {
    *(gbsafrc + i) = (double *)calloc(3 * sizeof(double));
  }

  vdwrad  = (double *)calloc((atm_cnt+1) * sizeof(double));
  lcpo_p1 = (double *)calloc((atm_cnt+1) * sizeof(double));
  lcpo_p2 = (double *)calloc((atm_cnt+1) * sizeof(double));
  lcpo_p3 = (double *)calloc((atm_cnt+1) * sizeof(double));
  lcpo_p4 = (double *)calloc((atm_cnt+1) * sizeof(double));

  ineighbor = (int **)calloc((atm_cnt+1) * sizeof(int *));
  for(i=1;i<=atm_cnt;i++)
  {
    *(ineighbor + i) = (int *)calloc(40 * sizeof(int));
  }

  j_vals = (int *)calloc((atm_cnt+1) * sizeof(int));
  k_vals = (int *)calloc((atm_cnt+1) * sizeof(int));
  bond_partner = (int *)calloc((atm_cnt+1) * sizeof(int));
  total_bond_partner = (int *)calloc((atm_cnt+1) * sizeof(int));

  /*
  num_reals = num_reals + size(lcpo_p1) + &
                          size(lcpo_p2) + &
                          size(lcpo_p3) + &
                          size(lcpo_p4) + &
                          size(vdwrad)*/
  num_reals = num_reals + 5 * atm_cnt;

  num_ints = num_ints + (1 + 1 + 40 + 1 + 1) * atm_cnt;

  //For each atom index we fing in the BONDS_WITHOUT_HYDROGEN array, add
  //1 to that index of bond_partner.taht way,we can tell how many atoms
  //are bonded to each atom

  for(i = 0; i <= nbona -1; i++)
  {
    atm_1 = gbl_bond[bonda_idx + i]%atm_i;
    atm_2 = gbl_bond[bonda_idx + i]%atm_j;
    bond_partner[atm_1] = bond_partner[atm_1] + 1;
    bond_partner[atm_2] = bond_partner[atm_2] + 1;
  }

  total_bond_partner = bond_partner;

  for(i = 1; i <= nbonh; i++)
  {
    atm_1 = gbl_bond[i]%atm_i;
    atm_2 = gbl_bond[i]%atm_j;
    total_bond_partner[atm_1] = total_bond_partner[atm_1] + 1;
    total_bond_partner[atm_2] = total_bond_partner[atm_2] + 1;
  }

  //Construct the parameters for the SA calculation; note taht the radii
  //stored in vdwrad are augmented by 1.4 Angstroms

  for(i = 1; i <= atm_cnt; i++)
  {
    if(loaded_atm_atomicnumber)
      atomicnumber = atm_atomicnumber[i];
    else
      call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber);

    isymbl = atm_isymbl[i];
    call upper(isymbl);
    num_bonds = bond_partner[i];
    total_bond_partner = total_bond_partner[i];

    if(atomicnumber == 6)
    {
      if(total_bonds == 4)
      {
        if(num_bonds == 1)
        {
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.77887d0
          lcpo_p2(i) = -0.28063d0
          lcpo_p3(i) = -0.0012968d0
          lcpo_p4(i) = 0.00039328d0
        }
        else if(num_bonds == 2)
        {
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.56482d0
          lcpo_p2(i) = -0.19608d0
          lcpo_p3(i) = -0.0010219d0
          lcpo_p4(i) = 0.0002658d0
        }
        else if(num_bonds == 3)
        {
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.23348d0
          lcpo_p2(i) = -0.072627d0
          lcpo_p3(i) = -0.00020079d0
          lcpo_p4(i) = 0.00007967d0
        }
        else if(num_bonds == 4)
        {
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.00000d0
          lcpo_p2(i) = 0.00000d0
          lcpo_p3(i) = 0.00000d0
          lcpo_p4(i) = 0.00000d0
        }
        else
        {
          write(6,*) 'Unusual nbond for sp3 C:', i, num_bonds, &
            ' Using default carbon LCPO parameters'
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.77887d0
          lcpo_p2(i) = -0.28063d0
          lcpo_p3(i) = -0.0012968d0
          lcpo_p4(i) = 0.00039328d0
        }
      }
      else
      {
        if(num_bonds == 2)
        {
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.51245d0
          lcpo_p2(i) = -0.15966d0
          lcpo_p3(i) = -0.00019781d0
          lcpo_p4(i) = 0.00016392d0
        }
        else if(num_bonds == 3)
        {
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.070344d0
          lcpo_p2(i) = -0.019015d0
          lcpo_p3(i) = -0.000022009d0
          lcpo_p4(i) = 0.000016875d0
        }
        else
        {
          write(6,*) 'Unusual nbond for C :', i, num_bonds, &
            ' Using default carbon LCPO parameters'
          vdwrad(i)  = 1.70d0 + 1.4d0
          lcpo_p1(i) = 0.77887d0
          lcpo_p2(i) = -0.28063d0
          lcpo_p3(i) = -0.0012968d0
          lcpo_p4(i) = 0.00039328d0
        }
      }
    }
    else if(atomicnumber == 8)
    {
      if(isymbl[]=="O")
      {
        vdwrad(i)  = 1.60d0 + 1.4d0
        lcpo_p1(i) = 0.68563d0
        lcpo_p2(i) = -0.1868d0
        lcpo_p3(i) = -0.00135573d0
        lcpo_p4(i) = 0.00023743d0
      }
      else if(isymbl[]=="O2")
      {
        vdwrad(i)  = 1.60d0 + 1.4d0
        lcpo_p1(i) = 0.88857d0
        lcpo_p2(i) = -0.33421d0
        lcpo_p3(i) = -0.0018683d0
        lcpo_p4(i) = 0.00049372d0
      }
      else
      {
        if(num_bonds == 1)
        {
          vdwrad(i)  = 1.60d0 + 1.4d0
          lcpo_p1(i) = 0.77914d0
          lcpo_p2(i) = -0.25262d0
          lcpo_p3(i) = -0.0016056d0
          lcpo_p4(i) = 0.00035071d0
        }
        else if(num_bonds == 2)
        {
          vdwrad(i)  = 1.60d0 + 1.4d0
          lcpo_p1(i) = 0.49392d0
          lcpo_p2(i) = -0.16038d0
          lcpo_p3(i) = -0.00015512d0
          lcpo_p4(i) = 0.00016453d0
        }
        else
        {
          write(6,*) 'Unusual nbond for O:', i, num_bonds, &
             ' Using default oxygen LCPO parameters'
          vdwrad(i)  = 1.60d0 + 1.4d0
          lcpo_p1(i) = 0.77914d0
          lcpo_p2(i) = -0.25262d0
          lcpo_p3(i) = -0.0016056d0
          lcpo_p4(i) = 0.00035071d0
        }
      }
    }
    else if(atomicnumber == 7)
    {
      if(isymbl == "N3")
      {
        if(num_bonds == 1)
        {
          vdwrad(i)  = 1.65d0 + 1.4d0
          lcpo_p1(i) = 0.078602d0
          lcpo_p2(i) = -0.29198d0
          lcpo_p3(i) = -0.0006537d0
          lcpo_p4(i) = 0.00036247d0
        }
        else if(num_bonds == 2)
        {
          vdwrad(i)  = 1.65d0 + 1.4d0
          lcpo_p1(i) = 0.22599d0
          lcpo_p2(i) = -0.036648d0
          lcpo_p3(i) = -0.0012297d0
          lcpo_p4(i) = 0.000080038d0
        }
        else if(num_bonds == 3)
      }
    }
  }


}

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
