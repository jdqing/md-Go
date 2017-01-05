/*
 * This header includes global constants and variables for this project.
 */
#ifndef _GLBCOVA_H
#define _GLBCOVA_H
#include "../includefi.h"

#define NPI 1 //the number of different interactions which require a list
				   //LJ (WCA) -> InInt = 0; Coulomb -> InInt = 1

extern const double pi ;
extern long int klok;
extern double Q;  //add by Meng

//Box parameters
extern double side_X, side_Y, side_Z, side_XYZ;
extern double side_hX, side_hY, side_hZ;
extern double BORDER_MIN;
extern double BORDER_MAX;

//General simulation parameters
extern double h;
extern double hsq2;
extern double T;//temperature

//Particle numbers
extern int TER0 ;
extern int TER1, TER2, TER3, TER4;
extern int N_amino ;  //the number of residues in the protein
extern int NBP_atom ; //the number of backbone atoms in the protein
                     //(QM: only heavy atoms and hydrogen atoms excluded)
extern int NB_atom ; //the total number of backbone atoms
                    //in the protein and crowders
extern int * NS_atom ; //the number of sites per backbone atom or crowder
                      //("NS_atom = m_crwd - 1" for crowders)
extern int NTP_atom ; //the total number of atoms in the protein
extern int NT_atom ; //the total number of atoms in the protein and crowders
extern int ** atom_key ;
extern char * amino_key ;
extern int * crowder_key ;
extern int * part_key ;
extern int * maxi_key ;
extern int * INDX , * JNDX ;
extern int * residue_mass ;
extern int ** residue_content ;
extern int * crowder_mass ;
extern int ** crowder_content ;
extern int * RIGID_SET , * RIGID_END ;

//Individual properties of atoms
extern const int ns ; //the number of different atom types
extern const int nps ; //the number of different pairs of atoms
extern const double MASS [] ; //atom masses: N, C, O, H, S, 1UBQ 8564.47
                             //(QM: mass of ubiquitin)
extern const double D [] ; //largest atom diameters: N, C, O, H, S, 1UBQ
extern double PD [];
extern double PD_sq [];
extern const double VISC [] ; //0.1 water viscosity
extern double K1 [], K2 [];
extern double SIGMA_FORCE []; //thermostat random force

////////////////////////////////////////////////////////////////////////////////
extern const int na ; //the number of different amber atom types
extern const int npa ; //the number of different pairs of amber atoms
extern double RLJ [];
extern double ELJ [];
extern double QC [];
extern double D_LJ []; //LJ pair diameters
extern double D2_LJ []; //LJ square pair diameters
extern double D2_LJ_CUTOFF []; //LJ square cutoff distances
extern double E_LJ []; //LJ energy prefactors
extern double F_LJ []; //LJ force prefactors

//Crowders (QM: for crowders)
extern const int n_crwd ; //the number of different crowder types
extern double cM_crwd [] ; //molar concentrations of crowders
extern double phi_crwd [] ; //volume fractions of crowders
extern double V_crwd []; //crowder volumes
extern int N_crwd []; //the number of crowders of each type
extern int m_crwd []; //the number of sites in a crowder of each type
extern double * RX_crwd []; //x coordinates of sites in a crowder of each type
extern double * RY_crwd []; //y coordinates of sites in a crowder of each type
extern double * RZ_crwd []; //z coordinates of sites in a crowder of each type
extern int * part_key_crwd []; //site (atom) "part_key" in a crowder of each type
extern int * maxi_key_crwd []; //site (atom) "maxi_key" in a crowder of each type
extern int RIGID_SET_crwd []; //the first site in a crowder rigid sub-unit
extern int RIGID_END_crwd []; //the last site in a crowder rigid sub-unit

//Coordinates, velocities and forces
////////////////////////////////////////////////////////////////////////////////
extern double * x ;
extern double * y ;
extern double * z ;

extern double * x_old ;
extern double * y_old ;
extern double * z_old ;

extern double * vx ;
extern double * vy ;
extern double * vz ;

extern double * fx ;
extern double * fy ;
extern double * fz ;

//Centre of mass motion
extern double * CMS ;
extern double * VCMS ;
extern double * RMASS ;
extern double * RVISC ;
extern double * RXYZ ;

//Rotational motion
extern double * AXES ;
extern double * W ;
extern double ** RX , ** RY , ** RZ ;
extern double * IR1 , * IR2 , * IR3 ;
extern double * AV , * BV , * CV ;
extern double * FV , * GV , * HV ;

////////////////////////////////////////////////////////////////////////////////
//Crowder constraints
extern const double r_a260_a260 ; extern const double k_a260_a260;
extern const double theta_a260_a260_a260 ; extern const double k_a260_a260_a260 ;

//Water constraints
extern const double r_OW_HW ; extern const double k_OW_HW ;
extern const double theta_HW_OW_HW ; extern const double k_HW_OW_HW ;

//Polymer constraints
extern const double r_C_N ; extern const double k_C_N ;
extern const double r_C_O ; extern const double k_C_O ;

////////////////////////////////////////////////////////////////////////////////
extern const double theta_C_N_CT ; extern const double k_C_N_CT ;
extern const double theta_CT_C_N ; extern const double k_CT_C_N ;
extern const double theta_CT_C_O ; extern const double k_CT_C_O ;
extern const double theta_N_C_O ; extern const double k_N_C_O ;

////////////////////////////////////////////////////////////////////////////////
/*
//Crowder constraints
extern const double r_a260_a260 ; extern const double k_a260_a260 ;
extern const double theta_a260_a260_a260 ; extern const double k_a260_a260_a260;

//Water constraints
extern const double r_OW_HW ; extern const double k_OW_HW ;
extern const double theta_HW_OW_HW ; extern const double k_HW_OW_HW ;

//Polymer constraints
////////////////////////////////////////////////////////////////////////////////
extern const double theta_C_N_CT ; extern const double k_C_N_CT = 50.0;
extern const double theta_CT_C_N ; extern const double k_CT_C_N = 70.0;
extern const double theta_CT_C_O ; extern const double k_CT_C_O = 80.0;

//const double theta_N_C_O ; const double k_N_C_O = 80.0;
extern const double theta_O_C_N ; extern const double k_O_C_N = 80.0;
extern const double theta_N_C_O ; extern const double k_N_C_O = 80.0;

//----------------------------------
extern const double theta_N_CT_C ; extern const double k_N_CT_C ;
extern const double theta_CT_N_CT ; extern const double k_CT_N_CT ;
extern const double theta_N_CT_CT ; extern const double k_N_CT_CT ;
extern const double theta_CT_CT_C ;extern const double k_CT_CT_C ;
extern const double theta_CT_CT_S ; extern const double k_CT_CT_S ;
extern const double theta_CT_CT_O ; extern const double k_CT_CT_O ;

extern const double r_C_N ;   extern const double k_C_N ;  // why not 490.0 ?
extern const double r_C_O ;   extern const double k_C_O ;  //why not 570.0
extern const double r_N_CT ;  extern const double k_N_CT ; //337.0;
extern const double r_CT_C ;  extern const double k_CT_C ; //317.
extern const double r_CT_CT ; extern const double k_CT_CT ; //310.
extern const double r_CT_O ; extern const double k_CT_O ; //320.
extern const double r_CT_S ; extern const double k_CT_S ; //320.

*/
/////////////Cells (for list population or pcfs calculation)////////////////////
extern int McX , McY , McZ , McT ; //cell numbers
extern int ** cell_content , * cell_mass ; //cell content and mass

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Lists//////////////////////////////////////////
extern int * list_content [][NPI];//[nps][npi]  ->  [][NPI]  multidimensional
                    //array must have bounds for all dimensions except the first
extern int  list_mass [][NPI];
extern double  R2_CUTOFF [][NPI];
extern double  R1_LIST [][NPI];
extern double  R2_LIST [][NPI];
extern double  PROGRESS [][NPI];
extern double  dr1 [][NPI];

////////////////////////////////////////////////////////////////////////////////
extern int * bond_list_content ;
extern int bond_list_mass ;
extern int * valence_list_content ;
extern int valence_list_mass ;
extern int * dihedral_list_content ;
extern int dihedral_list_mass ;

////////////////////////////////////////////////////////////////////////////////
extern int * native_list_content ;
extern int native_list_mass ;
extern double * D2_LJ_n ;
extern double * D2_LJ_CUTOFF_n ;
extern double * E_LJ_n ;
extern double * F_LJ_n ;
extern double LJ_SCALE_n;

////////////////////////////////////////////////////////////////////////////////
extern double k_RG;
extern double RG_0;
extern double k_NATIVE;
extern double Q0_NATIVE;
extern int N0_NATIVE;
extern double * THETA_NATIVE ;
extern int COUNT_NATIVE ;
extern double Q_NATIVE [101];

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////Observables///////////////////////////////////////
extern long int COUNT_Bfactor ;
extern double * RB ;
extern double * RB2 ;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////Cavity water//////////////////////////////////////
extern double xHOH [5] ;
extern double yHOH [5] ;
extern double zHOH [5] ;

//////////////added by jdq//////////////////////////////////////////////////////
//for rmsd
extern double * x0 ;
extern double * yy0 ;
extern double * z0 ;

extern double * x1 ;
extern double * yy1 ;
extern double * z1 ;
//
extern char **Flags2Atom;
//'b'->bonds;'a'->angles;'d'->dihedral;'c'->native contacts;
//'s'->same;'o'->others

#endif
