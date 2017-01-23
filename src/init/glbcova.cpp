/*
 * the initialization of global constants and variables.
 */
#include "glbcova.h"

const double pi = 3.14159265358979323846264338328;

long int klok;

double Q;  //add by Meng

//Box parameters
double side_X, side_Y, side_Z, side_XYZ;
double side_hX, side_hY, side_hZ;
double BORDER_MIN;
double BORDER_MAX;

//General simulation parameters
double h = 0.02;
double hsq2 = 0.5 * h * h;
double T;//temperature

//Particle numbers
int TER0 = 0;
int TER1, TER2, TER3, TER4;
int N_amino = 0;  //the number of residues in the protein
int NBP_atom = 0; //the number of backbone atoms in the protein
                 //(QM: only heavy atoms and hydrogen atoms excluded)
int NB_atom = 0; //the total number of backbone atoms in the protein and crowders
int * NS_atom = NULL; //the number of sites per backbone atom or crowder
                     //("NS_atom = m_crwd - 1" for crowders)
int NTP_atom = 0; //the total number of atoms in the protein
int NT_atom = 0; //the total number of atoms in the protein and crowders
int ** atom_key = NULL;
char * amino_key = NULL;
int * crowder_key = NULL;
int * part_key = NULL;
int * maxi_key = NULL;
int * INDX = NULL, * JNDX = NULL;
int * residue_mass = NULL;
int ** residue_content = NULL;
int * crowder_mass = NULL;
int ** crowder_content = NULL;
int * RIGID_SET = NULL, * RIGID_END = NULL;

//Individual properties of atoms
const int ns = 6; //the number of different atom types
const int nps = ns * (ns + 1) / 2; //the number of different pairs of atoms
const int npi = NPI;
const double MASS [ns] = { 14.01, 12.01, 16.00, 1.008, 32.06, 8564.47 };
              //atom masses: N, C, O, H, S, 1UBQ 8564.47 (QM: mass of ubiquitin)
const double D [ns] = { 3.648, 3.816, 3.442, 2.974, 4.0, 24.0 };
              //largest atom diameters: N, C, O, H, S, 1UBQ
double PD [nps];
double PD_sq [nps];

const double VISC0 [ns] = { 27.75 * D [0], 27.75 * D [1], 27.75 * D [2],
                            27.75 * D [3], 27.75 * D [4], 27.75 * D [5] }; //water viscosity
const double VISC [ns] = {  ORDERVISC *  VISC0[0],
                            ORDERVISC *  VISC0[1],
                            ORDERVISC *  VISC0[2],
                            ORDERVISC *  VISC0[3],
                            ORDERVISC *  VISC0[4],
                            ORDERVISC *  VISC0[5]}; //default: 0.1 water viscosity
  //0.1 water viscosity   to const.cpp .h

double K1 [ns], K2 [ns];
double SIGMA_FORCE [ns]; //thermostat random force

////////////////////////////////////////////////////////////////////////////////

const int na = 260; //the number of different amber atom types
const int npa = na * (na + 1) / 2; //the number of different pairs of amber atoms
double RLJ [na + 1];
double ELJ [na + 1];
double QC [na + 1];
double D_LJ [npa + 1]; //LJ pair diameters
double D2_LJ [npa + 1]; //LJ square pair diameters
double D2_LJ_CUTOFF [npa + 1]; //LJ square cutoff distances
double E_LJ [npa + 1]; //LJ energy prefactors
double F_LJ [npa + 1]; //LJ force prefactors

//Crowders (QM: for crowders)
const int n_crwd = 4; //the number of different crowder types
double cM_crwd [ n_crwd ] = { 0.0, 0.0, 0.0, 0.0 }; //molar concentrations of crowders
double phi_crwd [ n_crwd ] = { 0.0, 0.0, 0.0, 0.0 }; //volume fractions of crowders
double V_crwd [ n_crwd ]; //crowder volumes
int N_crwd [ n_crwd ]; //the number of crowders of each type
int m_crwd [ n_crwd ]; //the number of sites in a crowder of each type
double * RX_crwd [ n_crwd ]; //x coordinates of sites in a crowder of each type
double * RY_crwd [ n_crwd ]; //y coordinates of sites in a crowder of each type
double * RZ_crwd [ n_crwd ]; //z coordinates of sites in a crowder of each type
int * part_key_crwd [ n_crwd ]; //site (atom) "part_key" in a crowder of each type
int * maxi_key_crwd [ n_crwd ]; //site (atom) "maxi_key" in a crowder of each type
int RIGID_SET_crwd [ n_crwd ]; //the first site in a crowder rigid sub-unit
int RIGID_END_crwd [ n_crwd ]; //the last site in a crowder rigid sub-unit

//Coordinates, velocities and forces

////////////////////////////////////////////////////////////////////////////////
double * x = NULL;
double * y = NULL;
double * z = NULL;

double * x_old = NULL;
double * y_old = NULL;
double * z_old = NULL;

double * vx = NULL;
double * vy = NULL;
double * vz = NULL;

double * fx = NULL;
double * fy = NULL;
double * fz = NULL;

//Centre of mass motion
double * CMS = NULL;
double * VCMS = NULL;
double * RMASS = NULL;
double * RVISC = NULL;
double * RXYZ = NULL;

//Rotational motion
double * AXES = NULL;
double * W = NULL;
double ** RX = NULL, ** RY = NULL, ** RZ = NULL;
double * IR1 = NULL, * IR2 = NULL, * IR3 = NULL;
double * AV = NULL, * BV = NULL, * CV = NULL;
double * FV = NULL, * GV = NULL, * HV = NULL;

////////////////////////////////////////////////////////////////////////////////

//Crowder constraints

const double r_a260_a260 = D [5];
const double k_a260_a260 = 310.0 * (1.526 * 1.526) / (r_a260_a260 * r_a260_a260);

const double theta_a260_a260_a260 = 120.0 / 180.0 * pi;
const double k_a260_a260_a260 = 30.0;

//Water constraints
const double r_OW_HW = 0.9572; const double k_OW_HW = 553.0;
const double theta_HW_OW_HW = 104.52 / 180.0 * pi; const double k_HW_OW_HW = 100.0;

//Polymer constraints
const double r_C_N = 1.335; const double k_C_N = 70.0;
const double r_C_O = 1.229; const double k_C_O = 70.0;

/////////////////////////////////////
const double theta_C_N_CT = 121.9 / 180.0 * pi; const double k_C_N_CT = 50.0;
const double theta_CT_C_N = 116.6 / 180.0 * pi; const double k_CT_C_N = 70.0;
const double theta_CT_C_O = 120.4 / 180.0 * pi; const double k_CT_C_O = 80.0;
const double theta_N_C_O =  122.9 / 180.0 * pi; const double k_N_C_O = 80.0;

////////////////////////////////////////////////////////////////////////////////
/*
//Crowder constraints

const double r_a260_a260 = D [5];
const double k_a260_a260 = 310.0 * (1.526 * 1.526) / (r_a260_a260 * r_a260_a260);

const double theta_a260_a260_a260 = 120.0 / 180.0 * pi;
const double k_a260_a260_a260 = 30.0;

//Water constraints
const double r_OW_HW = 0.9572; const double k_OW_HW = 553.0;
const double theta_HW_OW_HW = 104.52 / 180.0 * pi; const double k_HW_OW_HW = 100.0;

//Polymer constraints
////////////////////////////////////////////////////////////////////////////////
const double theta_C_N_CT = 121.9 / 180.0 * pi; const double k_C_N_CT = 50.0;
const double theta_CT_C_N = 116.6 / 180.0 * pi; const double k_CT_C_N = 70.0;
const double theta_CT_C_O = 120.4 / 180.0 * pi; const double k_CT_C_O = 80.0;

//const double theta_N_C_O =  122.9 / 180.0 * pi; const double k_N_C_O = 80.0;
const double theta_O_C_N =  122.9 / 180.0 * pi; const double k_O_C_N = 80.0;
const double theta_N_C_O =  122.9 / 180.0 * pi; const double k_N_C_O = 80.0;

//----------------------------------
const double theta_N_CT_C = 110.10 / 180.0 * pi; const double k_N_CT_C = 63.0;
const double theta_CT_N_CT = 118.0 / 180.0 * pi; const double k_CT_N_CT = 50.0;
const double theta_N_CT_CT = 109.7 / 180.0 * pi; const double k_N_CT_CT = 80.0;
const double theta_CT_CT_C = 111.1 / 180.0 * pi;const double k_CT_CT_C = 63.0;
const double theta_CT_CT_S = 108.6 / 180.0 * pi; const double k_CT_CT_S = 50;
const double theta_CT_CT_O = 109.5 / 180.0 * pi; const double k_CT_CT_O = 50;

const double r_C_N = 1.335;   const double k_C_N = 490.0;  // why not 490.0 ?
const double r_C_O = 1.229;   const double k_C_O = 570.0;  //why not 570.0
const double r_N_CT = 1.449;  const double k_N_CT = 337.0; //337.0;
const double r_CT_C = 1.522;  const double k_CT_C = 317.0; //317.
const double r_CT_CT = 1.526; const double k_CT_CT = 310.; //310.
const double r_CT_O = 1.410; const double k_CT_O = 320.; //320.
const double r_CT_S = 1.410; const double k_CT_S = 320.; //320.

*/
/////////////Cells (for list population or pcfs calculation)////////////////////
int McX = 0, McY = 0, McZ = 0, McT = 0; //cell numbers
int ** cell_content = NULL, * cell_mass = NULL; //cell content and mass

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Lists//////////////////////////////////////////
int * list_content [nps][NPI];
int list_mass [nps][NPI];
double R2_CUTOFF [nps][NPI];
double R1_LIST [nps][NPI];
double R2_LIST [nps][NPI];
double PROGRESS [nps][NPI];
double dr1 [nps][NPI];

////////////////////////////////////////////////////////////////////////////////
int * bond_list_content = NULL;
int bond_list_mass = 0;
int * valence_list_content = NULL;
int valence_list_mass = 0;
int * dihedral_list_content = NULL;
int dihedral_list_mass = 0;

////////////////////////////////////////////////////////////////////////////////
int * native_list_content = NULL;
int native_list_mass = 0;
double * D2_LJ_n = NULL;
double * D2_LJ_CUTOFF_n = NULL;
double * E_LJ_n = NULL;
double * F_LJ_n = NULL;
double LJ_SCALE_n;

////////////////////////////////////////////////////////////////////////////////
double k_RG;
double RG_0;
double k_NATIVE;
double Q0_NATIVE;
int N0_NATIVE;
double * THETA_NATIVE = NULL;
int COUNT_NATIVE = 0;
double Q_NATIVE [101];

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////Observables///////////////////////////////////////
long int COUNT_Bfactor = 0;
double * RB = NULL;
double * RB2 = NULL;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////Cavity water//////////////////////////////////////
double xHOH [5] = { 13.985, 15.089, 13.874, 12.628, 12.001 };
double yHOH [5] = { -15.661, -13.777, -17.351, -16.070, -20.207 };
double zHOH [5] = { 16.965, 18.567, 14.737, 13.031, 13.453 };

//////////////added by jdq//////////////////////////////////////////////////////
//for rmsd
double * x0 = NULL;
double * yy0 = NULL;
double * z0 = NULL;

double * x1 = NULL;
double * yy1 = NULL;
double * z1 = NULL;
//
char **Flags2Atom;
//'b'->bonds;'a'->angles;'d'->dihedral;'c'->native contacts;
//'s'->same;'o'->others
