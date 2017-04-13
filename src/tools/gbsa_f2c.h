#ifndef _GBSA_F2C_H
#define _GBSA_F2C_H


#include "../init/init.h"

void MTM_init(int atm_cnt,const double Concent);

void MTM_ene( int atm_cnt);

void MTM_lcpo(int max_count, int max_i, double frespa);

#endif
