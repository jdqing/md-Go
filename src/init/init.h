/*
 * This is the init part and ending part of the project.
 */
#ifndef _INIT_H
#define _INIT_H

#include "glbcova.h"
#include "../box/boxoperation.h"
#include "../fun/force.h"
#include "../list/listroutine.h"
#include "../inout/readin.h"
#include "../tools/progress-bar.h"


int init_Go(const char *cfgfile);

int end_Go();

#endif
