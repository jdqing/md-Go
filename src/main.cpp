/*
 * The main.c of the molecular dynamics simulation using all-atom Go model.
 */
#include "includefi.h"
#include "init/glbcova.h"
#include "init/init.h"

#include "estimate/estimateGo.h"
#include "run/runGo.h"
#include "datadeal/dealfor.h"


int main(int argc, char *argv[])
{
	printf("////////////////////////////Go_Begin/////////////////////////////\n");
	printf("////Compiling Date: %s\n",__DATE__);
	printf("////Compiling Time: %s\n",__TIME__);
	printf("////////MTM_yes(1)/no(0): %d\n",MTM_yn);
	printf("/////////////////////////////////////////////////////////////////\n");

	char cfgfile[50];//the config file name.
	if(argc == 2){
		strcpy(cfgfile,argv[1]);
		printf("////Select the config file : %s\n",cfgfile);
	}
	else{
		printf("////Just accept 1 parameter. e.g. \"./a.out gomodel.cfg\"\n");
		return 0;
	}

	if(!init_Go(cfgfile))//initialize the program.
	{
		printf("//////ERROR: init_Go not succeed\n");
		return 0;
	}

	// estimate_Go();//estimate the parameters before running.

	// printf("LJ_SCALE_n=%lf\n",LJ_SCALE_n );

	run_Go();//begin run the MD.

	// con_Go();//continue if broken by some reason.

	// deal_for("RMSD");//deal the data for some results after running. e.g. Cv, RMSD...

	end_Go();

	return 0;
}
