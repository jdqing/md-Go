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
	init_Go();//initialize the program.

	estimate_Go();//estimate the parameters before running.

	run_Go();//begin run the MD.

	con_Go();//continue if broken by some reason.

	deal_for("RMSD");//deal the data for some results after running. e.g. Cv, RMSD...

	printf("%10.8lf\n",pi );
	return 0;
}
