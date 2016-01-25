#include "velext.h"
#include "vl_calls.h"


VLst *VL_init(int dim,int ver,char typ,char mfree) {
	VLst   *vlst;

  /* default values */
  vlst = (VLst*)calloc(1,sizeof(VLst));
  memset(&vlst->mesh,0,sizeof(Mesh));

  /* solution structure */
  memset(&vlst->sol,0,sizeof(Sol));
  vlst->sol.cl    = (Cl*)calloc(VL_CL,sizeof(Cl));
  vlst->sol.mat   = (Mat*)calloc(VL_MAT,sizeof(Mat));
  vlst->sol.res   = VL_RES;
  vlst->sol.nit   = VL_MAXIT;

  /* global parameters */
  vlst->info.dim    = dim;
  vlst->info.ver    = ver;
  vlst->info.verb   = '1';
  vlst->info.zip    = 0;
  vlst->info.ls     = 0;
  vlst->info.mfree  = mfree;

  /* init timer */
  tminit(vlst->info.ctim,TIMEMAX);
  chrono(ON,&vlst->info.ctim[0]);

  return(vlst);
}

/* free global data structure */
int VL_stop(VLst *vlst) {
	char   stim[32];

	/* release memory */
  free(vlst->sol.u);
	free(vlst->sol.cl);
	free(vlst->sol.mat);

  chrono(OFF,&vlst->info.ctim[0]);
  if ( vlst->info.verb != '0' ) {
    printim(vlst->info.ctim[0].gdif,stim);
    fprintf(stdout,"\n  Cumulative time: %s sec.\n",stim);
  }

	return(1);
}


int VL_velext(VLst *vlst) {
  int   ier;

  if ( vlst->info.dim == 2)
		ier = velex1_2d(vlst);
	else
		ier = velex1_3d(vlst);
  if ( ier < 1 )  return(ier);

//if ( !sctove(&mesh,&sol))   return(1);

	return(ier);	
}
