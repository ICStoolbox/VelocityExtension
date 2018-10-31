#include "velext.h"
#include "vl_calls.h"


VLst *VL_init(int dim,int ver,char mfree) {
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
  vlst->sol.nbcl = 0;
  vlst->sol.nmat = 0;

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


/* set params (facultative):  verb= '-|0|+',  zip = 0 / 1 */
void VL_setPar(VLst *vlst,char verb,int zip) {
  vlst->info.verb = verb;
  vlst->info.zip  = zip;
}


/* handle boundary conditions:
  typ= Dirichlet, Load
  ref= integer
  att= char 'v', 'f', 'n'
  elt= enum VL_ver, VL_edg, VL_tri, VL_tet */
int VL_setBC(VLst *vlst,int typ,int ref,char att,int elt,double *u) {
  Cl    *pcl;
  int    i;

  pcl = &vlst->sol.cl[vlst->sol.nbcl];
  pcl->typ = typ;
  pcl->ref = ref;
  pcl->att = att;
  pcl->elt = elt;

  if ( pcl->typ == Dirichlet ) {
    if ( !strchr("fv",pcl->att) ) {
      fprintf(stdout,"\n # wrong format: %c\n",pcl->att);
      return(0);
    }
  }
  else if ( pcl->typ == Neumann ) {
    if ( !strchr("fnv",pcl->att) ) {
      if ( vlst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: %c\n",pcl->att);
      return(0);
    }
    if ( (pcl->elt == VL_ver) && (pcl->att == 'n') ) {
      if ( vlst->info.verb != '0' )  fprintf(stdout,"\n # condition not allowed: %c\n",pcl->att);
      return(0);
    }
  }

  if ( pcl->att == 'v' ) {
    for (i=0; i<vlst->info.dim; i++)  pcl->u[i] = u[i];
  }
  else if ( pcl->att == 'n' ) {
    pcl->u[0] = u[0];
  }

  if ( vlst->sol.nbcl == VL_CL-1 )  return(0);
  vlst->sol.nbcl++;

  return(1);
}


/* specify elasticity Lame coefficients */
int VL_setLame(VLst *vlst,int ref,double alpha) {
  Mat    *pm;

  if ( vlst->sol.nmat == VL_MAT-1 )  return(0);

  pm = &vlst->sol.mat[vlst->sol.nmat];
  pm->ref    = ref;
  pm->alpha = alpha;

  vlst->sol.nmat++;

  return(1);
}


/* Construct solution */
int VL_newSol(VLst *vlst) {

  vlst->sol.u = (double*)calloc(vlst->info.dim*vlst->info.np,sizeof(double));
  assert(vlst->sol.u);
  return(1);
}


/* Add element u[dim] to solution at ip */
int VL_addSol(VLst *vlst,int ip,double *u) {
  memcpy(&vlst->sol.u[vlst->info.dim*(ip-1)],u,vlst->info.dim*sizeof(double));
  return(1);
}


/* construct mesh */
int VL_mesh(VLst *vlst,int np,int na,int nt,int ne) {
  if ( !vlst )  return(0);

  vlst->info.np = np;
  vlst->info.na = na;
  vlst->info.nt = nt;
  vlst->info.ne = ne;

  /* bound on number of nodes */
  vlst->mesh.point = (pPoint)calloc(vlst->info.np+1,sizeof(Point));
  assert(vlst->mesh.point);

  if ( vlst->info.na ) {
    vlst->mesh.edge  = (pEdge)calloc(vlst->info.na+1,sizeof(Edge));
    assert(vlst->mesh.edge);
  }
  if ( vlst->info.nt ) {
    vlst->mesh.tria  = (pTria)calloc(vlst->info.nt+1,sizeof(Tria));
    assert(vlst->mesh.tria);
  }
  if ( vlst->info.ne ) {
    vlst->mesh.tetra  = (pTetra)calloc(vlst->info.ne+1,sizeof(Tetra));
    assert(vlst->mesh.tetra);
  }

  return(1);
}


/* insert mesh elements into structure */
int VL_addVer(VLst *vlst,int idx,double *c,int ref) {
  pPoint   ppt;
  int      i;

  assert(idx > 0 && idx <= vlst->info.np);
  ppt = &vlst->mesh.point[idx];
  for (i=0; i<vlst->info.dim; i++)
    ppt->c[i] = c[i];
  ppt->ref = ref;

  return(1);
}


int VL_addEdg(VLst *vlst,int idx,int *v,int ref) {
  pEdge   pe;

  assert(idx > 0 && idx <= vlst->info.na);
  pe = &vlst->mesh.edge[idx];
  memcpy(&pe->v[0],&v[0],2*sizeof(int));
  pe->ref = ref;

  return(1);
}


int VL_addTri(VLst *vlst,int idx,int *v,int ref) {
  pTria   pt;

  assert(idx > 0 && idx <= vlst->info.nt);
  pt = &vlst->mesh.tria[idx];
  memcpy(&pt->v[0],&v[0],3*sizeof(int));
  pt->ref = ref;

  return(1);
}


int VL_addTet(VLst *vlst,int idx,int *v,int ref) {
  pTetra   pt;

  assert(idx > 0 && idx <= vlst->info.ne);
  pt = &vlst->mesh.tetra[idx];
  memcpy(&pt->v[0],&v[0],4*sizeof(int));
  pt->ref = ref;

  return(1);
}


/* initialize solution vector or Dirichlet conditions
   return: 1 if completion
           0 if no vertex array allocated
          -1 if previous data stored in struct. */
int VL_iniSol(VLst *vlst,double *u) {
  if ( !vlst->info.np )  return(0);

  /* no data already allocated */
  if ( !vlst->sol.u ) {
    vlst->sol.u  = (double*)u;
    return(1);
  }
  /* resolve potential conflict */
  else {
    free(vlst->sol.u);
    vlst->sol.u  = (double*)u;
    return(-1);
  }
}


/* return pointer to solution (Warning: starts at address 0) */
double *VL_getSol(VLst *vlst) {
  return(vlst->sol.u);
}


int VL_velext(VLst *vlst) {
  Cl   *pcl;
  int   i,ier;

  for (i=0; i<vlst->sol.nbcl; i++) {
    pcl = &vlst->sol.cl[i];
    vlst->sol.cltyp |= pcl->typ;
    vlst->sol.clelt |= pcl->elt;
  }

  if ( vlst->info.dim == 2)
    ier = velex1_2d(vlst);
  else
    ier = velex1_3d(vlst);
  if ( ier < 1 )  return(ier);

//if ( !sctove(&mesh,&sol))   return(1);

  return(ier);
}

