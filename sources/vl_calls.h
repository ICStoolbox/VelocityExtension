#ifndef __VE_CALLS_H
#define __VE_CALLS_H


enum {Dirichlet=1, Neumann=2, Gravity=4};
enum {VL_ver=1,VL_edg=2,VL_tri=4,VL_tet=8};


/* data structure */
typedef struct _VLst VLst;


/* prototypes */
VLst   *VL_init(int dim, int ver,char mfree);
int     VL_stop(VLst *vlst);
void    VL_setPar(VLst *vlst,char verb,int zip);
int     VL_setBC(VLst *vlst,int typ,int ref,char att,int elt,double *u);
int     VL_setLame(VLst *vlst,int ref,double alpha);
int     VL_newSol(VLst *vlst);
int     VL_addSol(VLst *vlst,int ip,double *u);
int     VL_mesh(VLst *vlst,int np,int na,int nt,int ne);
int     VL_addVer(VLst *vlst,int idx,double *c,int ref);
int     VL_addEdg(VLst *vlst,int idx,int *v,int ref);
int     VL_addTri(VLst *vlst,int idx,int *v,int ref);
int     VL_addTet(VLst *vlst,int idx,int *v,int ref);
int     VL_iniSol(VLst *vlst,double *u);
double *VL_getSol(VLst *vlst);
int     VL_velext(VLst *vlst);

//int   LS_setCoef(VLst *vlst,int ref,double alpha); // does not appear anymore in any *.c file

#endif

