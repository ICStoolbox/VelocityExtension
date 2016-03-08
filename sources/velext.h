#ifndef _VELEXT_H
#define _VELEXT_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "vl_calls.h"

#define VL_VER   "2.0a"
#define VL_REL   "Jan 19, 2016"
#define VL_CPY   "(C) Copyright 2014- , ICS-SU"

#define VL_ALPHA     1.0
#define VL_MAT       50
#define VL_CL        50
#define VL_RES       1.0e-6
#define VL_MAXIT     1000
#define VL_TGV       1.e+30
#define VL_EPSD      1.e-30
#define VL_EPSA      1.e-200
#define VL_LONMAX    1024

#define VL_MAX(a,b) (((a) < (b)) ? (b) : (a))
#define VL_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define VL_MIN3(a,b,c) ( (a) < (b) ? ((a)<(c) ? (a) : (c)) : ((b)<(c) ? (b) : (c)) )
#define VL_MAX3(a,b,c) ( (a) > (b) ? ((a)>(c) ? (a) : (c)) : ((b)>(c) ? (b) : (c)) )


/* data structures */
typedef struct {
  double    c[3];
  int       s,ref,new,mark;
} Point;
typedef Point * pPoint;

typedef struct {
	int        v[2],ref;
} Edge;
typedef Edge * pEdge;

typedef struct {
  int       v[3],adj[3],ref,mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int       v[4],adj[4],ref,mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  int       dim,ver;
  int       np,npi,na,nai,nt,nti,ne,nei,mark;
  char      verb,zip,ls,mfree;
  mytime    ctim[TIMEMAX];
} Info;

typedef struct {
  int      mark;
  char    *name;
  Point   *point;
	Edge    *edge;
  Tria    *tria;
  Tetra   *tetra;
} Mesh;
typedef Mesh * pMesh;

/* boundary conds. */
typedef struct {
  double   u[3];
  int      ref;
  char     typ,elt,att;
} Cl;
typedef Cl * pCl;

typedef struct {
  double  alpha;
  int     ref;
} Mat;
typedef Mat * pMat;

typedef struct {
  double   *u,*chi,gr[3],alpha,res,tim;
  int       nbcl,nmat,nit,nt;
  char     *namein,*nameout,*namepar,*namechi,cltyp,clelt;
  Cl       *cl;
  Mat      *mat;
} Sol;
typedef Sol * pSol;

struct _VLst {
  Mesh    mesh;
  Sol     sol;
  Info    info;
};


/* prototypes */
int  loadMesh(VLst *vlst);
int  loadSol(VLst *vlst);
int  loadChi(VLst *vlst);
int  saveSol(VLst *vlst);
int  pack_2d(VLst *vlst);
int  pack_3d(VLst *vlst);
int  unpack(VLst *vlst);
int  hashel_2d(VLst *vlst);
int  hashel_3d(VLst *vlst);
int  boulep_2d(pMesh mesh,int start,int ip,int *list);
int  boulep_3d(pMesh mesh,int start,int ip,int *list);
pCl  getCl(pSol sol,int ref,int elt);
int  getMat(pSol sol,int ref,double *alpha);
int  velex1_2d(VLst *vlst);
int  velex1_3d(VLst *vlst);
int  invmat_2d(double m[4],double mi[4]);
int  invmat_3d(double m[9],double mi[9]);

/*
pCl   getCl(pSol sol,int ref,int elt,char typ);

int   velex1_3d(pMesh ,pSol );
int   gradLS_2d(pMesh ,pSol ,int ,char ,double *);
int   gradLS_3d(pMesh ,pSol ,int ,char ,double *);
int   sctove_2d(pMesh ,pSol );
int   sctove_3d(pMesh ,pSol );
*/

#endif
