#include "velext.h"
#include "sparse.h"


/* find boundary conds in list */
pCl getCl(pSol sol,int ref,int elt) {
  pCl     pcl;
  int     i;

  for (i=0; i<sol->nbcl; i++) {
    pcl = &sol->cl[i];
    if ( (pcl->ref == ref) && (pcl->elt == elt) )  return(pcl);
  }
  return(0);
}


/* retrieve physical properties in list */
int getMat(pSol sol,int ref,double *alpha) {
  pMat   pm;
  int    i;

  *alpha = sol->alpha;
  if ( sol->nmat == 0 )  return(1);
  for (i=0; i<sol->nmat; i++) {
    pm = &sol->mat[i];
    if ( pm->ref == ref ) {
      *alpha = pm->alpha;
      return(1);
    }
  }
  *alpha = VL_ALPHA;

  return(0);
}


/* triangle area */
static inline double area_2d(double *a,double *b,double *c) {
  return(0.5 * ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])));
}


int invmat_2d(double m[4],double mi[4]) {
  double   det;

  det = m[0]*m[3] - m[2]*m[1];
  if ( fabs(det) < VL_EPSA )  return(0);
  det = 1.0 / det;
	mi[0] =  m[3]*det;
	mi[1] = -m[1]*det;
	mi[2] = -m[2]*det;
	mi[3] =  m[0]*det;

	return(1);
}


/* set TGV to diagonal coefficient when Dirichlet */
static int setTGV_2d(VLst *vlst,pCsr A) {
  pCl      pcl;
	pEdge    pa;
  pPoint   ppt;
  int      k;

  /* Set Dirichlet's boundary for the state system */
  if ( vlst->sol.clelt & VL_ver ) {
    for (k=1; k<=vlst->info.np; k++) {
      ppt = &vlst->mesh.point[k];
      pcl = getCl(&vlst->sol,ppt->ref,VL_ver);
      if ( pcl && pcl->typ == Dirichlet ) {
        /* set value for scalar or vector field */
        if ( vlst->info.ls )
          csrSet(A,k-1,k-1,VL_TGV);
        else {
          csrSet(A,2*(k-1)+0,2*(k-1)+0,VL_TGV);
          csrSet(A,2*(k-1)+1,2*(k-1)+1,VL_TGV);
        }
			}
    }
  }
	if ( vlst->sol.clelt & VL_edg )	{
    for (k=1; k<=vlst->info.na; k++) {
      pa = &vlst->mesh.edge[k];
			pcl = getCl(&vlst->sol,pa->ref,VL_edg);
      if ( pcl && pcl->typ == Dirichlet ) {
        /* set value for scalar or vector field */
        if ( vlst->info.ls ) {
          csrSet(A,pa->v[0]-1,pa->v[0]-1,VL_TGV);
          csrSet(A,pa->v[1]-1,pa->v[1]-1,VL_TGV);
        }
        else {
          csrSet(A,2*(pa->v[0]-1)+0,2*(pa->v[0]-1)+0,VL_TGV);
          csrSet(A,2*(pa->v[0]-1)+1,2*(pa->v[0]-1)+1,VL_TGV);
          csrSet(A,2*(pa->v[1]-1)+0,2*(pa->v[1]-1)+0,VL_TGV);
          csrSet(A,2*(pa->v[1]-1)+1,2*(pa->v[1]-1)+1,VL_TGV);
        }
      }
		}
	}

	return(1);
}


static pCsr matA1_2d(VLst *vlst) {
  pCsr     A;
  pTria    pt;
  pPoint   p0,p1,p2;
  double   Dp[3][2],m[4],im[4],Gr[3][2],alpha,vol,kij,term0,termG;
  int      nr,nc,nbe,k,ni,nj,il,ic;
  char     i,j;

	/* memory allocation (rough estimate) */
	nr  = nc = vlst->info.np;
  nbe = 10*vlst->info.np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  /* Dp */
  Dp[0][0] = -1.0; Dp[1][0] = 1.0; Dp[2][0] = 0.0;
  Dp[0][1] = -1.0; Dp[1][1] = 0.0; Dp[2][1] = 1.0; 

  /* Fill stiffness matrix of Laplace problem */
  for (k=1; k<=vlst->info.nt; k++) {
    pt = &vlst->mesh.tria[k];
    if ( !getMat(&vlst->sol,pt->ref,&alpha) )  continue;

    p0 = &vlst->mesh.point[pt->v[0]];
    p1 = &vlst->mesh.point[pt->v[1]];
    p2 = &vlst->mesh.point[pt->v[2]];

    m[0] = p1->c[0]-p0->c[0];  m[1] = p1->c[1]-p0->c[1];
    m[2] = p2->c[0]-p0->c[0];  m[3] = p2->c[1]-p0->c[1];
    if ( !invmat_2d(m,im) )  return(0);

    /* volume of element k */
    vol = area_2d(p0->c,p1->c,p2->c);

    /* Gradients of shape functions : Gr[i] = im*Dp[i] */
    for (i=0; i<3; i++) {
      Gr[i][0] = im[0]*Dp[i][0] + im[1]*Dp[i][1];
      Gr[i][1] = im[2]*Dp[i][0] + im[3]*Dp[i][1];
    }

    /* Flow local stiffness matrix into global one */
    for (i=0; i<3; i++) {
      for (j=i; j<3; j++) {
        if ( j==i )
          term0 = vol / 6.0;
        else
          term0 = vol / 12.0;

        termG = vol * (Gr[i][0]*Gr[j][0] + Gr[i][1]*Gr[j][1]);   
        kij = term0 + alpha * termG;
        ni  = pt->v[i]; 
        nj  = pt->v[j];
        if ( ni < nj ) {
          il = ni-1;
          ic = nj-1;
        }
        else {
          il = nj-1; 
          ic = ni-1;
        }
        csrPut(A,il,ic,kij);
      }
    }
  }
  setTGV_2d(vlst,A);
  csrPack(A);

  if ( vlst->info.verb == '+' )
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nr,nc,100.0*A->nbe/nr/nc);

  return(A);
}


/* build stiffness matrix (vector case) */
static pCsr matA2_2d(VLst *vlst) {
  pCsr     A;
  pTria    pt;
  pPoint   p0,p1,p2;
  double   m[4],im[4],Gr[3][3],alpha,cof,vol,kij,term0,termG;
  int      nr,nc,nbe,k,il,ic;
  char     i,j;

  /* memory allocation (rough estimate) */
  nr  = nc = 2*vlst->info.np;
  nbe = 10*vlst->info.np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  /* Fill stiffness matrix of Laplace problem */
  for (k=1; k<=vlst->info.nt; k++) {
    pt = &vlst->mesh.tria[k];
    if ( !getMat(&vlst->sol,pt->ref,&alpha) )  continue;

    p0 = &vlst->mesh.point[pt->v[0]];
    p1 = &vlst->mesh.point[pt->v[1]];
    p2 = &vlst->mesh.point[pt->v[2]];
    
    m[0] = p1->c[0]-p0->c[0];  m[1] = p1->c[1]-p0->c[1];
    m[2] = p2->c[0]-p0->c[0];  m[3] = p2->c[1]-p0->c[1];

    /* volume of element k */
    vol = area_2d(p0->c,p1->c,p2->c);
    cof = 1.0 / (4.0*vol);

    Gr[0][0] = cof*(m[2]-m[0])*(m[2]-m[0])+ cof*(m[3]-m[1])*(m[3]-m[1]);
    Gr[1][1] = cof*(m[2]*m[2]+m[3]*m[3]);
    Gr[2][2] = cof*(m[0]*m[0]+m[1]*m[1]);
    Gr[0][1] = -cof*(m[2]-m[0])*(m[2])-cof*(m[3]-m[1])*(m[3]);
    Gr[0][2] = cof*(m[2]-m[0])*(m[0])+cof*(m[3]-m[1])*(m[1]);
    Gr[1][2] = -cof*(m[2])*(m[0])-cof*(m[3])*(m[1]);

    /* Flow local stiffness matrix into global one */
    for (i=0; i<3; i++) {
      for (j=i; j<3; j++) {
        if ( j==i )
          term0 = vol / 6.0;
        else
          term0 = vol / 12.0;
       
        /* termG = vol * (Gr[i][0]*Gr[j][0] + Gr[i][1]*Gr[j][1]); */
        termG = Gr[i][j];
        kij = term0 + alpha * termG;
        il  = 2*(pt->v[i]-1);
        ic  = 2*(pt->v[j]-1);
        if ( i == j ) {
          csrPut(A,il+0,ic+0,kij);
          csrPut(A,il+1,ic+1,kij);
        }
        else {
          csrPut(A,il+0,ic+0,kij);
          csrPut(A,ic+0,il+0,kij);
          csrPut(A,il+1,ic+1,kij);
          csrPut(A,ic+1,il+1,kij);
        }
      }
    }
  }
  setTGV_2d(vlst,A);
  csrPack(A);

  if ( vlst->info.verb == '+' )
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nr,nc,100.0*A->nbe/nr/nc);

  return(A);
}


/* build right hand side vector and set boundary conds. */
static double *rhsF_2d(VLst *vlst) {
  pPoint   ppt;
  pEdge    pa;
  pTria    pt;
  pCl      pcl;
  double  *F,*vp,*a,*b,*c,area;
  int      i,k,ig,nc1,nc2,nc3;

  if ( vlst->info.verb == '+' )  fprintf(stdout,"     boundary conditions: ");
  if ( vlst->info.ls )
    F = (double*)calloc(vlst->info.np,sizeof(double));
  else
    F = (double*)calloc(2*vlst->info.np,sizeof(double));
  assert(F);

  nc1 = nc2 = nc3 = 0;

  /* gravity as external force */
  if ( vlst->sol.cltyp & Gravity ) {
    for (k=1; k<=vlst->info.nt; k++) {
      pt = &vlst->mesh.tria[k];

      /* measure of K */
      a = &vlst->mesh.point[pt->v[0]].c[0]; 
      b = &vlst->mesh.point[pt->v[1]].c[0]; 
      c = &vlst->mesh.point[pt->v[2]].c[0]; 
      area = area_2d(a,b,c) / 3.0;
      if ( vlst->info.ls ) {
        for (i=0; i<3; i++)
          F[pt->v[i]-1] += area * vlst->sol.gr[0];
      }
      else {
        for (i=0; i<3; i++) {
          F[2*(pt->v[i]-1)+0] += area * vlst->sol.gr[0];
          F[2*(pt->v[i]-1)+1] += area * vlst->sol.gr[1];
        }
      }
      nc1++;
    }
  }

  /* nodal boundary conditions */
  if ( vlst->sol.clelt & VL_ver ) {
	  for (k=1; k<=vlst->info.np; k++) {
	    ppt = &vlst->mesh.point[k];
			pcl = getCl(&vlst->sol,ppt->ref,VL_ver);
	    if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        if ( vlst->info.ls ) {
          vp = pcl->att == 'f' ? &vlst->sol.u[k-1] : &pcl->u[0];
          F[k-1] = VL_TGV * vp[0];
        }
        else {
          vp = pcl->att == 'f' ? &vlst->sol.u[2*(k-1)] : &pcl->u[0];
          F[2*(k-1)+0] = VL_TGV * vp[0];
          F[2*(k-1)+1] = VL_TGV * vp[1];
	      }
      }
      nc2++;
		}
	}

  /* external load along boundary edges */
  if ( vlst->sol.clelt & VL_edg ) {
    for (k=1; k<=vlst->info.na; k++) {
      pa  = &vlst->mesh.edge[k];
      pcl = getCl(&vlst->sol,pa->ref,VL_edg);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        if (vlst->info.ls ) {
          vp = pcl->att == 'f' ? &vlst->sol.u[k-1] : &pcl->u[0];
          F[pa->v[0]-1] = VL_TGV * vp[0];
          F[pa->v[1]-1] = VL_TGV * vp[0];
        }
        else {
          vp = pcl->att == 'f' ? &vlst->sol.u[2*(k-1)] : &pcl->u[0];
          F[2*(pa->v[0]-1)+0] = VL_TGV * vp[0];
          F[2*(pa->v[0]-1)+1] = VL_TGV * vp[1];
          F[2*(pa->v[1]-1)+0] = VL_TGV * vp[0];
          F[2*(pa->v[1]-1)+1] = VL_TGV * vp[1];
        }
        nc3++;
      }
    }
  }
  if ( vlst->info.verb == '+' ) {
    if ( nc1 > 0 )  fprintf(stdout," %d gravity",nc1);
    if ( nc2 > 0 )  fprintf(stdout," %d nodal",nc2);
    if ( nc3 > 0 )  fprintf(stdout," %d load",nc3);
    fprintf(stdout,"\n");
  }

	return(F);
}


/* solve Helmholz */
int velex1_2d(VLst *vlst) {
  pCsr     A;
  double  *F;
  int      ier;
  char     stim[32];

  /* -- Part I: matrix assembly */
  if ( vlst->info.verb != '0' )  fprintf(stdout,"    Matrix and right-hand side assembly\n");

  /* allocating memory (for dylib) */
  if ( !vlst->sol.u ) {
    vlst->sol.u  = (double*)calloc(vlst->info.dim*vlst->info.np,sizeof(double));
    assert(vlst->sol.u);
  }

  /* build matrix */
  A = vlst->info.ls ? matA1_2d(vlst) : matA2_2d(vlst);
  F = rhsF_2d(vlst);

  /* free mesh structure + boundary conditions */
  if ( vlst->info.mfree ) {
		free(vlst->mesh.tria);
    if ( !vlst->info.zip )  free(vlst->mesh.point);
	}

  /* -- Part II: Laplace problem solver */
  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Solving linear system:");  fflush(stdout);
    ier = csrPrecondGrad(A,vlst->sol.u,F,&vlst->sol.res,&vlst->sol.nit,1);
    if ( ier <= 0 )
      fprintf(stdout,"\n # convergence problem: %d (%f,%d)\n",ier,vlst->sol.res,vlst->sol.nit);
    else
      fprintf(stdout," %E in %d iterations\n",vlst->sol.res,vlst->sol.nit);
	}
  else {
    ier = csrPrecondGrad(A,vlst->sol.u,F,&vlst->sol.res,&vlst->sol.nit,1);
  }

  /* free memory */
  csrFree(A);
	free(F);

  return(ier > 0);  
}

