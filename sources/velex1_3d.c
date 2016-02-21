#include "velext.h"
#include "sparse.h"


/* tetrahedron volume */
static double volume(double *a,double *b,double *c,double *d) {
  double    bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,vol;

  bx = b[0] - a[0];  by = b[1] - a[1];  bz = b[2] - a[2];
  cx = c[0] - a[0];  cy = c[1] - a[1];  cz = c[2] - a[2];
  dx = d[0] - a[0];  dy = d[1] - a[1];  dz = d[2] - a[2];

  vx  = cy*dz - cz*dy;
  vy  = cz*dx - cx*dz;
  vz  = cx*dy - cy*dx; 
  vol = fabs(bx*vx + by*vy + bz*vz) / 6.0;

  return(vol);
}


/* compute 3d triangle area */
static double area_3d(double *a,double *b,double *c) {
  double    aa,bb,cc,ux,uy,uz,vx,vy,vz,aire;

  ux = b[0] - a[0];  uy = b[1] - a[1];  uz = b[2] - a[2];
  vx = c[0] - a[0];  vy = c[1] - a[1];  vz = c[2] - a[2];

  aa = uy*vz - uz*vy;
  bb = uz*vx - ux*vz;
  cc = ux*vy - uy*vx;
  aire = sqrt(aa*aa + bb*bb + cc*cc);

  return(aire /2.0);
}


/* invert 3x3 non-symmetric matrix */
int invmat_3d(double m[9],double mi[9]) {
  double  aa,bb,cc,det,vmin,vmax,maxx;
  int     k;

  /* check ill-conditionned matrix */
  vmin = vmax = fabs(m[0]);
  for (k=1; k<9; k++) {
    maxx = fabs(m[k]);
    if ( maxx < vmin )  vmin = maxx;
    else if ( maxx > vmax )  vmax = maxx;
  }
  if ( vmax == 0.0 )  return(0);

  /* compute sub-dets */
  aa = m[4]*m[8] - m[5]*m[7];
  bb = m[5]*m[6] - m[3]*m[8];
  cc = m[3]*m[7] - m[4]*m[6];
  det = m[0]*aa + m[1]*bb + m[2]*cc;
  if ( fabs(det) < VL_EPSD )  return(0);
  det = 1.0 / det;

  mi[0] = aa*det;
  mi[3] = bb*det;
  mi[6] = cc*det;
  mi[1] = (m[2]*m[7] - m[1]*m[8])*det;
  mi[4] = (m[0]*m[8] - m[2]*m[6])*det;
  mi[7] = (m[1]*m[6] - m[0]*m[7])*det;
  mi[2] = (m[1]*m[5] - m[2]*m[4])*det;
  mi[5] = (m[2]*m[3] - m[0]*m[5])*det;
  mi[8] = (m[0]*m[4] - m[1]*m[3])*det;

  return(1);
}


/* set TGV to diagonal coefficient when Dirichlet */
static int setTGV_3d(VLst *vlst,pCsr A) {
  pCl      pcl;
	pTria    pt;
  pPoint   ppt;
  int      k,ig;
	char     i;

  /* Set Dirichlet's boundary for the state system */
  if ( vlst->sol.clelt & VL_ver ) {
    /* at mesh nodes */
    for (k=1; k<=vlst->info.np; k++) {
      ppt = &vlst->mesh.point[k];
      pcl = getCl(&vlst->sol,ppt->ref,VL_ver);
      if ( pcl && pcl->typ == Dirichlet ) {
        /* set value for scalar or vector field */
        if ( vlst->info.ls ) {
          csrSet(A,k-1,k-1,VL_TGV);
        }
        else {
          ig = 3*(k-1);
          csrSet(A,ig+0,ig+0,VL_TGV);
          csrSet(A,ig+1,ig+1,VL_TGV);
          csrSet(A,ig+2,ig+2,VL_TGV);
        }
      }
		}
	}

	return(1);
}


static pCsr matA1_3d(VLst *vlst) {
  pCsr     A;
  pTetra   pt;
  pPoint   p0,p1,p2,p3;
  double   Dp[4][3],m[9],im[9],Gr[4][3],alpha,vol,kij,term0,termG;
  int      nr,nc,nbe,k,ip0,ip1,ip2,ip3,ni,nj,il,ic;
  char     i,j;

	/* memory allocation (rough estimate) */
	nr  = nc = vlst->info.np;
  nbe = 20*vlst->info.np;
  A   = csrNew(nr,nc,nbe,CS_UT+CS_SYM);

  /* Dp */
  Dp[0][0] = -1.0 ; Dp[1][0] = 1.0 ; Dp[2][0] = 0.0 ; Dp[3][0] = 0.0;
  Dp[0][1] = -1.0 ; Dp[1][1] = 0.0 ; Dp[2][1] = 1.0 ; Dp[3][1] = 0.0; 
  Dp[0][2] = -1.0 ; Dp[1][2] = 0.0 ; Dp[2][2] = 0.0 ; Dp[3][2] = 1.0;

  /* Fill stiffness matrix of Laplace problem */
  for (k=1; k<=vlst->info.ne; k++) {
    pt = &vlst->mesh.tetra[k];
    if ( !getMat(&vlst->sol,pt->ref,&alpha) )  continue;

    p0 = &vlst->mesh.point[pt->v[0]];
    p1 = &vlst->mesh.point[pt->v[1]];
    p2 = &vlst->mesh.point[pt->v[2]];
    p3 = &vlst->mesh.point[pt->v[3]];

    m[0] = p1->c[0]-p0->c[0];  m[1] = p1->c[1]-p0->c[1];  m[2] = p1->c[2]-p0->c[2];
    m[3] = p2->c[0]-p0->c[0];  m[4] = p2->c[1]-p0->c[1];  m[5] = p2->c[2]-p0->c[2];
    m[6] = p3->c[0]-p0->c[0];  m[7] = p3->c[1]-p0->c[1];  m[8] = p3->c[2]-p0->c[2];
    if ( !invmat_3d(m,im) )  return(0);

    /* volume of element k */
    vol = volume(p0->c,p1->c,p2->c,p3->c);

    /* Gradients of shape functions : Gr[i] = im*Dp[i] */
    for (i=0; i<4; i++) {
      Gr[i][0] = im[0]*Dp[i][0] + im[1]*Dp[i][1] + im[2]*Dp[i][2];
      Gr[i][1] = im[3]*Dp[i][0] + im[4]*Dp[i][1] + im[5]*Dp[i][2];
      Gr[i][2] = im[6]*Dp[i][0] + im[7]*Dp[i][1] + im[8]*Dp[i][2];
    }

    /* Flow local stiffness matrix into global one */
    for (i=0; i<4; i++) {
      for (j=i; j<4; j++) {
        if ( j==i )
          term0 = vol / 10.0;  
        else
          term0 = vol / 20.0;

        termG = vol * (Gr[i][0]*Gr[j][0] + Gr[i][1]*Gr[j][1] + Gr[i][2]*Gr[j][2]);   
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

  setTGV_3d(vlst,A);
	csrPack(A);

  if ( vlst->info.verb == '+' )
    fprintf(stdout,"     %dx%d matrix, %.2f sparsity\n",nr,nc,100.0*A->nbe/nr/nc);

  return(A);
}


static pCsr matA2_3d(VLst *vlst) {
  pCsr    A;
  return(A);
}


/* build right hand side vector and set boundary conds. */
static double *rhsF_3d(VLst *vlst) {
	pTetra   pt;
  pTria    ptt;
  pPoint   ppt;
  pCl      pcl;
  double  *F,*vp,w[3],*a,*b,*c,*d,vol;
  int      k,il,nc;
  char     i;

  if ( vlst->info.verb == '+' )  fprintf(stdout,"     boundary conditions: ");
  if ( vlst->info.ls )
    F = (double*)calloc(vlst->info.np,sizeof(double));
  else
    F = (double*)calloc(3*vlst->info.np,sizeof(double));
  assert(F);

  /* gravity as external force */
  if ( vlst->sol.cltyp & Gravity ) {
    nc = 0;
    for (k=1; k<=vlst->info.ne; k++) {
      pt = &vlst->mesh.tetra[k];

      /* measure of K */
      a = &vlst->mesh.point[pt->v[0]].c[0]; 
      b = &vlst->mesh.point[pt->v[1]].c[0]; 
      c = &vlst->mesh.point[pt->v[2]].c[0];
      d = &vlst->mesh.point[pt->v[3]].c[0];
      vol = volume(a,b,c,d) / 4.0;
      if ( vlst->info.ls ) {
        for (i=0; i<4; i++)
          F[pt->v[i]-1] += vol * vlst->sol.gr[0];
      }
      else {
        for (i=0; i<4; i++) {
          F[3*(pt->v[i]-1)+0] += vol * vlst->sol.gr[0];
          F[3*(pt->v[i]-1)+1] += vol * vlst->sol.gr[1];
          F[3*(pt->v[i]-1)+2] += vol * vlst->sol.gr[2];
        }
      }
      nc++;
    }
    if ( vlst->info.verb == '+' )  fprintf(stdout,"     %d gravity values assigned\n",nc);
  }

  /* nodal boundary conditions for interface */
  if ( vlst->sol.clelt & VL_ver ) {
	  for (k=1; k<=vlst->info.np; k++) {
	    ppt = &vlst->mesh.point[k];
      if ( !ppt->ref )  continue;
			pcl = getCl(&vlst->sol,ppt->ref,VL_ver);
	    if ( !pcl )  continue;
      if ( vlst->info.ls ) {
        vp = pcl->att == 'f' ? &vlst->sol.u[k-1] : &pcl->u[0];
        F[k-1] = VL_TGV * vp[0];
      }
      else {
        vp = pcl->att == 'f' ? &vlst->sol.u[2*(k-1)] : &pcl->u[0];
        F[3*(k-1)+0] = VL_TGV * vp[0];
        F[3*(k-1)+1] = VL_TGV * vp[1];
        F[3*(k-1)+2] = VL_TGV * vp[2];
	    }
		}
	}

  /* external load along boundary triangles */
  if ( vlst->sol.clelt & VL_tri ) {
    nc = 0;
    for (k=1; k<=vlst->info.nt; k++) {
      ptt = &vlst->mesh.tria[k];
      pcl = getCl(&vlst->sol,ptt->ref,VL_tri);
      if ( !pcl )  continue;
      else if ( pcl->typ == Dirichlet ) {
        if ( vlst->info.ls ) {
          vp = pcl->att == 'f' ? &vlst->sol.u[k-1] : &pcl->u[0];
          w[0] = VL_TGV * vp[0];
          F[ptt->v[0]-1] = w[0];
          F[ptt->v[1]-1] = w[0];
          F[ptt->v[2]-1] = w[0];
        }
        else {
          vp = pcl->att == 'f' ? &vlst->sol.u[3*(k-1)] : &pcl->u[0];
          w[0] = VL_TGV * vp[0];
          w[1] = VL_TGV * vp[1];
          w[2] = VL_TGV * vp[2];
          for (i=0; i<3; i++) {
            F[2*(ptt->v[i]-1)+0] = w[0];
            F[2*(ptt->v[i]-1)+1] = w[1];
            F[2*(ptt->v[i]-1)+2] = w[2];
          }
        }
        nc++;
      }
    }
    if ( vlst->info.verb == '+' && nc > 0 )  fprintf(stdout,"     %d load values\n",nc);
  }

  return(F);
}


/* solve Helmholz */
int velex1_3d(VLst *vlst) {
  pCsr     A;
  double  *F,err;
  int      nit,ier;
  char     stim[32];

  /* -- Part I: matrix assembly */
  if ( vlst->info.verb != '0' )  fprintf(stdout,"    Matrix and right-hand side assembly\n");
	
  /* build matrix */
  A = vlst->info.ls ? matA1_3d(vlst) : matA2_3d(vlst);
  F = rhsF_3d(vlst);

  /* free mesh structure + boundary conditions */
  if ( vlst->info.mfree ) {
		free(vlst->mesh.tria);
    free(vlst->mesh.tetra);
    if ( !vlst->info.zip )  free(vlst->mesh.point);
	}

  /* -- Part II: Laplace problem solver */
  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Solving linear system:");  fflush(stdout);
    ier = csrPrecondGrad(A,vlst->sol.u,F,&vlst->sol.res,&vlst->sol.nit,1);
    if ( ier <= 0 )
      fprintf(stdout,"\n # convergence problem: %d\n",ier);
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


