#include "velext.h"


/* compactify mesh structure */
int pack_3d(VLst *vlst) {
  pTetra    pei,pef;
  pTria     pti,ptf;
  double    alpha;
  int      *perm,i,k,nf;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=vlst->info.ne; k++) {
    pei = &vlst->mesh.tetra[k];
    if ( getMat(&vlst->sol,pei->ref,&alpha) ) {
      nf++;
      for (i=0; i<4; i++)  vlst->mesh.point[pei->v[i]].old = 1;
    }
  }
  if ( nf == vlst->info.ne )  return(-1);

  /* store permutations */
  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  vlst->info.zip = 1;
  perm = (int*)calloc(vlst->info.np+1,sizeof(int));
  assert(perm);

  /* compress and realloc vertices */
  vlst->info.npi = vlst->info.np;
  nf = 0;
  for (k=1; k<=vlst->info.np; k++) {
    if ( vlst->mesh.point[k].old > 0 ) {
      nf++;
      if ( nf < k )
        memcpy(&vlst->mesh.point[nf],&vlst->mesh.point[k],sizeof(Point));
      vlst->mesh.point[nf].old = k;
      perm[k] = nf;
    }
  }
  for (k=nf+1; k<=vlst->info.np; k++)  vlst->mesh.point[k].old = 0;
  vlst->info.np = nf;

  /* compress and renum tetrahedra */
  vlst->info.nei = vlst->info.ne;
  nf = 0;
  for (k=1; k<=vlst->info.ne; k++) {
    pei = &vlst->mesh.tetra[k];
    if ( getMat(&vlst->sol,pei->ref,&alpha) ) {
      nf++;
      if ( nf < k )
        memcpy(&vlst->mesh.tetra[nf],&vlst->mesh.tetra[k],sizeof(Tetra));
      pef = &vlst->mesh.tetra[nf];
      for (i=0; i<4; i++)  pef->v[i] = perm[pef->v[i]];
    }
  }
  vlst->info.ne = nf;

  /* renum triangles */
  vlst->info.nti = vlst->info.nt;
  nf = 0;
  for (k=1; k<=vlst->info.nt; k++) {
    pti = &vlst->mesh.tria[k];
    if ( perm[pti->v[0]] && perm[pti->v[1]] && perm[pti->v[2]] ) {
      nf++;
      if ( nf < k )
        memcpy(&vlst->mesh.tria[nf],&vlst->mesh.tria[k],sizeof(Tria));
      ptf = &vlst->mesh.tria[nf];
      ptf->v[0] = perm[ptf->v[0]];
      ptf->v[1] = perm[ptf->v[1]];
      ptf->v[2] = perm[ptf->v[2]];
    }
  }
  vlst->info.nt = nf;

  /* compress solution (data) */
  if ( vlst->sol.u ) {
		for (k=1; k<=vlst->info.np; k++) {
			vlst->sol.u[3*(k-1)+0] = vlst->sol.u[3*(perm[k]-1)+0];
			vlst->sol.u[3*(k-1)+1] = vlst->sol.u[3*(perm[k]-1)+1];
			vlst->sol.u[3*(k-1)+2] = vlst->sol.u[3*(perm[k]-1)+2];
		}
  }
  free(perm);

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",vlst->info.np);
    if ( vlst->info.na )  fprintf(stdout,", %d edges",vlst->info.na);
    if ( vlst->info.nt )  fprintf(stdout,", %d triangles",vlst->info.nt);
    if ( vlst->info.ne )  fprintf(stdout,", %d tetrahedra",vlst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


int pack_2d(VLst *vlst) {
  pTria     pt;
  pEdge     pa;
  double    alpha;
  int      *perm,i,k,nf,id;

  nf = 0;
  for (k=1; k<=vlst->info.nt; k++) {
    pt = &vlst->mesh.tria[k];
    if ( getMat(&vlst->sol,pt->ref,&alpha) ) {
      nf++;
      for (i=0; i<3; i++)  vlst->mesh.point[pt->v[i]].old = pt->v[i];
    }
  }
  if ( nf == vlst->info.nt )  return(-1);

  /* store permutations */
  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  vlst->info.zip = 1;
  perm = (int*)calloc(vlst->info.np+1,sizeof(int));
  assert(perm);

  /* compress and realloc vertices */
  nf = 0;
  for (k=1; k<=vlst->info.np; k++) {
    perm[k] = k;
    id = vlst->mesh.point[k].old;
    vlst->mesh.point[k].old = k;
    if ( id > 0 ) {
      nf++;
      /* swap k and nf */
      if ( nf < k ) {
        memcpy(&vlst->mesh.point[0],&vlst->mesh.point[nf],sizeof(Point));
        memcpy(&vlst->mesh.point[nf],&vlst->mesh.point[k],sizeof(Point));
        memcpy(&vlst->mesh.point[k],&vlst->mesh.point[0],sizeof(Point));
        perm[k] = nf;
      }
    }
  }
  vlst->info.np = nf;

  /* compress and renum triangles */
  nf = 0;
  for (k=1; k<=vlst->info.nt; k++) {
    pt = &vlst->mesh.tria[k];
    if ( getMat(&vlst->sol,pt->ref,&alpha) ) {
      nf++;
      if ( nf < k ) {
        memcpy(&vlst->mesh.tria[0],&vlst->mesh.tria[nf],sizeof(Tria));
        memcpy(&vlst->mesh.tria[nf],&vlst->mesh.tria[k],sizeof(Tria));
        memcpy(&vlst->mesh.tria[k],&vlst->mesh.tria[0],sizeof(Tria));
      }
    }
  }
  for (k=1; k<=vlst->info.nt; k++) {
    pt = &vlst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = perm[pt->v[i]];
  }
  vlst->info.nt = nf;

  /* renum edges */
  nf = 0;
  for (k=1; k<=vlst->info.na; k++) {
    pa = &vlst->mesh.edge[k];
    if ( perm[pa->v[0]] <= vlst->info.np && perm[pa->v[1]] <= vlst->info.np ) {
      nf++;
      if ( nf < k )
        memcpy(&vlst->mesh.edge[nf],&vlst->mesh.edge[k],sizeof(Edge));
      pa = &vlst->mesh.edge[nf];
      pa->v[0] = perm[pa->v[0]];
      pa->v[1] = perm[pa->v[1]];
    }
  }
  vlst->info.na = nf;

  /* compress solution (data) */
  if ( vlst->sol.u ) {
		for (k=1; k<=vlst->info.npi; k++) {
			vlst->sol.u[2*(perm[k]-1)+0] = vlst->sol.u[2*(k-1)+0];
			vlst->sol.u[2*(perm[k]-1)+1] = vlst->sol.u[2*(k-1)+1];
		}
  }
  free(perm);

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",vlst->info.np);
    if ( vlst->info.na )  fprintf(stdout,", %d edges",vlst->info.na);
    if ( vlst->info.nt )  fprintf(stdout,", %d triangles",vlst->info.nt);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* restore solution at initial vertices */
int unpack(VLst *vlst) {
  pPoint  ppt;
	double *unew;
  int     k;
	char    i;

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Uncompressing data: ");
    fflush(stdout);
  }
	unew = (double*)calloc(vlst->info.dim*vlst->info.npi,sizeof(double));
	assert(unew);

  for (k=1; k<=vlst->info.npi; k++) {
    ppt = &vlst->mesh.point[k];
	  for (i=0; i<vlst->info.dim; i++)
      unew[vlst->info.dim*(ppt->old-1)+i] = vlst->sol.u[vlst->info.dim*(k-1)+i];
  }
	free(vlst->sol.u);
	vlst->sol.u   = unew;
  vlst->info.np = vlst->info.npi;
  vlst->info.na = vlst->info.nai;
  vlst->info.nt = vlst->info.nti;
  vlst->info.ne = vlst->info.nei;

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"%d data vectors\n",vlst->info.np);
  }
  return(1);
}

