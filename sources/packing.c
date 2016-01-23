#include "velext.h"


/* compactify mesh structure */
int pack_3d(VLst *vlst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  double    alpha;
  int      *perm,i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=vlst->info.ne; k++) {
    pe = &vlst->mesh.tetra[k];
    if ( getMat(&vlst->sol,pe->ref,&alpha) ) {
      nf++;
      for (i=0; i<4; i++)  vlst->mesh.point[pe->v[i]].old = pe->v[i];
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

  /* compress and renum vertices */
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

  /* compress and renum tetrahedra */
  nf = 0;
  for (k=1; k<=vlst->info.ne; k++) {
    pe = &vlst->mesh.tetra[k];
    if ( getMat(&vlst->sol,pe->ref,&alpha) ) {
      nf++;
      if ( nf < k ) {
        memcpy(&vlst->mesh.tria[0],&vlst->mesh.tria[nf],sizeof(Tetra));
        memcpy(&vlst->mesh.tria[nf],&vlst->mesh.tria[k],sizeof(Tetra));
        memcpy(&vlst->mesh.tria[k],&vlst->mesh.tria[0],sizeof(Tetra));
      }
    }
  }
  for (k=1; k<=vlst->info.ne; k++) {
    pe = &vlst->mesh.tetra[k];
    for (i=0; i<4; i++)  pe->v[i] = perm[pe->v[i]];
  }
  vlst->info.ne = nf;

  /* compress solution (data) */
  if ( vlst->sol.u ) {
		for (k=1; k<=vlst->info.npi; k++) {
			vlst->sol.u[2*(perm[k]-1)+0] = vlst->sol.u[2*(k-1)+0];
			vlst->sol.u[2*(perm[k]-1)+1] = vlst->sol.u[2*(k-1)+1];
			vlst->sol.u[2*(perm[k]-1)+2] = vlst->sol.u[2*(k-1)+2];
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

  /* check if compression needed */
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

  /* compress and renum vertices */
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

