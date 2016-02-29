#include "velext.h"


/* compactify mesh structure */
int pack_3d(VLst *vlst) {
  pTetra    pe;
  pTria     pt;
  pEdge     pa;
  double    alpha,w[3];
  int       i,k,nf,id;

  /* check if compression needed */
  nf  = 0;
  for (k=1; k<=vlst->info.nei; k++) {
    pe = &vlst->mesh.tetra[k];
    if ( getMat(&vlst->sol,pe->ref,&alpha) ) {
      nf++;
      for (i=0; i<4; i++)  vlst->mesh.point[pe->v[i]].new = 1;
    }
  }
  if ( nf == vlst->info.nei )  return(-1);

  /* store permutations */
  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  vlst->info.zip = 1;

  /* compress and renum vertices */
  nf = vlst->info.npi;
  k  = 1;
  while ( k < nf ) {
    if ( vlst->mesh.point[k].new == 0 ) {
      while ( (vlst->mesh.point[nf].new == 0) && (k < nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&vlst->mesh.point[0],&vlst->mesh.point[nf],sizeof(Point));
        memcpy(&vlst->mesh.point[nf],&vlst->mesh.point[k],sizeof(Point));
        memcpy(&vlst->mesh.point[k],&vlst->mesh.point[0],sizeof(Point));
        /* swap solution too */
        if ( vlst->sol.u ) {
          memcpy(&w,&vlst->sol.u[3*(nf-1)],3*sizeof(double));
          memcpy(&vlst->sol.u[3*(nf-1)],&vlst->sol.u[3*(k-1)],3*sizeof(double));
          memcpy(&vlst->sol.u[3*(k-1)],&w,3*sizeof(double));
        }
        vlst->mesh.point[k].new  = nf;
        vlst->mesh.point[nf].new = k;
      }
      nf--;
    }
    k++;
  }
  vlst->info.np = nf;

  /* compress and renum tetrahedra */
  for (k=1; k<=vlst->info.nei; k++) {
    pe = &vlst->mesh.tetra[k];
    for (i=0; i<4; i++)  pe->v[i] = vlst->mesh.point[pe->v[i]].new;
  }
  nf = vlst->info.nei;
  k  = 1;
  while ( k < nf ) {
    pe = &vlst->mesh.tetra[k];
    for (i=0; i<4; i++)  
      if ( pe->v[i] > vlst->info.np || pe->v[i] == 0 )  break;
    if ( i < 4 ) {
      do {
        pe = &vlst->mesh.tetra[nf];
        for (i=0; i<4; i++)
          if ( pe->v[i] > vlst->info.np || pe->v[i] == 0 )  break;
        if ( i == 4 )  break;
        nf --;
      }
      while ( k < nf );
      /* put nf into k */
      memcpy(&vlst->mesh.tetra[k],&vlst->mesh.tetra[nf],sizeof(Tetra));
      nf--;
    }
    k++;
  }
  vlst->info.ne = nf;

  /* renum triangles */
  for (k=1; k<=vlst->info.nti; k++) {
    pt = &vlst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = vlst->mesh.point[pt->v[i]].new;
  }
  nf = vlst->info.nti;
  k  = 1;
  while ( k < nf ) {
    pt = &vlst->mesh.tria[k];
    for (i=0; i<3; i++)
      if ( pt->v[i] > vlst->info.np || pt->v[i] == 0 )  break;
    if ( i < 3 ) {
      do {
        pt = &vlst->mesh.tria[nf];
        for (i=0; i<3; i++)
          if ( pt->v[i] > vlst->info.np || pt->v[i] == 0 )  break;
        if ( i == 3 )  break;
        nf--;
      }
      while ( k < nf );
      /* put nf in k */
      memcpy(&vlst->mesh.tria[k],&vlst->mesh.tria[nf],sizeof(Tria));
      nf--;
    }
    k++;
  }
  vlst->info.nt = nf;

  /* renum edges */
  for (k=1; k<=vlst->info.na; k++) {
    pa = &vlst->mesh.edge[k];
    for (i=0; i<2; i++)  pa->v[i] = vlst->mesh.point[pa->v[i]].new;
  }
  nf = vlst->info.nai;
  k  = 1;
  while ( k < nf ) {
    pa = &vlst->mesh.edge[k];
    if ( (pa->v[0] > vlst->info.np || pa->v[0] == 0) ||
         (pa->v[1] > vlst->info.np || pa->v[1] == 0) ) {
      do {
        pa = &vlst->mesh.edge[nf];
        if ( (pa->v[0] > vlst->info.np || pa->v[0] == 0 ||
              pa->v[1] > vlst->info.np || pa->v[1] == 0) )  break;
        nf--;
      }
      while ( k < nf );
      /* put nf in k */
      memcpy(&vlst->mesh.edge[k],&vlst->mesh.edge[nf],sizeof(Edge));
      nf--;
    }
    k++;
  }
  vlst->info.na = nf;

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"%d vertices",vlst->info.np);
    if ( vlst->info.na )  fprintf(stdout,", %d edges",vlst->info.na);
    if ( vlst->info.nt )  fprintf(stdout,", %d triangles",vlst->info.nt);
    if ( vlst->info.ne )  fprintf(stdout,", %d tetrahedra",vlst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* mesh renumbering and packing */
int pack_2d(VLst *vlst) {
  pTria     pt;
  pEdge     pa;
  double    alpha,w[2];
  int       i,k,nf,id;

  /* check if compression needed */
  nf = 0;
  for (k=1; k<=vlst->info.nti; k++) {
    pt = &vlst->mesh.tria[k];
    if ( getMat(&vlst->sol,pt->ref,&alpha) ) {
      nf++;
      for (i=0; i<3; i++)  vlst->mesh.point[pt->v[i]].new = 1;
    }
  }
  if ( nf == vlst->info.nti )  return(-1);

  /* store permutations */
  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Compressing mesh: ");
    fflush(stdout);
  }
  vlst->info.zip = 1;

  /* compress and renum vertices */
  nf = vlst->info.npi;
  k  = 1;
  while ( k < nf ) {
    if ( vlst->mesh.point[k].new == 0 ) {
      while ( (vlst->mesh.point[nf].new == 0) && (k < nf) )  nf--;
      if ( k < nf ) {
        /* swap k and nf */
        memcpy(&vlst->mesh.point[0],&vlst->mesh.point[nf],sizeof(Point));
        memcpy(&vlst->mesh.point[nf],&vlst->mesh.point[k],sizeof(Point));
        memcpy(&vlst->mesh.point[k],&vlst->mesh.point[0],sizeof(Point));
        /* swap solution too */
        if ( vlst->sol.u ) {
          memcpy(&w,&vlst->sol.u[2*(nf-1)],2*sizeof(double));
          memcpy(&vlst->sol.u[2*(nf-1)],&vlst->sol.u[2*(k-1)],2*sizeof(double));
          memcpy(&vlst->sol.u[2*(k-1)],&w,2*sizeof(double));
        }
        vlst->mesh.point[k].new  = nf;
        vlst->mesh.point[nf].new = k;
      }
      nf--;
    }
    k++;
  }
  vlst->info.np = nf;

  /* compress and renum triangles */
  for (k=1; k<=vlst->info.nti; k++) {
    pt = &vlst->mesh.tria[k];
    for (i=0; i<3; i++)  pt->v[i] = vlst->mesh.point[pt->v[i]].new;
  }
  nf = vlst->info.nti;
  k  = 1;
  while ( k < nf ) {
    pt = &vlst->mesh.tria[k];
    for (i=0; i<3; i++)
      if ( pt->v[i] > vlst->info.np || pt->v[i] == 0 )  break;
    if ( i < 3 ) {
      do {
        pt = &vlst->mesh.tria[nf];
        for (i=0; i<3; i++)
          if ( pt->v[i] > vlst->info.np || pt->v[i] == 0 )  break;
        if ( i < 3 )  break;
        nf--;
      }
      while ( k < nf );
      /* put nf into k */
      memcpy(&vlst->mesh.tria[k],&vlst->mesh.tria[nf],sizeof(Tria));
      nf--;
    }
    k++;
  }
  vlst->info.nt = nf;

  /* compress and renum edges */
  for (k=1; k<=vlst->info.na; k++) {
    pa = &vlst->mesh.edge[k];
    for (i=0; i<2; i++)  pa->v[i] = vlst->mesh.point[pa->v[i]].new;
  }
  nf = vlst->info.nai;
  k  = 1;
  while ( k < nf ) {
    pa = &vlst->mesh.edge[k];
    if ( (pa->v[0] > vlst->info.np || pa->v[0] == 0) ||
         (pa->v[1] > vlst->info.np || pa->v[1] == 0) ) {
      do {
        pa = &vlst->mesh.edge[nf];
        if ( (pa->v[0] > vlst->info.np || pa->v[0] == 0 ||
              pa->v[1] > vlst->info.np || pa->v[1] == 0) )  break;
        nf--;
      }
      while ( k < nf );
      /* put nf in k */
      memcpy(&vlst->mesh.edge[k],&vlst->mesh.edge[nf],sizeof(Edge));
      nf--;
    }
    k++;
  }
  vlst->info.na = nf;

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
	double  w[3];
  int     k,dim;
	char    i;

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"    Uncompressing data: ");
    fflush(stdout);
  }

  dim = vlst->info.dim;
  for (k=1; k<=vlst->info.np; k++) {
    ppt = &vlst->mesh.point[k];
    if ( ppt->new != k ) {
      memcpy(&w,&vlst->sol.u[dim*(k-1)+0],dim*sizeof(double));
      memcpy(&vlst->sol.u[dim*(k-1)+0],&vlst->sol.u[dim*(ppt->new-1)+0],dim*sizeof(double));
      memcpy(&vlst->sol.u[dim*(ppt->new-1)+0],&w,dim*sizeof(double));
    }
  }
  vlst->info.np = vlst->info.npi;
  vlst->info.na = vlst->info.nai;
  vlst->info.nt = vlst->info.nti;
  vlst->info.ne = vlst->info.nei;

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout,"%d data vectors\n",vlst->info.np);
  }

  return(1);
}

