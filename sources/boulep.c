#include "velext.h"


/* find all vertices connected to P
   in:  start : tetrahedron containing p 
        ip    : index of p in start
        list  : dynamic list structure (allocated outside)
   out: list  : list of vertices */
int boulep_3d(pMesh mesh,int start,int ip,int *list) {
  pTetra  pt,pt1;
  pPoint  ppt;
  int     adj,i,j,kk,indp,iel,iadr,ilist,iilist,nump;
  int     llist[VL_LONMAX+2];
  static int base = 0;

  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  ppt  = &mesh->point[nump];

  /* store initial tetra */
  base++;
  pt->mark = base;
  ilist    = 0;
  iilist   = 1;
  llist[1] = 4*start + ip;

  /* store first 3 vertices */
	list[0] = nump;
  for (j=0; j<4; j++) {
    if ( j != ip ) {
      ilist++;
      list[ilist] = pt->v[j];
      mesh->point[pt->v[j]].mark = base;
    }
  }

  /* store 3 neighbors sharing P */
  iadr = (start-1)*4 + 1;
  for (i=0; i<4; i++) {
    if ( i == ip )  continue;
    adj = pt->adj[i] / 4;
    if ( adj ) {
      pt1 = &mesh->tetra[adj];
      if ( pt1->mark != base ) {
        /* not stored yet */
        pt1->mark = base;
        for (j=0; j<4; j++)
          if ( pt1->v[j] == nump )  break;
        iilist++;
        llist[iilist] = 4*adj + j;
      }
    }
  }
  if ( iilist < 2 )  return(ilist);

  /* explore list of neighbors */
  indp = 2;
  do {
    iel  = llist[indp] / 4;
    pt   = &mesh->tetra[iel];
    iadr = 4*(iel-1) + 1;
    for (i=0; i<4; i++) {
      if ( pt->v[i] == nump )  continue;
      adj = pt->adj[i] / 4;
      if ( adj ) {
        pt1 = &mesh->tetra[adj];
        if ( pt1->mark != base ) {
          pt1->mark = base;
          for (j=0; j<4; j++)
            if ( pt1->v[j] == nump )  break;
          iilist++;
          llist[iilist] = 4*adj + j;
        }
      }
    }
    /* overflow */
    if ( iilist > VL_LONMAX-3 )  return(-ilist);
  }
  while ( ++indp <= iilist );

  /* store vertices from tetra list */
  for (i=2; i<=iilist; i++) {
    kk = llist[i] / 4;
    pt = &mesh->tetra[kk];
    for (j=0; j<4; j++) {
      if ( pt->v[j] != nump ) {
        ppt = &mesh->point[ pt->v[j] ];
        if ( ppt->mark < base ) {
          ilist++;
          list[ilist] = pt->v[j];
          ppt->mark = base;
        }
      }
    }
  }

  return(ilist);
}


/* return all vertices connected to ip, list[0] = ip */
int boulep_2d(pMesh mesh,int start,int ip,int *list) {
  pTria    pt;
  int      k,ilist;
  char     i,i1,i2;

  pt = &mesh->tria[start];
  list[0] = pt->v[ip];
  ilist   = 0;

  /* store neighbors */
  k  = start;
  i  = ip;
  i1 = (i+1) % 3;
  i2 = (i1+1) % 3;
  do {
    if ( ilist > VL_LONMAX-2 )  return(-ilist);
    ilist++;
    list[ilist] = pt->v[i2];

    k  = pt->adj[i1] / 3;
    i2 = pt->adj[i1] % 3;
    i1 = (i2+2) % 3;
    pt = &mesh->tria[k];
  }
  while ( k && k != start );
  if ( k > 0 )  return(ilist);

  /* reverse loop */
  k  = start;
  i  = ip;
  pt = &mesh->tria[k];
  i1 = (i+1) % 3;
  i2 = (i1+1) % 3;
  do {
    if ( ilist > VL_LONMAX-2 )  return(-ilist);
    ilist++;
    list[ilist] = pt->v[i1];

    k  = pt->adj[i2] / 3;
    i1 = pt->adj[i2] % 3;
    i2 = (i1+2) % 3;
    pt = &mesh->tria[k];
  }
  while ( k > 0 );

  return(ilist);
}




