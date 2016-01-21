#include "velext.h"
#include "vl_calls.h"
#include "libmesh5.h"


int loadMesh(VLst *vlst) {
  pPoint     ppt;
  pTetra     pt;
  pTria      pt1;
  pEdge      pa;
  float      fp1,fp2,fp3;
  int        k,inm;
  char      *ptr,data[128];

  strcpy(data,vlst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if ( !(inm = GmfOpenMesh(data,GmfRead,&vlst->info.ver,&vlst->info.dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&vlst->info.ver,&vlst->info.dim)) ) {
        fprintf(stderr," # %s file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&vlst->info.ver,&vlst->info.dim)) ) {
    fprintf(stderr," # %s file not found.\n",data);
    return(0);
  }

  if ( vlst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  vlst->info.np = GmfStatKwd(inm,GmfVertices);
	vlst->info.na = GmfStatKwd(inm,GmfEdges);
  vlst->info.nt = GmfStatKwd(inm,GmfTriangles);
  vlst->info.ne = GmfStatKwd(inm,GmfTetrahedra);

  if ( !vlst->info.np ) {
    if ( vlst->info.verb != '0' )  fprintf(stdout,"\n # missing data\n");
    return(0);
  }
	vlst->info.npi = vlst->info.np;
  vlst->info.nai = vlst->info.na;
	vlst->info.nti = vlst->info.nt;
	vlst->info.nei = vlst->info.ne;

  /* memory alloc */
  vlst->mesh.point = (Point*)calloc(vlst->info.np+1,sizeof(Point));
  assert(vlst->mesh.point);
  if ( vlst->info.nt > 0 ) {
    vlst->mesh.tria  = (Tria*)calloc(vlst->info.nt+1,sizeof(Tria));
    assert(vlst->mesh.tria);
  }
  if ( vlst->info.ne > 0 ) {
    vlst->mesh.tetra  = (Tetra*)calloc(vlst->info.ne+1,sizeof(Tetra));
    assert(vlst->mesh.tetra);
  }

  /* 2d mesh */
  if ( vlst->info.dim == 2 ) {
    GmfGotoKwd(inm,GmfVertices);
    for (k=1; k<=vlst->info.np; k++) {
      ppt = &vlst->mesh.point[k];
      if ( vlst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->ref);
    }
    /* read mesh edges */
    if ( vlst->info.na > 0 ) {
      vlst->mesh.edge = (Edge*)calloc(vlst->info.na+1,sizeof(Edge));
      assert(vlst->mesh.edge);
    }
    GmfGotoKwd(inm,GmfEdges);
    for (k=1; k<=vlst->info.na; k++) {
      pa = &vlst->mesh.edge[k];
      GmfGetLin(inm,GmfEdges,&pa->v[0],&pa->v[1],&pa->ref);
    }
    /* read triangles */
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=vlst->info.nt; k++) {
      pt1 = &vlst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    }
  }
  /* 3d mesh */
  else {
    GmfGotoKwd(inm,GmfVertices);
    for (k=1; k<=vlst->info.np; k++) {
      ppt = &vlst->mesh.point[k];
      if ( vlst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ppt->ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ppt->ref);
    }
    /* read triangles */
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=vlst->info.nt; k++) {
      pt1 = &vlst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&pt1->ref);
    }
    /* read tetrahedra */
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=vlst->info.ne; k++) {
      pt = &vlst->mesh.tetra[k];
  		pt->mark =0;
      GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&pt->ref);
    }
  }
  GmfCloseMesh(inm);

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout," %d vertices",vlst->info.np);
    if ( vlst->info.na )  fprintf(stdout,", %d edges",vlst->info.na);
    if ( vlst->info.nt )  fprintf(stdout,", %d triangles",vlst->info.nt);
    if ( vlst->info.ne )  fprintf(stdout,", %d tetrahedra",vlst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* load initial solution (vector field) */
int loadSol(VLst *vlst) {
  double       bufd[GmfMaxTyp];
  float        buf[GmfMaxTyp];
  int          i,k,dim,ver,np,inm,type,size,offset,typtab[GmfMaxTyp];
  char        *ptr,data[128];

	if ( !vlst->sol.namein )  return(-1);
  strcpy(data,vlst->sol.namein);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);

  if ( dim != vlst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&offset,&typtab);
  if ( !np || typtab[0] != GmfVec || np != vlst->info.np )  return(-1);

  if ( vlst->info.verb != '0' )  fprintf(stdout,"    %s :",data);

  /* read sol: assume velocity if 1st field */
  GmfGotoKwd(inm,GmfSolAtVertices);
  if ( ver == GmfFloat ) {
	  for (k=0; k<vlst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,&buf);
      for (i=0; i<dim; i++)
        vlst->sol.u[dim*k+i] = buf[i];
    }
  }
  else {
	  for (k=0; k<vlst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,bufd);
      for (i=0; i<dim; i++) {
        vlst->sol.u[dim*k+i] = bufd[i];
			}
    }
  }
  GmfCloseMesh(inm);

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout," %d data vectors\n",vlst->info.np);
  }

  return(1);
}


/* load characteristic function */
int loadChi(VLst *vlst) {
  double       bufd[GmfMaxTyp];
  float        buf[GmfMaxTyp];
  int          k,inm,np,ver,dim,type,size,typtab[GmfMaxTyp];
  char        *ptr,data[128];

	if ( !vlst->sol.namein )  return(-1);
  strcpy(data,vlst->sol.namein);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".chi.sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".chi.solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".chi.solb");
      *ptr = '\0';
      strcat(data,".chi.sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);

  if ( dim != vlst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&size,typtab);
  if ( !np || typtab[0] != GmfSca || np != vlst->info.np )  return(-1);

  if ( vlst->info.verb != '0' )  fprintf(stdout,"    %s :",data);

  /* read characteristic function */
  GmfGotoKwd(inm,GmfSolAtVertices);
  if ( ver == GmfFloat ) {
    for (k=0; k<vlst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,buf);
      vlst->sol.chi[k] = buf[0];
    }
	}
  else {
    for (k=0; k<vlst->info.np; k++) {
      GmfGetLin(inm,GmfSolAtVertices,bufd);
      vlst->sol.chi[k] = bufd[0];
    }
  }
  GmfCloseMesh(inm);

  if ( vlst->info.verb != '0' ) {
    fprintf(stdout," %d data scalar\n",vlst->info.np);
  }

  return(1);
}


int saveSol(VLst *vlst) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          i,k,outm,type,typtab[GmfMaxTyp];
  char        *ptr,data[128];

  strcpy(data,vlst->sol.nameout);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".new.solb");
  }
  else {
    ptr = strstr(data,".sol");
    if ( !ptr )  strcat(data,".new.sol");
  }

  if ( !(outm = GmfOpenMesh(data,GmfWrite,vlst->info.ver,vlst->info.dim)) ) {
    fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( vlst->info.verb != '0' )  fprintf(stdout,"    %s:",data);
  type = 1;
  typtab[0] = GmfVec;

  /* write sol */
  GmfSetKwd(outm,GmfSolAtVertices,vlst->info.np,type,typtab);
  if ( vlst->info.ver == GmfFloat ) {
    for (k=0; k<vlst->info.np; k++) {
      for (i=0; i<vlst->info.dim; i++)
		    fbuf[i] = vlst->sol.u[vlst->info.dim*k+i];
      GmfSetLin(outm,GmfSolAtVertices,fbuf);
    }
	}
  else {
    for (k=0; k<vlst->info.np; k++) {
      for (i=0; i<vlst->info.dim; i++)
			  dbuf[i] = vlst->sol.u[vlst->info.dim*k+i];
      GmfSetLin(outm,GmfSolAtVertices,dbuf);
    }
  }
  GmfCloseMesh(outm);

  if ( vlst->info.verb != '0' )  fprintf(stdout," %d data vectors\n",vlst->info.np);

  return(1);
}
