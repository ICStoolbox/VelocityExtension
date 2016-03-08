/*
 * main program file for velext
 * (C) Copyright 2014- , ICS-SU
 *
 * This file is part of velext.
 *
 * velext is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * velext is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with velext.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "velext.h"
#include "vl_calls.h"


static void excfun(int sigid) {
  fprintf(stdout,"\n Unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout,"  Abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout,"  Code error...\n");  break;
    case SIGFPE:
      fprintf(stdout,"  Floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout,"  Illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout,"  Segmentation fault.\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout,"  Programm killed.\n");  break;
  }
  exit(1);
}


static void usage(char *prog) {
  fprintf(stdout,"usage: %s [+/-v | -h] [-a val] [-n nit] [-r res] data_mesh[.mesh] [-c chi[.sol]] [-p param_file[.velext]] [-s data_vel[.sol]] [-o output[.sol]]\n",prog);
  fprintf(stdout,"\nOptions and flags:\n\
  --help       show the syntax and exit.\n\
  --version    show the version and date of release and exit.\n\n\
  -a val       diffusion coefficient\n\
  -n nit       number of iterations max for convergence\n\
  -r res       value of the residual (Krylov space) for convergence\n\
  -v           suppress any message (for use with function call).\n\
  +v           increase the verbosity level for output.\n\n\
  source.mesh    name of the mesh file\n\
  chi.sol        characteristic function (scalar)\n\
  param.velext   name of file containing elasticity parameters\n\
  data.sol       name of file containing the initial solution or boundary conditions\n\
  output.sol     name of the output file (velocity field)\n");
  exit(1);
}


/* parse command line */
static int parsar(int argc,char *argv[],VLst *vlst) {
  int      i;
  char    *ptr,*data;
  
  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i] == '+') ) {
      switch(argv[i][1]) {
      case '-':
        if ( !strcmp(argv[i],"--help") )
          usage(argv[0]);
        else if ( !strcmp(argv[i],"--version") ) {
          fprintf(stdout,"%s: version: %s release: %s\n",argv[0],VL_VER,VL_REL);
          exit(1);
        }
        break;
      case 'h':
      case '?':
        usage(argv[0]);
        break;
			case 'a':
        if ( ++i < argc && isdigit(argv[i][0]) )
          vlst->sol.alpha = atof(argv[i]);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
				break;
      case 'c':
        if ( ++i < argc ) {
          vlst->sol.namechi = argv[i];
          ptr = strstr(vlst->sol.namechi,".sol");
          if ( !ptr )  strcat(vlst->sol.namechi,".sol");
          vlst->info.ls = 1;
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;  
      case 'i':
        if ( ++i < argc ) {
          vlst->mesh.name = argv[i];
          ptr = strstr(vlst->mesh.name,".mesh");
          if ( !ptr )  strcat(vlst->mesh.name,".mesh");
        }
        else {
          fprintf(stdout,"%s: missing input file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'n':
        if ( ++i < argc && isdigit(argv[i][0]) )
          vlst->sol.nit = atoi(argv[i]);
        else {
          fprintf(stdout,"%s: missing input file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'o':
        if ( ++i < argc ) {
          vlst->sol.nameout = argv[i];
          ptr = strstr(vlst->sol.nameout,".sol");
          if ( !ptr )  strcat(vlst->sol.nameout,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'p':
        if ( ++i < argc ) {
          vlst->sol.namepar = argv[i];
          ptr = strstr(vlst->sol.namepar,".elas");
          if ( !ptr )  strcat(vlst->sol.namepar,".elas");
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'r':
        if ( ++i && isdigit(argv[i][0]) )
          vlst->sol.res = strtod(argv[i],NULL);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
        break;
      case 's':
        if ( ++i < argc ) {
          vlst->sol.namein = argv[i];
          ptr = strstr(vlst->sol.namein,".sol");
          if ( !ptr )  strcat(vlst->sol.namein,".sol");
        }
        else {
          fprintf(stdout,"%s: missing data file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'v':
        if ( !strcmp(argv[i],"-v") )
          vlst->info.verb = '0';
        else if ( !strcmp(argv[i],"+v") )
          vlst->info.verb = '+';
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stdout,"%s: illegal option -- %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( vlst->mesh.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        vlst->mesh.name = data;
      }
      else {
        fprintf(stdout,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( vlst->mesh.name == NULL ) {
    fprintf(stderr,"%s: missing argument\n",argv[0]);
    usage(argv[0]);
  }

  return(1);
}


static int parsop(VLst *vlst) {
  Cl         *pcl;
  Mat        *pm;
  int         i,j,npar,ncld,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  if ( !vlst->sol.namepar ) {
    strcpy(data,vlst->mesh.name);
    ptr = strstr(data,".mesh");
    if ( ptr )  *ptr = '\0';
    strcat(data,".velext");
    in = fopen(data,"r");
    if ( !in ) {
      sprintf(data,"%s","DEFAULT.velext");
      in = fopen(data,"r");
    }
  }
  else {
    strcpy(data,vlst->sol.namepar);
    ptr = strstr(data,".velext");
    if ( !ptr )  strcat(data,".velext");
    in = fopen(data,"r");
  }
  if ( !in ) {
    if ( vlst->info.verb != '0' )  fprintf(stdout," # parameter file %s not found\n",data); 
    return(0);
  }
  if ( vlst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* read flow parameters */
  vlst->sol.nbcl = 0;
  npar = 0;
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for condition type */
    if ( !strcmp(data,"dirichlet")) {
      fscanf(in,"%d",&ncld);
      npar++;
      for (i=vlst->sol.nbcl; i<vlst->sol.nbcl+ncld; i++) {
        pcl = &vlst->sol.cl[i];
        pcl->typ = Dirichlet;
        fscanf(in,"%d %s %c",&pcl->ref,buf,&pcl->att);

        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        pcl->att = tolower(pcl->att);
        if ( !strchr("fv",pcl->att) ) {
          if ( vlst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: [%s] %c\n",buf,pcl->att);
          return(0);
        }
        if ( pcl->att == 'v' ) {
          for (j=0; j<vlst->info.dim; j++)  fscanf(in,"%lf",&pcl->u[j]);
        }
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = VL_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = VL_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = VL_tri;
		  }
      vlst->sol.nbcl += ncld;
    }
    /* Neumann */
    else if ( !strcmp(data,"neumann") ) {
      fscanf(in,"%d",&ncld);
      npar++;
      for (i=vlst->sol.nbcl; i<vlst->sol.nbcl+ncld; i++) {
        pcl = &vlst->sol.cl[i];
        pcl->typ = Neumann;
        fscanf(in,"%d %s %c",&pcl->ref,buf,&pcl->att);
        
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        pcl->att = tolower(pcl->att);
        if ( !strchr("fnv",pcl->att) ) {
          if ( vlst->info.verb != '0' )  fprintf(stdout,"\n # wrong format: [%s] %c\n",buf,pcl->att);
          return(0);
        }
        if ( pcl->att == 'v' ) {
          for (j=0; j<vlst->info.dim; j++)  fscanf(in,"%lf",&pcl->u[j]);
        }
        else if ( pcl->att == 'n' )  fscanf(in,"%lf ",&pcl->u[0]);
        
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = VL_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = VL_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = VL_tri;

        /* for the time being: no normal et vertices known */
        if ( (pcl->elt == VL_ver) && (pcl->att == 'n') ) {
          if ( vlst->info.verb != '0' )  fprintf(stdout,"\n # condition not allowed: [%s] %c\n",buf,pcl->att);
          return(0);
        }
      }
      vlst->sol.nbcl += ncld;
    }
    /* gravity or body force */
    else if ( !strcmp(data,"gravity") ) {
			npar++;
      vlst->sol.cltyp |= Gravity;
      for (j=0; j<vlst->info.dim; j++)
        fscanf(in,"%lf",&vlst->sol.gr[j]);
    }    
    else if ( !strcmp(data,"domain") ) {
      npar++;
      fscanf(in,"%d",&ncld);
      assert(ncld <= VL_MAT);
      vlst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &vlst->sol.mat[i];
        fscanf(in,"%d %lf",&pm->ref,&pm->alpha);
      }
    }
  }
  fclose(in);
  if ( (npar > 0) && (vlst->info.verb != '0') )  fprintf(stdout," %d parameters\n",npar);

  return(1);
}


int main(int argc,char *argv[]) {
  VLst    vlst;
  int     ier;
  char    stim[32];

  tminit(vlst.info.ctim,TIMEMAX);
  chrono(ON,&vlst.info.ctim[0]);

  /* interrupts */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);

  /* init structure */
  memset(&vlst.mesh,0,sizeof(Mesh));
  memset(&vlst.sol,0,sizeof(Sol));
	vlst.sol.cl  = (Cl*)calloc(VL_CL,sizeof(Cl));
	vlst.sol.mat = (Mat*)calloc(VL_MAT,sizeof(Cl));
  vlst.sol.res = VL_RES;
  vlst.sol.nit = VL_MAXIT;
  vlst.sol.nt  = 0;
  vlst.sol.tim = 0.0;
  vlst.sol.alpha = VL_ALPHA;

  /* global parameters */
  vlst.info.dim    = 3;
	vlst.info.ver    = 1;
  vlst.info.verb   = '1';
  vlst.info.zip    = 0;
  vlst.info.ls     = 0;
  vlst.info.mfree  = 1;
  
  /* parse command line */
  if ( !parsar(argc,argv,&vlst) )  return(1);

  /* chrono */
  chrono(ON,&vlst.info.ctim[1]);

  if ( vlst.info.verb != '0' ) {
    fprintf(stdout," - VELEXT, Release %s, %s\n   %s\n\n",VL_VER,VL_REL,VL_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* loading mesh */
  ier = loadMesh(&vlst);
	if ( ier <=0 )  return(1);

  /* parse parameters in file */
  if ( !parsop(&vlst) )  return(1);

  /* allocating memory */
  if ( !vlst.sol.u ) {
    vlst.sol.u = (double*)calloc(vlst.info.dim*vlst.info.npi,sizeof(double));
    assert(vlst.sol.u);
  }

  /* loading solution (or Dirichlet values) */
  if ( vlst.sol.namein ) {
    ier = loadSol(&vlst);
    if ( !ier )  return(1);
  }

  /* load characteristic function */
  if ( vlst.info.ls ) {
    vlst.sol.chi = (double*)calloc(vlst.info.dim*vlst.info.np,sizeof(double));
    assert(vlst.sol.chi);
    ier = loadChi(&vlst);
    if ( ier <= 0 ) {
      if ( vlst.info.verb != '0' )  fprintf(stdout," # missing or wrong file %s",vlst.sol.namechi);
      return(1);
    } 
  }

  /* packing mesh if needed */
  if ( vlst.sol.nmat ) {
    ier = vlst.info.dim == 2 ? pack_2d(&vlst) : pack_3d(&vlst);
    if ( ier == 0 ) {
      if ( vlst.info.verb != '0' )  fprintf(stdout," # Packing error.\n");
      return(1);
    }
  }
  
  /* build adjacency table */
  if ( vlst.info.ls )
    vlst.info.dim == 2 ? hashel_2d(&vlst) : hashel_3d(&vlst);

  chrono(OFF,&vlst.info.ctim[1]);
  printim(vlst.info.ctim[1].gdif,stim);
  if ( vlst.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);

  /* solve */
  chrono(ON,&vlst.info.ctim[2]);
  if ( vlst.info.verb != '0' )
    fprintf(stdout,"\n ** MODULE VELEXT: %s\n",VL_VER);

  ier = VL_velext(&vlst);
  if ( !ier )  return(1);

  chrono(OFF,&vlst.info.ctim[2]);
  if ( vlst.info.verb != '0' ) {
		printim(vlst.info.ctim[2].gdif,stim);
    fprintf(stdout," ** COMPLETED: %s\n\n",stim);
	}

  /* save file */
  if ( vlst.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&vlst.info.ctim[3]);
  if ( vlst.info.zip && !unpack(&vlst) )  return(1);

  if ( !vlst.sol.nameout ) {
    vlst.sol.nameout = (char *)calloc(128,sizeof(char));
    assert(vlst.sol.nameout);
    strcpy(vlst.sol.nameout,vlst.mesh.name);
  }

  ier = saveSol(&vlst);
	if ( !ier )   return(1);
  chrono(OFF,&vlst.info.ctim[3]);
  if ( vlst.info.verb != '0' ) {
    printim(vlst.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
	free(vlst.sol.u);
  if ( vlst.sol.chi )  free(vlst.sol.chi);
	free(vlst.sol.cl);
	free(vlst.sol.mat);

  chrono(OFF,&vlst.info.ctim[0]);
  if ( vlst.info.verb != '0' ) {
    printim(vlst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}




