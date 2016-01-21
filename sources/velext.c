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
  fprintf(stdout,"usage: %s [+/-v | -h | -c] [-a val] [-e err] [-n nit] data_mesh[.mesh] [-s data_vel[.sol]]\n",prog);
  exit(1);
}


/* parse command line */
static int parsar(int argc,char *argv[],VLst *vlst) {
  int      i;
  char    *ptr;
  
  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i] == '+') ) {
      switch(argv[i][1]) {
      case 'h':
      case '?':
        usage(argv[0]);
        break;
			case 'a':
        if ( !strcmp(argv[i],"-a") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            vlst->sol.alpha = atof(argv[i]);
          else
            i--;
        }
				break;
      case 'c':
		    if ( !strcmp(argv[i],"-c") ) {
		      vlst->info.ls = 1;
		    }
        break;  
      case 'e':
        if ( !strcmp(argv[i],"-e") ) {
          if ( ++i && isdigit(argv[i][0]) )
            vlst->sol.err = strtod(argv[i],NULL);
          else
            --i; 
        }
        break;
      case 'i':
        if ( !strcmp(argv[i],"-i") ) {
          ++i;
          vlst->mesh.name = argv[i];
        }
        break;
      case 'n':
        if ( !strcmp(argv[i],"-n") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            vlst->sol.nit = atoi(argv[i]);
          else
            i--;
        }
        break;
      case 's':
			  if ( !strcmp(argv[i],"-s") ) {
			    ++i;
			    vlst->sol.namein = argv[i];
			    ptr = strstr(vlst->sol.namein,".sol");
			    if ( !ptr )  strcat(vlst->sol.namein,".sol");
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
      if ( vlst->mesh.name == NULL )
        vlst->mesh.name = argv[i];
      else if ( vlst->sol.namein == NULL ) {
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
  float       fp1;
  int         i,j,npar,ncld,ret;
  char       *ptr,buf[256],data[256];
  FILE       *in;

  /* check for parameter file */
  strcpy(data,vlst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';
  strcat(data,".velext");
  in = fopen(data,"r");
  if ( !in ) {
    sprintf(data,"%s","DEFAULT.velext");
    in = fopen(data,"r");
    if ( !in )  return(1);
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
    if ( !strcmp(data,"dirichlet") || !strcmp(data,"neumann") ) {
      fscanf(in,"%d",&ncld);
      npar++;
      for (i=vlst->sol.nbcl; i<vlst->sol.nbcl+ncld; i++) {
        pcl = &vlst->sol.cl[i];
        if ( !strcmp(data,"dirichlet") )     pcl->typ = Dirichlet; 
        else if ( !strcmp(data,"neumann") )  pcl->typ = Neumann;

        /* check for entity */
        fscanf(in,"%d %s ",&pcl->ref,buf);
        for (j=0; j<strlen(buf); j++)  buf[j] = tolower(buf[j]);
        fscanf(in,"%c",&pcl->att);
        pcl->att = tolower(pcl->att);
        if ( (pcl->typ == Dirichlet) && (pcl->att != 'v' && pcl->att != 'f') ) {
          fprintf(stdout,"\n # wrong format: %s %c\n",buf,pcl->att);
          continue;
        }
        else if ( (pcl->typ == Load) && (pcl->att != 'v' && pcl->att != 'f' && pcl->att != 'n') ) {
          fprintf(stdout,"\n # wrong format: %s %c\n",buf,pcl->att);
          continue;
        }
        if ( !strcmp(buf,"vertices") || !strcmp(buf,"vertex") )          pcl->elt = VL_ver;
        else if ( !strcmp(buf,"edges") || !strcmp(buf,"edge") )          pcl->elt = VL_edg;
        else if ( !strcmp(buf,"triangles") || !strcmp(buf,"triangle") )  pcl->elt = VL_tri;
        else {
          fprintf(stdout,"  %%%% Wrong format: %s\n",buf);
          continue;
        }

        if ( pcl->att == 'v' ) {
          for (j=0; j<vlst->info.dim; j++) {
            fscanf(in,"%f ",&fp1);
            pcl->u[j] = fp1;
          }
        }
        else if ( pcl->att == 'n' ) {
          fscanf(in,"%f ",&fp1);
          pcl->u[0] = fp1;
        }
		  }
      vlst->sol.nbcl += ncld;
    }
    else if ( !strcmp(data,"domain") ) {
      npar++;
      fscanf(in,"%d",&ncld);
      vlst->sol.nmat = ncld;
      for (i=0; i<ncld; i++) {
        pm = &vlst->sol.mat[i];
        fscanf(in,"%d %lf",&pm->ref,&pm->alpha);
      }
    }
  }
  fclose(in);

  /* identify type of BC */
  for (i=0; i<vlst->sol.nbcl; i++) {
    pcl = &vlst->sol.cl[i];
    vlst->sol.cltyp |= pcl->elt;
  }

  if ( npar > 0 && vlst->info.verb != '0' ) {
    fprintf(stdout," %d conditions\n",npar);
  }

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
  vlst.sol.err = VL_RES;
  vlst.sol.nit = VL_MAXIT;

  /* global parameters */
  vlst.info.dim    = 3;
	vlst.info.ver    = 1;
  vlst.info.verb   = '1';
  vlst.info.zip    = 0;
  vlst.info.mfree  = 1;
  vlst.info.ls     = 0;
  
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
    if ( !ier )  return(1);
    else if ( ier < 0 ) {
      free(vlst.sol.chi);
      vlst.info.ls = 0;
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

  chrono(OFF,&vlst.info.ctim[0]);
  if ( vlst.info.verb != '0' ) {
	  printim(vlst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}




