#ifndef __VE_CALLS_H
#define __VE_CALLS_H

enum {Dirichlet=1, Neumann, Load};
enum {VL_ver=1,VL_edg=2,VL_tri=4,VL_tet=8};


/* data structure */
typedef struct _VLst VLst;

/* prototypes */
VLst *VL_init(int dim, int ver, char typ,char mfree);
int   VL_stop(VLst *vlst);

int   VL_velext(VLst *vlst);


#endif