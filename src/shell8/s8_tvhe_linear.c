#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  calculates metrics                                    m.gee 6/01    |
 *----------------------------------------------------------------------*/
/*
C.......................................................................
C!    METRIK DES SCHALENRAUMS 
c!    gmkovc wird hier erneut berechnet, um die Vernachlaessigung der
c!    in e3 quadratischen Anteile zu beruecksichtigen, d.h.
c!    gmkovc_ij ungleich gkovc_i*gkovc_j                                         .
C.......................................................................
*/
void s8_tvhe_linear(DOUBLE **gmkovr,
                    DOUBLE **gmkovc,
                    DOUBLE **gmkonr,
                    DOUBLE **gmkonc,
                    DOUBLE **gkovr,
                    DOUBLE **gkovc,
                    DOUBLE  *detr,
                    DOUBLE  *detc)
{
INT i,j,k;
DOUBLE heps[3][3];
DOUBLE skalar;
DOUBLE det_dummy;
#ifdef DEBUG 
dstrc_enter("s8_tvhe_linear");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      skalar=0.0;
      for (k=0; k<3; k++) skalar += gkovc[k][i]*gkovr[k][j];
      heps[i][j] = skalar;
   }
}
gmkovc[0][0] = 2.0*heps[0][0] - gmkovr[0][0];
gmkovc[1][1] = 2.0*heps[1][1] - gmkovr[1][1];
gmkovc[2][2] = 2.0*heps[2][2] - gmkovr[2][2];

gmkovc[0][1] = heps[0][1]+heps[1][0]-gmkovr[0][1];
gmkovc[0][2] = heps[0][2]+heps[2][0]-gmkovr[0][2];
gmkovc[1][2] = heps[1][2]+heps[2][1]-gmkovr[1][2];

gmkovc[1][0] = gmkovc[0][1];
gmkovc[2][0] = gmkovc[0][2];
gmkovc[2][1] = gmkovc[1][2];

math_array_copy(gmkovr,3,3,gmkonr);
math_inv3(gmkonr,&det_dummy);
if (det_dummy <= 0.0) det_dummy = 1.0E-08;
*detr = sqrt(det_dummy);

math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = 1.0E-08;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tvhe_linear */
#endif
