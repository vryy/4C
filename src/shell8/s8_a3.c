/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"

/*----------------------------------------------------------------------*
 | evaluation of directors at nodal points                m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8a3(ELEMENT   *ele,
          S8_DATA   *data,/* is' hier ueberfluessig, rausschmeissen!*/
          INT        option)
{
INT        i,j;
INT        ialpha;
INT        idim;
INT        inode;
INT        iel;
DOUBLE     r,s;
DOUBLE     gkov[3][3];
DOUBLE     a3[3];
DOUBLE     a3norm;
ARRAY      funct_a; DOUBLE  *funct;
ARRAY      deriv_a; DOUBLE **deriv;
DOUBLE   **a3ref;
DOUBLE    *thick;

#ifdef DEBUG
dstrc_enter("s8a3");
#endif
/*----------------------------------------------------------------------*/
iel = ele->numnp;
funct = amdef("funct",&funct_a,iel,1,"DV");
deriv = amdef("deriv",&deriv_a,2,iel,"DA");
a3ref = amdef("a3ref",&(ele->e.s8->a3ref),3,iel,"DA");
thick = amdef("thick",&(ele->e.s8->thick_node),iel,1,"DV");
for (i=0; i<iel; i++) thick[i] =  ele->e.s8->thick;
/*--------------------------------------------------- loop nodal points */
for (i=0; i<iel; i++)
{
   r = s8_local_coord_node(i,0,ele->distyp);
   s = s8_local_coord_node(i,1,ele->distyp);
   s8_funct_deriv(funct,deriv,r,s,ele->distyp,1);
/*-------------------------------------------------------------- a1, a2 */
   for (ialpha=0; ialpha<2; ialpha++)
   {
      for (idim=0; idim<3; idim++)
      {
         gkov[idim][ialpha]=0.0;
         for (inode=0; inode<iel; inode++)
         {
            gkov[idim][ialpha] +=

            deriv[ialpha][inode] * ele->node[inode]->x[idim];
         }
      }
   }
/*------------------------------------------------------------------ a3 */
a3[0] = gkov[1][0]*gkov[2][1] - gkov[2][0]*gkov[1][1];
a3[1] = gkov[2][0]*gkov[0][1] - gkov[0][0]*gkov[2][1];
a3[2] = gkov[0][0]*gkov[1][1] - gkov[1][0]*gkov[0][1];
a3norm = a3[0]*a3[0] + a3[1]*a3[1] + a3[2]*a3[2];
a3norm = sqrt(a3norm);
a3norm = 1.0/a3norm;
a3[0] *= a3norm;
a3[1] *= a3norm;
a3[2] *= a3norm;
for (j=0; j<3; j++) a3ref[j][i] = a3[j];
} /* end of loop over nodes */
/*----------------------------------------------------------------------*/
amdel(&funct_a);
amdel(&deriv_a);
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8a3 */


/*----------------------------------------------------------------------*
 | evaluation of directors at nodal points                m.gee 10/01   |
 *----------------------------------------------------------------------*/
void s8a3ref_extern(DOUBLE   *funct,
                       DOUBLE  **deriv,
                       DOUBLE   *thick,
                       DOUBLE  **a3ref,
                       ELEMENT  *ele)
{
INT        i,j;
INT        ialpha;
INT        idim;
INT        inode;
INT        iel;
DOUBLE     r,s;
DOUBLE     gkov[3][3];
DOUBLE     a3[3];
DOUBLE     a3norm;

#ifdef DEBUG
dstrc_enter("s8a3ref_extern");
#endif
/*----------------------------------------------------------------------*/
iel = ele->numnp;
/*--------------------------------------------------- loop nodal points */
for (i=0; i<iel; i++)
{
   r = s8_local_coord_node(i,0,ele->distyp);
   s = s8_local_coord_node(i,1,ele->distyp);
   s8_funct_deriv(funct,deriv,r,s,ele->distyp,1);
/*-------------------------------------------------------------- a1, a2 */
   for (ialpha=0; ialpha<2; ialpha++)
   {
      for (idim=0; idim<3; idim++)
      {
         gkov[idim][ialpha]=0.0;
         for (inode=0; inode<iel; inode++)
         {
            gkov[idim][ialpha] +=

            deriv[ialpha][inode] * ele->node[inode]->x[idim];
         }
      }
   }
/*------------------------------------------------------------------ a3 */
a3[0] = gkov[1][0]*gkov[2][1] - gkov[2][0]*gkov[1][1];
a3[1] = gkov[2][0]*gkov[0][1] - gkov[0][0]*gkov[2][1];
a3[2] = gkov[0][0]*gkov[1][1] - gkov[1][0]*gkov[0][1];
a3norm = a3[0]*a3[0] + a3[1]*a3[1] + a3[2]*a3[2];
a3norm = sqrt(a3norm);
a3norm = 1.0/a3norm;
a3[0] *= a3norm;
a3[1] *= a3norm;
a3[2] *= a3norm;
for (j=0; j<3; j++) a3ref[j][i] = a3[j];
} /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8a3ref_extern */




/*----------------------------------------------------------------------*
 | modified director bischoff style                       m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8averdir(DOUBLE **dir_list, INT numa3, DOUBLE *a3)
{
INT        i;
DOUBLE     averdir[3];
DOUBLE     davn[3];
DOUBLE     normal[3];
DOUBLE     lenght;
DOUBLE     denom;
DOUBLE     alpha;
#ifdef DEBUG
dstrc_enter("s8averdir");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- now loop over all directors */
averdir[0] = dir_list[0][0];
averdir[1] = dir_list[1][0];
averdir[2] = dir_list[2][0];
for (i=1; i<numa3; i++)
{
   /*------------------------------ make cross product of two directors */
   normal[0] = averdir[1]*dir_list[2][i] - averdir[2]*dir_list[1][i];
   normal[1] = averdir[2]*dir_list[0][i] - averdir[0]*dir_list[2][i];
   normal[2] = averdir[0]*dir_list[1][i] - averdir[1]*dir_list[0][i];
   lenght = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
   if (lenght <= EPS12)/*----------------------- directors are parallel */
   {
      davn[0] = 0.5*(averdir[0]+dir_list[0][i]);
      davn[1] = 0.5*(averdir[1]+dir_list[1][i]);
      davn[2] = 0.5*(averdir[2]+dir_list[2][i]);
   }
   else
   {
      denom =
                  (DSQR(dir_list[0][i])+DSQR(dir_list[2][i]))*DSQR(averdir[1])
                 +(-2.*dir_list[0][i]*averdir[0]*dir_list[1][i]-2.*dir_list[2][i]
                 *averdir[2]*dir_list[1][i])*averdir[1]+(DSQR(dir_list[2][i])
                 +DSQR(dir_list[1][i]))*DSQR(averdir[0])-2.*averdir[2]*averdir[0]
                 *dir_list[2][i]*dir_list[0][i]+(DSQR(dir_list[0][i])+DSQR(dir_list[1][i]))
                 *DSQR(averdir[2]);

      if (ABS(denom)<=EPS13) dserror("Making of mod. directors failed");

      alpha   =  (averdir[2]*dir_list[2][i]-DSQR(dir_list[0][i])+averdir[0]*dir_list[0][i]
                  -DSQR(dir_list[1][i])+dir_list[1][i]*averdir[1]-DSQR(dir_list[2][i]))/denom;

      davn[0] =-alpha*DSQR(averdir[1])*dir_list[0][i]+alpha*averdir[1]*averdir[0]*dir_list[1][i]
               +averdir[0]+alpha*averdir[2]*averdir[0]*dir_list[2][i]-alpha*DSQR(averdir[2])*dir_list[0][i];

      davn[1] =alpha*averdir[0]*averdir[1]*dir_list[0][i]+averdir[1]+alpha*averdir[2]*averdir[1]
               *dir_list[2][i]-alpha*DSQR(averdir[0])*dir_list[1][i]-alpha*DSQR(averdir[2])*dir_list[1][i];

      davn[2] =-alpha*DSQR(averdir[1])*dir_list[2][i]+alpha*averdir[1]*averdir[2]*dir_list[1][i]
               -alpha*DSQR(averdir[0])*dir_list[2][i]+alpha*averdir[0]*averdir[2]*dir_list[0][i]+averdir[2];
   }
   a3[0] = davn[0];
   a3[1] = davn[1];
   a3[2] = davn[2];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8averdir */
#endif
