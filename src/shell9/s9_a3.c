/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9a3:  which evaluates the directors at nodal points, and writes the
          normed director on 'ele->e.s9->a3ref'
 - s9averdir: which evaluates the director 'Bischoff style', which means
              the director is averaged if the shell has kinks. The shared
              director is later put back on 'ele->e.s9->a3ref'
 - s9a3ref_extern: which evaluates the directors at nodal points, norms
                   them, but these dirctors are only localy used for the
                   calculation of the element surface loads ('s9_load1').
                   Here no averaged director is used!


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief evaluates the director at nodal points

<pre>                     m.gee 6/01              modified by    sh 02/03
This routine evaluates the director at the nodal points of a shell9
element. The directors are normed to unit length an writen to the element
array 'ele->e.s9->a3ref'.
</pre>
\param  ELEMENT   *ele   (i->modified for a3ref) actual shell9 element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9init() [s9_init.c]

*----------------------------------------------------------------------*/
void s9a3(ELEMENT   *ele)
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9a3");
#endif
/*----------------------------------------------------------------------*/
iel = ele->numnp;
funct = amdef("funct",&funct_a,iel,1,"DV");
deriv = amdef("deriv",&deriv_a,2,iel,"DA");
a3ref = amdef("a3ref",&(ele->e.s9->a3ref),3,iel,"DA");
thick = amdef("thick",&(ele->e.s9->thick_node),iel,1,"DV");
for (i=0; i<iel; i++) thick[i] =  ele->e.s9->thick;
/*--------------------------------------------------- loop nodal points */
for (i=0; i<iel; i++)
{
   r = s9_local_coord_node(i,0,ele->distyp);
   s = s9_local_coord_node(i,1,ele->distyp);
   s9_funct_deriv(funct,deriv,r,s,ele->distyp,1);
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
} /* end of s9a3 */


/*!----------------------------------------------------------------------
\brief evaluates the director at nodal points

<pre>              m.gee 10/01             modified by           sh 02/03
This routine evaluates the director at the nodal points of a shell9
element. The directors are normed to unit length, but they are not writen
on the element array but only used for the calculation of the element
loads. Therefor no shared director is used.
NOTE: The routine is identical to 's9a3' except for the storage of a3ref!!
</pre>
\param  DOUBLE   **a3ref (0) director at nodal points
\param  ELEMENT   *ele   (i) actual shell9 element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9eleload() [s9_load1.c]

*----------------------------------------------------------------------*/
void s9a3ref_extern(DOUBLE   *funct,
                    DOUBLE  **deriv,
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9a3ref_extern");
#endif
/*----------------------------------------------------------------------*/
iel = ele->numnp;
/*--------------------------------------------------- loop nodal points */
for (i=0; i<iel; i++)
{
   r = s9_local_coord_node(i,0,ele->distyp);
   s = s9_local_coord_node(i,1,ele->distyp);
   s9_funct_deriv(funct,deriv,r,s,ele->distyp,1);
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
} /* end of s9a3ref_extern */




/*!----------------------------------------------------------------------
\brief makes the modified director bischoff style

<pre>              m.gee 6/01              modified by           sh 02/03
This routine evaluates the shared director at kink of adjacent shell9
elements. The modification is done according to Dis. Bischoff (p.129,
picture 8.2c). The shared director is later put back on the element array
'ele->e.s9->a3ref'
</pre>
\param  DOUBLE **dir_list (i) list of different directors at this point
\param  INT      numa3    (i) number of directors at this point
\param  DOUBLE  *a3       (o) shared director

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9init() [s9_init.c]

*----------------------------------------------------------------------*/
void s9averdir(DOUBLE **dir_list, INT numa3, DOUBLE *a3)
{
INT        i;
DOUBLE     averdir[3];
DOUBLE     davn[3];
DOUBLE     normal[3];
DOUBLE     lenght;
DOUBLE     denom;
DOUBLE     alpha;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9averdir");
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
   lenght = DSQR(normal[0])+DSQR(normal[1])+DSQR(normal[2]);
   if (lenght <= EPS12)/*----------------------- directors are parallel */
   {
      davn[0] = 0.5*(averdir[0]+dir_list[0][i]);
      davn[1] = 0.5*(averdir[1]+dir_list[1][i]);
      davn[2] = 0.5*(averdir[2]+dir_list[2][i]);
   }
   else  /*----------------------- averaging of non-parallel directors */
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
/*  make normation of director in s9_init.c -> elements could have different thicknesses !!!*/
/*   a3[0] = davn[0]/h2;*/
/*   a3[1] = davn[1]/h2;*/
/*   a3[2] = davn[2]/h2;*/
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9averdir */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
