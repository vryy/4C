/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_tvmr: which calculates all the metrics in the middle surface of each 
            kinematic layer -> akov, akon, amkov, amkon, ...


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
\brief calculates the metrics in the middle surface of each kinematic layer                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the metrics (akov, akon, amkov, amkon) in the 
middle of each kinematic layer.
</pre>
\param  DOUBLE   **x       (i)  coordinates at nodal points
\param  DOUBLE  ***a3      (i)  director of each kinematic layer
\param  DOUBLE  ***akov    (o)  kovariant basis vectors in reference layer of each kinematic layer
\param  DOUBLE  ***akon    (o)  kontravariant basis vectors in reference layer of each kinematic layer
\param  DOUBLE  ***amkov   (o)  kovariant metric in reference layer of each kinematic layer
\param  DOUBLE  ***amkon   (o)  kontravariant metric in reference layer of each kinematic layer
\param  DOUBLE    *det     (o)  det of akon
\param  DOUBLE    *funct   (i)  shape functions at GP
\param  DOUBLE   **deriv   (i)  shape function derivatives at GP
\param  INT        iel     (i)  number of nodes to this element
\param  DOUBLE  ***a3kvp   (o)  partial derivatives of a3_L for each kinematic layer
\param  INT        num_klay(i)  number of kin layers to this element  

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug()       [s9_static_keug.c]
                             s9_stress()           [s9_stress.c]
                             s9_ans_colloqpoints() [s9_ans.c]
                             s9jaco()              [s9_jaco.c]

*----------------------------------------------------------------------*/
void s9_tvmr(DOUBLE    **x,
             DOUBLE   ***a3,
             DOUBLE   ***akov,
             DOUBLE   ***akon,
             DOUBLE   ***amkov,
             DOUBLE   ***amkon,
             DOUBLE    **akovh,
             DOUBLE    **akonh,
             DOUBLE    **amkovh,
             DOUBLE    **amkonh,
             DOUBLE     *det,
             DOUBLE     *funct,
             DOUBLE    **deriv,
             INT         iel,
             DOUBLE   ***a3kvp,
             INT         num_klay)
{
INT    i,j,k,klay,idim,ialpha,inode;
DOUBLE det_dummy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_tvmr");
#endif
/*----------------------------------------------------------------------*/
/*-------------------- loop over all layers ----------------------------*/
for (klay=0; klay<num_klay; klay++) 
{  
 /*------------------------ interpolation --------- a1,a2 (kov.) ------*/
  for (idim=0; idim<3; idim++)
  {
    akov[idim][0][klay] = 0.0;
    akov[idim][1][klay] = 0.0;
    for (inode=0; inode<iel; inode++)
    {
      akov[idim][0][klay] += deriv[0][inode] * x[idim][inode];
      akov[idim][1][klay] += deriv[1][inode] * x[idim][inode];
    }
  }
 /*------------------------ interpolation ------------------ a3 ------*/
  for (idim=0; idim<3; idim++)
  {
    akov[idim][2][klay] = 0.0;
    for (inode=0; inode<iel; inode++)
    {
      akov[idim][2][klay] += funct[inode] * a3[idim][inode][klay];
    }
  }
 /*--------------- kontravariant basis vectors g1,g2,g3 (inverse of kov) */
    /* copy akov[3][3][klay] -> akovh[3][3] */
  for (i=0; i<3; i++) for (j=0; j<3; j++) akovh[i][j] = akov[i][j][klay];

    math_array_copy(akovh,3,3,akonh);
    math_inv3(akonh,det);
    math_tran(akonh,3);

    /* copy akonh[3][3] -> akon[3][3][klay] */
  for (i=0; i<3; i++) for (j=0; j<3; j++) akon[i][j][klay] = akonh[i][j];
 /*--------------------------------------------- kovariant metrik tensor */
  for (i=0; i<3; i++)
  {
     for (j=i; j<3; j++)
     {
        amkov[i][j][klay]=0.0;
        for (k=0; k<3; k++) 
        amkov[i][j][klay] += akov[k][i][klay]*akov[k][j][klay];
     }
  }   
        amkov[1][0][klay] = amkov[0][1][klay];
        amkov[2][0][klay] = amkov[0][2][klay];
        amkov[2][1][klay] = amkov[1][2][klay];
 /*----------------------------------------- kontravariant metrik tensor */
    /* copy amkov[3][3][klay] -> amkovh[3][3] */
  for (i=0; i<3; i++) for (j=0; j<3; j++) amkovh[i][j] = amkov[i][j][klay];

    math_array_copy(amkovh,3,3,amkonh);
    math_inv3(amkonh,&det_dummy);

    /* copy amkonh[3][3] -> amkon[3][3][klay] */
  for (i=0; i<3; i++) for (j=0; j<3; j++) amkon[i][j][klay] = amkonh[i][j];
 /*------------------------------------------- partial derivatives of a3 */
  for (ialpha=0; ialpha<2; ialpha++)
  {
     for (idim=0; idim<3; idim++)
     {
        a3kvp[idim][ialpha][klay]=0.0;
        for (inode=0; inode<iel; inode++)
        {
          a3kvp[idim][ialpha][klay] += deriv[ialpha][inode] * a3[idim][inode][klay];
        }
     }
  }

} /*---------------- end loop over all layers ---------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tvmr */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
