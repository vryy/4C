/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_call_mat' which selects proper material
law

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief selects proper material law

<pre>                                                              fh 10/02
This routine selects the proper material law for a spatial 1d-beam-element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (o)  actual material
\param *eps      DOUBLE    (i)  strain vector
\param **bop     DOUBLE    (i)  B-Operator matrix
\param **D       DOUBLE    (o)  constitutive matrix
\param  *stress  DOUBLE    (o)  stress at actual gauss point
\param ip        INT       (i)  actual integration point
\param istore    INT       (i)  flag for storing of new stresses
\param newval    INT       (o)  flag for evaluating of new stresses


\warning There is nothing special in this routine
\return void
\sa calling:   b3_mat_linel() , b3_mat_plast_mises()
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_call_mat(ELEMENT   *ele,
                 MATERIAL  *mat,
                 DOUBLE    *eps,
		 DOUBLE   **bop,
                 DOUBLE   **d,
                 DOUBLE    *stress,
		 INT        ip,
		 INT        istore,
		 INT        newval)
{

const INT numdf=6;
INT init=0;

#ifdef DEBUG
dstrc_enter("b3_call_mat");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_stvenant:/*--------------------------------- linear elastic ---*/
    b3_mat_linel(ele,mat->m.stvenant->youngs,
                 mat->m.stvenant->possionratio,d);
    /*---------- calculate internal forces of actual gauss point--------*/
    if (newval==1) math_matvecdense(stress,d,eps,numdf,numdf,0,1.);
  break;

  case m_pl_mises:/*--------------------------- von Mises plasticity ---*/
    b3_mat_plast_mises(mat->m.pl_mises->youngs,
                       mat->m.pl_mises->possionratio,
		       mat->m.pl_mises->ALFAT,
		       mat->m.pl_mises->Sigy,
		       mat->m.pl_mises->Hard,
		       mat->m.pl_mises->GF,
		       ele,eps,ip,stress,d,istore,newval,init);
  break;
  default:
    dserror(" unknown type of material law");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_call_mat */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
