/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_intforce: which calculates the internal forces: 
                R = int (Bt * stress_r)

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief make internal forces                                      

<pre>                     m.gee 11/01             modified by    sh 11/02
This routine calculates the internal forces:
      R = int (Bt * stress_r)
</pre>
\param  double *intforce  (i/o) internal forces to be modified (intforce +=)
\param  double *stress_r   (i)  stress resultants (integrated over thickness)
\param  double **bop       (i)  B-Operator matrix
\param  int      iel       (i)  number of nodes to this element
\param  int      numdf     (i)  number of degrees of freedom to one node
\param  int      nstress_r (i)  number of stress resultants (shell9 -> 12)
\param  double   weight    (i)  weight at GP

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_intforce(double *intforce, double *stress_r, double **bop,
                 int iel, int numdf, int nstress_r, double weight)
{
int nd;

#ifdef DEBUG 
dstrc_enter("s9_intforce");
#endif
/*----------------------------------------------------------------------*/
nd = iel*numdf;
/*
make intforce[nd] = transposed(bop[nstress_r][nd]) * stress_r[nstress_r]
*/
math_mattrnvecdense(intforce,bop,stress_r,nd,nstress_r,1,weight);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_intforce */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
