/*!----------------------------------------------------------------------
\file
\brief ale control part of fsi-problems 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
#include "fsi_prototypes.h"                   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;  
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief  solving for mesh displacements 

<pre>                                                          ck 06/03
control of fsi ale part; program continues depending on ALE_TYP in input

</pre>

\param  *actfield   FIELD          (i)     ale field 			  
\param   mctrl      INT            (i)     control flag		  
\param   numfa      INT            (i)     number of ale field	  
\warning 
\return void                                               

\sa   calling: fsi_ale_lin(), fsi_ale_nln(), fsi_ale_2step(), 
               fsi_ale_spring(), fsi_ale_laplace()
      called by: fluid_mf()

*----------------------------------------------------------------------*/
void fsi_ale(
               FIELD            *actfield,
               INT               mctrl
	    )
{
#ifdef D_FSI
ALE_DYNAMIC *adyn;

adyn = alldyn[genprob.numaf].adyn;

#ifdef DEBUG 
dstrc_enter("fsi_ale");
#endif
switch (adyn->typ)
{
/*---------------------------------------- purely linear calculation ---*/
   case classic_lin:
      fsi_ale_lin(actfield,mctrl);
   break;
/*------------------- incremental calculation stiffened with min J_e ---*/
   case min_Je_stiff:
      fsi_ale_nln(actfield,mctrl);
   break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in 
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
   case two_step:
      fsi_ale_2step(actfield,mctrl);
   break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for 
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
   case springs:
       fsi_ale_spring(actfield,mctrl);
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
   case laplace:
       fsi_ale_laplace(actfield,mctrl);
    break;
/*-------------------------------------------- Large Amplitude Sloshing */
/*  just a interpolation for large amplitude sloshing                    */
   case LAS:
      fsi_ale_LAS(actfield,mctrl);
   break;
/*---------------------------------------------------------- default ---*/
   default:
      dserror("unknown ale typ for fsi");
   break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
#endif
return;
} /* end of fsi_ale */
#endif
/*! @} (documentation module close)*/
