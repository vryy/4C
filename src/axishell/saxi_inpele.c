/*!----------------------------------------------------------------------
\file
\brief contains the routine 'saxiinp' which reads the axisymmetric 
       shell element

*----------------------------------------------------------------------*/
#ifdef D_AXISHELL
#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"

/*! 
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief reads an axisymmetric shell element from the input file

<pre>                                                              mn 05/03 
This routine reads an axisymmetric shell element from the input file

</pre>
\param *ele  ELEMENT  (o)   the element

\warning buffer[50] is not needed locally
\return void                                               
\sa caling:    ---; 
    called by: inp_struct_field()

*----------------------------------------------------------------------*/
void saxi_inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("saxi_inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.saxi = (AXISHELL*)CCACALLOC(1,sizeof(AXISHELL));
if (ele->e.saxi==NULL) dserror("Allocation of element failed");
/*---------------------------------------------- read elements topology */
frchk("LINE2",&ierr);
if (ierr==1) 
{
   ele->distyp = line2;
   ele->numnp=2;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("LINE2",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("LINE3",&ierr);
if (ierr==1) 
{
   ele->distyp = line3;
   ele->numnp=3;
   ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("LINE3",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of AXISHELL element failed");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of saxi_inp */
/*----------------------------------------------------------------------*/
#endif /*D_AXISHELL*/
/*! @} (documentation module close)*/
