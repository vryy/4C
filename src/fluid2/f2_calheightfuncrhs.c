/*!----------------------------------------------------------------------
\file
\brief time RHS for fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
/*!--------------------------------------------------------------------- 
\brief galerkin part of time forces for vertical height function

<pre>                                                         genk 05/03

In this routine the galerkin part of the time forces for vel dofs
is calculated:

      /
   + |  psi * phi     d_gamma_FS
    /  

                   /
   - (1-THETA)*dt |  psi * Ux * phi,x     d_gamma_FS
                 /  

                   /
   + (1-THETA)*dt |  psi * Uy      d_gamma_FS
                 /  
   
</pre>
\param   *eforce      DOUBLE           (i/o)  element force vector
\param   *velint      DOUBLE           (i)    vel. at integr. point
\param    phint       DOUBLE           (i)    height func at integr. p.
\param	 *funct       DOUBLE           (i)    nat. shape functions      
\param	  phiderxy    DOUBLE           (i)    global derivative of height func
\param	  fac	      DOUBLE           (i)    weighting factor
\param	  iel	      INT              (i)    num. of nodes in ele
\param	  ngnode      INT              (i)    num. of nodes on edge
\param   *iedgnod     INT              (i)    edge node numbers
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calrhs_vhf(
                     DOUBLE          *eforce,    
		     DOUBLE          *velint,  /* at n+1 */
		     DOUBLE          *vel2int, /* at n */
		     DOUBLE           phiint,    
	   	     DOUBLE          *funct,    
                     DOUBLE         **derxy,
		     DOUBLE           phiderxng,  /* at n+1 */ 
		     DOUBLE           phiderxn,   /* at n */
		     DOUBLE           fac,      
		     INT              iel,
		     INT              ngnode,
		     INT             *iedgnod       
                    ) 
{
INT    irn,irow;	      
INT    nd;                 
DOUBLE c,facr,facl;
DOUBLE tau_supg,tau_dc;
FLUID_DYNAMIC   *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_calrhs_vhf");
#endif

fdyn = alldyn[genprob.numff].fdyn;

facr = fac * fdyn->thsr;
facl = fac * fdyn->thsl;	      
nd = NUMDOF_FLUID2*iel;

/*----------------------------------------------------------------------*
   Calculate iteration force vector:
               /
   + THETA*dt |  psi * Ux * phi,x     d_gamma_FS
             /  
 *----------------------------------------------------------------------*/ 
for(irn=0;irn<ngnode;irn++)
{
   irow = nd+iedgnod[irn];
   eforce[irow]	 += funct[irn]*velint[0]*phiderxng*facl;
} /* end loop over irow */


/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
      /
   + |  psi * phi     d_gamma_FS
    /  
 *----------------------------------------------------------------------*/ 
for(irn=0;irn<ngnode;irn++)
{
   irow = nd+iedgnod[irn];
   eforce[irow]	 += funct[irn]*phiint*fac;
} /* end loop over irow */


/*----------------------------------------------------------------------*
   Calculate time force vector:
                   /
   - (1-THETA)*dt |  psi * Ux_old * phi,x     d_gamma_FS
                 /  
 *----------------------------------------------------------------------*/ 
for(irn=0;irn<ngnode;irn++)
{
   irow = nd+iedgnod[irn];
   eforce[irow]	 -= funct[irn]*vel2int[0]*phiderxn*facr;
} /* end loop over irow */
	      
/*----------------------------------------------------------------------*
   Calculate time force vector:
                   /
   + (1-THETA)*dt |  psi * Uy_old      d_gamma_FS
                 /  
 *----------------------------------------------------------------------*/ 
for(irn=0;irn<ngnode;irn++)
{
   irow = nd+iedgnod[irn];
   eforce[irow]	 += funct[irn]*vel2int[1]*facr;
} /* end loop over irow */

/*------------------------------------------------------- STABILISATION */
if (fdyn->hf_stab>0)
{
   tau_supg = fdyn->tau[3];
   tau_dc   = fdyn->tau[4];   

/*----------------------------------------------------------------------*
   Calculate iteration force vector:
               /
   + THETA*dt |  tau_SUPG * Ux * psi,x * Ux * phi,x     d_gamma_FS
             /  
 *----------------------------------------------------------------------*/ 
   c = fac*tau_supg*velint[0]*velint[0]*phiderxng;
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 += funct[irn]*facl;
   } /* end loop over irow */

/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /
 *----------------------------------------------------------------------*/ 
   c = fac*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 += derxy[0][irn]*phiint*c;
   } /* end loop over irow */

/*----------------------------------------------------------------------*
   Calculate time force vector:
   (this part exists twice, because of the linearisation)
                 /
 -(1-THETA)*dt  |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS   
                /
 *----------------------------------------------------------------------*/ 
   c=TWO*facr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 -= derxy[0][irn]*vel2int[0]*phiderxn*c;
   } /* end loop over irow */	      

/*----------------------------------------------------------------------*
   Calculate time force vector:
                  /
 +(1-THETA)*dt   |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS   
                /
 *----------------------------------------------------------------------*/ 
   c=facr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 += derxy[0][irn]*vel2int[1]*c;
   } /* end loop over irow */ 
} /* endif fdyn->hf_stab>0 */

/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgaltf_vhf */


/*!--------------------------------------------------------------------- 
\brief stabilisation part of time forces for vertical height function

<pre>                                                         genk 05/03

In this routine the galerkin part of the time forces for vel dofs
is calculated:

    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /

                 /
 -(1-THETA)*dt  |  tau_DC * psi,x * phi,x    d_gamma_FS
               /

                 /
 -(1-THETA)*dt  |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS   
                /

                 /
 +(1-THETA)*dt  |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS   
                /
   
</pre>
\param   *eforce      DOUBLE           (i/o)  element force vector
\param   *velint      DOUBLE           (i)    vel. at integr. point
\param    phint       DOUBLE           (i)    height func at integr. p.
\param	 *funct       DOUBLE           (i)    nat. shape functions      
\param	  phiderxy    DOUBLE           (i)    global derivative of height func
\param	  fac	      DOUBLE           (i)    weighting factor
\param	  iel	      INT              (i)    num. of nodes in ele
\param	  ngnode      INT              (i)    num. of nodes on edge
\param   *iedgnod     INT              (i)    edge node numbers
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calstabtf_vhf(
                     ELEMENT         *ele,
                     DOUBLE          *eforce,    
		     DOUBLE          *velint,
		     DOUBLE           phiint,    
	   	     DOUBLE          *funct,    
                     DOUBLE         **derxy,
		     DOUBLE           phiderx,   
		     DOUBLE           fac,      
		     INT              iel,
		     INT              ngnode,
		     INT             *iedgnod       
                    ) 
{
INT    irn,irow;	      
INT    nd;                 
DOUBLE facr;
DOUBLE c;
DOUBLE tau_dc, tau_supg;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_calstabtf_vhf");
#endif

fdyn = alldyn[genprob.numff].fdyn;

facr = fac * fdyn->thsr;	      
nd = NUMDOF_FLUID2*iel;

#if 1
if (iel==4)
{
   dserror("heightfunc stabilisation");
   tau_supg = fdyn->tau[4];
/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
    /
   |  tau_SUPG * Ux * psi,x * phi    d_gamma_FS
  /
 *----------------------------------------------------------------------*/ 
   c = fac*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 += derxy[0][irn]*phiint*c;
   } /* end loop over irow */

/*----------------------------------------------------------------------*
   Calculate time force vector:
                 /
 -(1-THETA)*dt  |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS   
                /
 *----------------------------------------------------------------------*/ 
   c=facr*tau_supg*velint[0]*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 -= derxy[0][irn]*phiderx*c;
   } /* end loop over irow */	      
   
/*----------------------------------------------------------------------*
   Calculate time force vector:
                  /
 +(1-THETA)*dt   |   tau_SUPG * Ux * psi,x * Uy   d_gamma_FS   
                /
 *----------------------------------------------------------------------*/ 
   c=facr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 += derxy[0][irn]*velint[1]*c;
   } /* end loop over irow */ 
} /* endif (ele->e.f2->ifssu!=0) */

if (iel==4)
{
   dserror("heightfunc stabilisation");
   tau_dc = fdyn->tau[3];
/*----------------------------------------------------------------------*
   Calculate time force vector:
                 /
 -(1-THETA)*dt  |  tau_DC * psi,x * phi,x    d_gamma_FS
               /
 *----------------------------------------------------------------------*/ 
   c=facr*tau_dc;
   for(irn=0;irn<ngnode;irn++)
   {
      irow = nd+iedgnod[irn];
      eforce[irow]	 -= derxy[0][irn]*phiderx*c;
   } /* end loop over irow */
}/*endif (ele->e.f2->ifsdc!=0) */	      

#endif
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calstabtf_vhf */

/*!--------------------------------------------------------------------- 
\brief galerkin part of iteration forces for vertical height function

<pre>                                                         genk 05/03

In this routine the galerkin part of the time forces for vel dofs
is calculated:

               /
     THETA*dt |  psi * Ux * phi,x     d_gamma_FS
             /  

   
</pre>
\param   *eforce      DOUBLE           (i/o)  element force vector
\param	 *funct       DOUBLE           (i)    nat. shape functions      
\param   *velint      DOUBLE           (i)    vel. at integr. point
\param	  phiderx     DOUBLE           (i)    global derivative of height func
\param	  fac	      DOUBLE           (i)    weighting factor
\param	  iel	      INT              (i)    num. of nodes in ele
\param	  ngnode      INT              (i)    num. of nodes on edge
\param   *iedgnod     INT              (i)    edge node numbers
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calgalif_vhf(
                     DOUBLE          *eforce,    
	   	     DOUBLE          *funct,    
		     DOUBLE          *velint,    
		     DOUBLE           phiderx,   
		     DOUBLE           fac,      
		     INT              iel,
		     INT              ngnode,
		     INT             *iedgnod       
                    ) 
{
INT    irn,irow;	      
INT    nd;                 
DOUBLE facl;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_calgalif_vhf");
#endif

fdyn = alldyn[genprob.numff].fdyn;

facl = fac * fdyn->thsl;	      
nd = NUMDOF_FLUID2*iel;



/*----------------------------------------------------------------------*
   Calculate iteration force vector:
               /
     THETA*dt |  psi * Ux * phi,x     d_gamma_FS
             /  
 *----------------------------------------------------------------------*/ 
for(irn=0;irn<ngnode;irn++)
{
   irow = nd+iedgnod[irn];
   eforce[irow]	 += funct[irn]*velint[0]*phiderx*facl;
} /* end loop over irow */
	      
/*----------------------------------------------------------------------*/ 
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_calgalif_vhf */

#endif
/*! @} (documentation module close)*/
 
