/*!----------------------------------------------------------------------
\file
\brief element integration for heightfunction

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
#include "../fluid_full/fluid_prototypes.h"
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
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*!--------------------------------------------------------------------- 
\brief control routine for element integration of fluid2

<pre>                                                         genk 06/03

stabilisation parameter for height function
			     
</pre>
\param      *ele        ELEMENT           actual element
\param     **xyze       DOUBLE            element co-ord.
\param     **evlng      DOUBLE            element vels
\param      *funct      DOUBLE            shape funcs
\param      *velint     DOUBLE            vel at intpoint
\param       phiintn    DOUBLE            phi at time (n)
\param       phintng    DOUBLE            phi at time (ng)
\param       phiderx    DOUBLE            deriv. of phi w.r.t x
\param       e1         DOUBLE            integr. weight
\param      *iedgnod    INT               edge nodes
\param       ngnode     INT               num of nodes at actual edge
\param       typ        DIS_TYP           discr. typ
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_stabpar_hfsep(
                        ELEMENT          *ele,
                        DOUBLE          **xyze,
                        DOUBLE          **evelng,
                        DOUBLE           *funct,
                        DOUBLE           *velint,
                        DOUBLE            phiintn,
                        DOUBLE            phiintng,
                        DOUBLE            phiderx,
                        DOUBLE            e1,
                        INT              *iedgnod,
                        INT               ngnode,
                        DIS_TYP           typ 
                     )
{
DOUBLE x0,x1,h,velno;
DOUBLE phidot,res;
const DOUBLE C = ONE;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_stabpar_hfsep");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

/*--------------------------- compute tau_SUPG according to SOULAIMANI:
                h
tau_SUPG = -------------,    where
            2 sqrt(|Ux|)
     
       h is the measure of the surface element size                     */

dsassert(ngnode==2,"stabilisation of hf only for two noded edges!\n");        
/*---------------------------------------------------- get element size */
x0 = xyze[0][iedgnod[0]];
x1 = xyze[0][iedgnod[1]];
h  = FABS(x0-x1);

if (fdyn->hf_stab==1) /* stabparameter at element centre */
{
   dserror("hf stabilisation at element centre not possible!\n");
   /*---------------------------- get shape functions at element centre */
   f2_degrectri(funct,NULL,e1,typ,0);
   /*----------------------------------- get velocity at element centre */
   f2_edgeveci(velint,funct,evelng,ngnode,iedgnod);
}

/*---------------------------------------------------- norm of velocity */
velno = FABS(velint[0]); 
velno = sqrt(velno);
velno = DMAX(velno,EPS6);

phidot = (phiintng-phiintn)/fdyn->dta;


/*--------------------------- compute tau_SUPG according to SOULAIMANI:
                h
tau_SUPG = -------------,    where
            2 sqrt(|Ux|)
     
       h is the measure of the surface element size                     */
fdyn->tau[3]=h/(TWO*velno);

/*------------------------------------- compute tau_DC according to BEHR:

                  
tau_DC = C * h * | phidot + Ux*phi,x - Uy |
                 
*/
res = phidot + velint[0]*phiderx - velint[1];
res = FABS(res);
fdyn->tau[4] = C*h*res; 


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_stabpar_hfsep */


/*!--------------------------------------------------------------------- 
\brief control routine for element integration of fluid2

<pre>                                                         genk 06/03

element integration of height function evaluation
			     
</pre>
\param   *ele        ELEMENT           actual element
\param   *data       FLUID_DATA        integr. data
\param   *funct      DOUBLE            shape funcs
\param  **deriv      DOUBLE            nat. deriv. of shape funcs
\param  **xjm        DOUBLE            jacobian matrix
\param  **xyze       DOUBLE            element co-ord.
\param    ngnode     INT               num nodes at actual edge
\param    nil        INT               num of integr. points along edge
\param   *iedgnod    INT               edge nodes
\param   *velint     DOUBLE            vel at int. point
\param   *vel2int    DOUBLE            another vel at int. point
\param  **evelng     DOUBLE            nodal vels at n+g
\param  **eveln      DOUBLE            nodal vels at n
\param   *ephing     DOUBLE            height func at n+g
\param   *ephin      DOUBLE            height func at n
\param  **derxy      DOUBLE            global derivs w.r.t. x
\param    typ        DIS_TYP           discr. typ
\param  **estif      DOUBLE            element stiffness matrix
\param   *eiforce    DOUBLE            element RHS 
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calint_hfsep(
                     ELEMENT           *ele,
                     FLUID_DATA        *data,
                     DOUBLE            *funct,
                     DOUBLE           **deriv,
                     DOUBLE           **xjm,
                     DOUBLE           **xyze,
                     INT                ngnode,
                     INT                nil,
                     INT               *iedgnod,
                     DOUBLE            *velint,
                     DOUBLE            *vel2int,
                     DOUBLE           **evelng,
                     DOUBLE           **eveln,
                     DOUBLE            *ephing,
                     DOUBLE            *ephin,
                     DOUBLE           **derxy,
                     DIS_TYP            typ,
                     DOUBLE           **estif,
                     DOUBLE            *eiforce
                     )
{
INT k;
INT    lr,node;
DOUBLE phiderx,e1,facr,fac,det;
DOUBLE phiintn, phiintng;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_calint_hfsep");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

for (lr=0;lr<nil;lr++)
{
   /*------------- get values of  shape functions and their derivatives */
   e1	= data->qxg[lr][nil-1];
   facr = data->qwgt[lr][nil-1];
   f2_degrectri(funct,deriv,e1,typ,1);
   
   /*------------------------------------- compute jacobian determinant */
   /*------------------------ integration is performed along the x-axis */
   xjm[0][0] = ZERO ;
   for (k=0; k<ngnode; k++) /* loop all nodes of the element */
   {
      node=iedgnod[k];
      xjm[0][0] += deriv[0][k] * xyze[0][node] ;
   } /* end loop over iel */
   det = FABS(xjm[0][0]);
   fac = det*facr;

   /*--------------------------------------- compute global derivatives */ 
   f2_edgegder(deriv,derxy,xjm,ngnode);   
   /*------------------------- get velocity at integration point at n+1 */
   f2_edgeveci(velint,funct,evelng,ngnode,iedgnod);
   /*--------------------------- get velocity at integration point at n */
   f2_edgeveci(vel2int,funct,eveln,ngnode,iedgnod);
   /*-------------------- get global derivative of height function at n */      
   phiderx = f2_phider(derxy,ephin,ngnode,iedgnod);
   /*-------------------- get height function at integration point at n */
   phiintn = f2_edgescali(funct,ephin,iedgnod,ngnode);
   /*------------------ get height function at integration point at n+1 */
   phiintng = f2_edgescali(funct,ephing,iedgnod,ngnode);

   /*-------------------------- get stab parameter at integration point */
   if (fdyn->hf_stab==2)
      f2_stabpar_hfsep(ele,xyze,NULL,NULL,velint,phiintn,phiintng,
                       phiderx,e1,iedgnod,ngnode,typ);

   /*---------------- compute Galerkin & stabilisation part of matrices */  
   f2_calmat_vhf_sep(ele,estif,ngnode,funct,derxy,velint,fac);
   /*----------------------------------------------- compute RHS at n+1 */
   f2_caliterhs_vhf_sep(ele,eiforce,ngnode,funct,derxy,velint,fac);
   /*------------------------------------------------- compute RHS at n */
   f2_caltimerhs_vhf_sep(ele,eiforce,ngnode,funct,derxy,velint,
                         vel2int,phiintn,phiderx,fac);
} /* end if loop over integration points */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_calint_hfsep */


/*!--------------------------------------------------------------------- 
\brief control routine for element integration of fluid2

<pre>                                                         genk 06/03

evaluate stiffness matrix of vertical heightfunction
			     
</pre>
\param   *ele        ELEMENT        actuale element
\param  **estif      DOUBLE         element stiffness matrix
\param    ngnode     INT            num nodes at actual edge
\param   *funct      DOUBLE         shape funcs
\param  **derxy      DOUBLE         global derivs. w.r.t x
\param   *velint     DOUBLE         vel at integr. point
\param    fac        DOUBLE         integration factor
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calmat_vhf_sep(
		       ELEMENT                *ele,
		       DOUBLE                **estif,
		       INT                     ngnode,
		       DOUBLE                 *funct,
		       DOUBLE                **derxy,
		       DOUBLE                 *velint,
		       DOUBLE                  fac
	              )
{
INT irn,icn;
DOUBLE c;
DOUBLE tau_supg, tau_dc;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_calmat_vhf_sep");
#endif

fdyn = alldyn[genprob.numff].fdyn;

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Mpsipsi:
    /
   |  psi * phi    d_gamma_FS
  /
 *----------------------------------------------------------------------*/
for(icn=0;icn<ngnode;icn++)
{
   for(irn=0;irn<ngnode;irn++)
   {
      estif[irn][icn]	  += funct[irn]*funct[icn]*fac;
   } /* end loop over irn */
} /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kpsipsi:
               /
   + THETA*dt |  psi * Ux * phi,x    d_gamma_FS
             /
 *----------------------------------------------------------------------*/
c = fac*fdyn->thsl;
for(icn=0;icn<ngnode;icn++)
{
   for(irn=0;irn<ngnode;irn++)
   {
      estif[irn][icn] += funct[irn]*velint[0]*derxy[0][icn]*c;
   } /* end loop over irn */
} /* end loop over icn */

/*------------------------------------------------------- STABILISATION */
if (fdyn->hf_stab>0)
{
   tau_supg = fdyn->tau[3];
   tau_dc   = fdyn->tau[4];
/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Mpsipsi:
     /
 +  |   tau_SUPG * Ux * psi,x * phi  d_gamma_FS   
   /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]	  += derxy[0][irn]*funct[icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */


/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
             /
 + THETA*dt |   tau_SUPG * Ux * psi,x * Ux * phi,x  d_gamma_FS   
           /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsl*tau_supg*velint[0]*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]     += derxy[0][irn]*derxy[0][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
             /
 + THETA*dt |   tau_DC * psi,x * phi,x  d_gamma_FS   
           /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsl*tau_dc;
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]     += derxy[0][irn]*derxy[0][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_calmat_vhf_sep */


/*!--------------------------------------------------------------------- 
\brief Iteration RHS for vertical heightfunction seperate

<pre>                                                         genk 06/03

evaluate RHS for height function
			     
</pre>
\param   *ele        ELEMENT           the actual element
\param   *eforce     DOUBLE            element RHS
\param    ngnode     INT               num nodes on actual edge
\param   *funct      DOUBLE            shape funcs
\param  **derxy      DOUBLE            global deriv w.r.t to x
\param   *velint     DOUBLE            vel at integr. point
\param    fac        DOUBLE            integr. factor
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_caliterhs_vhf_sep(
 	 	           ELEMENT                *ele,
		           DOUBLE                 *eforce,
		           INT                     ngnode,
		           DOUBLE                 *funct,
                           DOUBLE                **derxy,
		           DOUBLE                 *velint,
			   DOUBLE                  fac
	                  )
{
INT irn;
DOUBLE c;
DOUBLE tau_supg;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_caliterhs_vhf_sep");
#endif

fdyn    = alldyn[genprob.numff].fdyn;

	      
/*----------------------------------------------------------------------*
   Calculate iteration force vector:
                 /
   + THETA*dt   |  psi * Uy     d_gamma_FS
               /  
*----------------------------------------------------------------------*/ 
c = fac*fdyn->thsl;
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    += funct[irn]*velint[1]*c;
} /* end loop over irow */

/*------------------------------------------------------- STABILISATION */
if (fdyn->hf_stab>0)
{
   tau_supg = fdyn->tau[3];   
/*----------------------------------------------------------------------*
   Calculate stabilisation part of iteration force vector:
                  /
   + THETA*dt    | tau_SUPG * Ux * psi,x * Uy     d_gamma_FS
                /  
  *----------------------------------------------------------------------*/ 
   c = fac*fdyn->thsl*tau_supg*velint[0]*velint[1];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    += derxy[0][irn]*c;
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_caliterhs_vhf_sep */

/*!--------------------------------------------------------------------- 
\brief Time RHS for vertical heightfunction seperate

<pre>                                                         genk 06/03

evaluate RHS for heightfunction
			     
</pre>
\param   *ele        ELEMENT           the actual element
\param   *eforce     DOUBLE            element RHS
\param    ngnode     INT               num nodes on actual edge
\param   *funct      DOUBLE            shape funcs
\param  **derxy      DOUBLE            global deriv w.r.t to x
\param   *velint     DOUBLE            vel at integr. point
\param   *vel2int    DOUBLE            vel at integr. point
\param    phiint     DOUBLE            phi at integr. point
\param    phiderx    DOUBLE            deriv. of phi w.r.t. x
\param    fac        DOUBLE            integr. factor
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_caltimerhs_vhf_sep(
                           ELEMENT                *ele,
                           DOUBLE                 *eforce,
                           INT                     ngnode,
                           DOUBLE                 *funct,
                           DOUBLE                **derxy,
                           DOUBLE                 *velint, /* at n+1 */
                           DOUBLE                 *vel2int,/* at n */
                           DOUBLE                  phiint,
                           DOUBLE                  phiderx,
                           DOUBLE                  fac
                          )
{
INT irn;
DOUBLE c;
DOUBLE tau_supg,tau_dc;
FLUID_DYNAMIC *fdyn;


#ifdef DEBUG 
dstrc_enter("f2_caltimerhs_vhf_sep");
#endif

fdyn    = alldyn[genprob.numff].fdyn;


/*----------------------------------------------------------------------*
   Calculate intertia forces of time force vector:
      /
   + |  psi * phi     d_gamma_FS
    /  
 *----------------------------------------------------------------------*/ 
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    += funct[irn]*phiint*fac;
} /* end loop over irow */	 

	      
/*----------------------------------------------------------------------*
   Calculate time force vector:
                    /
   + (1-THETA)*dt  |  psi * (Uy)_old      d_gamma_FS
                  /  
 *----------------------------------------------------------------------*/ 
c = fac*fdyn->thsr;
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    += funct[irn]*vel2int[1]*c; 
} /* end loop over irow */

/*----------------------------------------------------------------------*
   Calculate time force vector:
                    /
   - (1-THETA)*dt  |  psi * (Ux)_old  phi,x_old    d_gamma_FS
                  /  
 *----------------------------------------------------------------------*/ 
c =  fac*fdyn->thsr;
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    -= funct[irn]*vel2int[0]*phiderx*c;
} /* end loop over irow */


/*------------------------------------------------------- STABILISATION */
if (fdyn->hf_stab>0)
{
   tau_supg = fdyn->tau[3];
   tau_dc   = fdyn->tau[4];   
/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
     /
 +  |   tau_SUPG * Ux * psi,x * phi_old  d_gamma_FS   
   /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]     += derxy[0][irn]*phiint*c;
   } /* end loop over irn */

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 + (1-THETA)*dt   |   tau_SUPG * Ux * psi,x * (Uy)_old  d_gamma_FS   
                 /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]     += derxy[0][irn]*vel2int[1]*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 - (1-THETA)*dt   |   tau_SUPG * Ux * psi,x * (Ux)_old * phi,x_old  d_gamma_FS   
                 /
 *----------------------------------------------------------------------*/
   c= fac*fdyn->thsr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    -= derxy[0][irn]*vel2int[0]*phiderx*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 - (1-THETA)*dt   |   tau_DC * psi,x * phi,x_old  d_gamma_FS   
                 /
 *----------------------------------------------------------------------*/
   c= fac*fdyn->thsr*tau_dc;
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    -= derxy[0][irn]*phiderx*c;
   }
}



/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_caltimerhs_vhf_sep */

#endif
