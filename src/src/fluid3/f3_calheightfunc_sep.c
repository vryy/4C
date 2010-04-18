/*!----------------------------------------------------------------------
\file
\brief element integration for heightfunction

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup FLUID3
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
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

<pre>                                                         genk 04/04

element integration of height function evaluation

</pre>
\param   *ele        ELEMENT           actual element
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
void f3_calint_hfsep(
                     ELEMENT           *ele,
                     DOUBLE            *funct,
                     DOUBLE           **deriv,
                     DOUBLE           **xjm,
                     DOUBLE           **wa2,
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
INT    icode=2;
INT    lr,ls;
INT    intc=0;
DOUBLE facr=0.0,facs=0.0;
DOUBLE phiderxy[2],e1,e2,fac,det;
DOUBLE phiintn, phiintng;
FLUID_DYNAMIC *fdyn;
FLUID_DATA    *data;

#ifdef DEBUG
dstrc_enter("f3_calint_hfsep");
#endif

fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

for (lr=0;lr<nil;lr++)
{
for (ls=0;ls<nil;ls++)
{
   /*--------------- get values of  shape functions and their derivatives */
   switch(typ)
   {
   case quad4: case quad8: case quad9:   /* hex element -->quad - surface */
      e1   = data->qxg[lr][nil-1];
      facr = data->qwgt[lr][nil-1];
      e2   = data->qxg[ls][nil-1];
      facs = data->qwgt[ls][nil-1];
      f3_rec(funct,deriv,NULL,e1,e2,typ,icode);
   break;
   case tri3: case tri6:   /* tet element --> tri - surface */
      dserror("hf tri free surface not implemented yet!\n");
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      f3_tri(funct,deriv,NULL,e1,e2,typ,icode);
   break;
   default:
      dserror("typ unknown!");
   } /* end switch(typ) */
   /*------------------------------- compute jacobian determinant */
   f3_edgejaco(xyze,deriv,xjm,&det,iedgnod,ngnode,ele);
   fac = det*facr*facs;
   /*--------------------------------------- compute global derivatives */
   f3_edgegder(derxy,deriv,xjm,det,ngnode);
   /*------------------------- get velocity at integration point at n+1 */
   f3_edgeveli(velint,funct,evelng,ngnode,iedgnod);
   /*--------------------------- get velocity at integration point at n */
   f3_edgeveli(vel2int,funct,evelng,ngnode,iedgnod);
   /*-------------------- get global derivative of height function at n */
   f3_phider(phiderxy,derxy,ephin,ngnode,iedgnod);
   /*-------------------- get height function at integration point at n */
   phiintn = f3_phii(funct,ephin,iedgnod,ngnode);
   /*------------------ get height function at integration point at n+1 */
   phiintng= f3_phii(funct,ephing,iedgnod,ngnode);
   /*-------------------------- get stab parameter at integration point */
   if (fdyn->hf_stab>0)
   f3_stabpar_hfsep(ele,wa2,xjm,xyze,velint,phiintn,
                    phiintng,phiderxy,iedgnod,ngnode,typ);
   /*-------------------------------- compute Galerkin part of matrices */
   f3_calmat_vhf_sep(ele,estif,ngnode,funct,derxy,velint,fac);
   /*----------------------------------------------- compute RHS at n+1 */
   f3_caliterhs_vhf_sep(ele,eiforce,ngnode,funct,derxy,velint,fac);
   /*------------------------------------------------- compute RHS at n */
   f3_caltimerhs_vhf_sep(ele,eiforce,ngnode,funct,derxy,velint,
                         vel2int,phiintn,phiderxy,fac);
} /* end loop over ls */
} /* end loop over lr */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_calint_hfsep */

/*!---------------------------------------------------------------------
\brief control routine for element integration of fluid2

<pre>                                                         genk 04/04

stabilisation parameter for height function

</pre>
\param      *ele        ELEMENT           actual element
\param     **xyze       DOUBLE            element co-ord.
\param     **evlng      DOUBLE            element vels
\param      *funct      DOUBLE            shape funcs
\param      *velint     DOUBLE            vel at intpoint
\param       phiintn    DOUBLE            phi at time (n)
\param       phintng    DOUBLE            phi at time (ng)
\param       phiderxy   DOUBLE            deriv. of phi w.r.t xy
\param      *iedgnod    INT               edge nodes
\param       ngnode     INT               num of nodes at actual edge
\param       typ        DIS_TYP           discr. typ
\return void

------------------------------------------------------------------------*/
void f3_stabpar_hfsep(
                        ELEMENT          *ele,
                        DOUBLE          **deriv,
                        DOUBLE          **xjm,
                        DOUBLE          **xyze,
                        DOUBLE           *velint,
                        DOUBLE            phiintn,
                        DOUBLE            phiintng,
                        DOUBLE           *phiderxy,
                        INT              *iedgnod,
                        INT               ngnode,
                        DIS_TYP           typ
                     )
{
DOUBLE h,velno;
DOUBLE phidot,res;
DOUBLE area;
DOUBLE funct[MAXNOD_F3];
const DOUBLE C = ONE;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("f3_stabpar_hfsep");
#endif

dsassert(ngnode==4,"stabilisation of hf only for two noded edges!\n");
fdyn   = alldyn[genprob.numff].fdyn;
/*---------------------------------------------------- get element size */
if (fdyn->hf_stab==1) /* stabparameter at element centre */
{
   dserror("hf stabilisation at element centre not possible!\n");
}

/*------------------------------- compute area using 1 integrationpoint */
switch(typ)
{
case quad4:    /* --> quad - element */
   f3_rec(funct,deriv,NULL,ZERO,ZERO,typ,2);
break;
default:
      dserror("typ unknown!");
} /* end switch(typ) */

f3_edgejaco(xyze,deriv,xjm,&area,iedgnod,ngnode,ele);

/*---------------------------------------------- determine element size */
h = sqrt(area);

/*---------------------------------------------------- norm of velocity */
velno = DSQR(velint[0])+DSQR(velint[1]);
velno = sqrt(velno);
velno = DMAX(velno,EPS6);

phidot = (phiintng-phiintn)/fdyn->dta;


/*--------------------------- compute tau_SUPG according to SOULAIMANI:
                h
tau_SUPG = -------------,    where
            2 sqrt(|Us|)

       h is the measure of the surface element size                     */
fdyn->tau[3]=h/(TWO*velno);

/*------------------------------------- compute tau_DC according to BEHR:


tau_DC = C * h * | phidot + Ux*phi,x + Uy*phi,y - Uz |

*/
res = phidot + velint[0]*phiderxy[0] + velint[1]*phiderxy[1] - velint[2];
res = FABS(res);
fdyn->tau[4] = C*h*res;

#if 0
printf("ELEMENT   %d    tau_DC = %lf \n",ele->Id,fdyn->tau[4]);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_stabpar_hfsep */


/*!---------------------------------------------------------------------
\brief control routine for element integration of fluid2

<pre>                                                         genk 04/04

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
void f3_calmat_vhf_sep(
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
dstrc_enter("f3_calmat_vhf_sep");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

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

/*----------------------------------------------------------------------*
   Calculate full Galerkin part of matrix Kpsipsi:
               /
   + THETA*dt |  psi * Uy * phi,y    d_gamma_FS
             /
 *----------------------------------------------------------------------*/
c = fac*fdyn->thsl;
for(icn=0;icn<ngnode;icn++)
{
   for(irn=0;irn<ngnode;irn++)
   {
      estif[irn][icn] += funct[irn]*velint[1]*derxy[1][icn]*c;
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
   Calculate full stabilisation part of matrix Mpsipsi:
     /
 +  |   tau_SUPG * Uy * psi,y * phi  d_gamma_FS
   /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[1];
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]	  += derxy[1][irn]*funct[icn]*c;
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
 + THETA*dt |   tau_SUPG * Uy * psi,y * Ux * phi,x  d_gamma_FS
           /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsl*tau_supg*velint[1]*velint[0];
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]     += derxy[1][irn]*derxy[0][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
             /
 + THETA*dt |   tau_SUPG * Ux * psi,x * Uy * phi,y  d_gamma_FS
           /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsl*tau_supg*velint[0]*velint[1];
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]     += derxy[0][irn]*derxy[1][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
             /
 + THETA*dt |   tau_SUPG * Uy * psi,y * Uy * phi,y  d_gamma_FS
           /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsl*tau_supg*velint[1]*velint[1];
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]     += derxy[1][irn]*derxy[1][icn]*c;
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

/*----------------------------------------------------------------------*
   Calculate full stabilisation part of matrix Kpsipsi:
             /
 + THETA*dt |   tau_DC * psi,y * phi,y  d_gamma_FS
           /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsl*tau_dc;
   for(icn=0;icn<ngnode;icn++)
   {
      for(irn=0;irn<ngnode;irn++)
      {
         estif[irn][icn]     += derxy[1][irn]*derxy[1][icn]*c;
      } /* end loop over irn */
   } /* end loop over icn */
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_calmat_vhf_sep */

/*!---------------------------------------------------------------------
\brief Iteration RHS for vertical heightfunction seperate

<pre>                                                         genk 04/04

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
void f3_caliterhs_vhf_sep(
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
DOUBLE  c;
DOUBLE tau_supg;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("f3_caliterhs_vhf_sep");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

/*----------------------------------------------------------------------*
   Calculate iteration force vector:
                 /
   + THETA*dt   |  psi * Uz     d_gamma_FS
               /
*----------------------------------------------------------------------*/
c = fac*fdyn->thsl;
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    += funct[irn]*velint[2]*c;
} /* end loop over irow */

/*------------------------------------------------------- STABILISATION */
if (fdyn->hf_stab>0)
{
   tau_supg = fdyn->tau[3];
/*----------------------------------------------------------------------*
   Calculate stabilisation part of iteration force vector:
                  /
   + THETA*dt    | tau_SUPG * Ux * psi,x * Uz     d_gamma_FS
                /
  *----------------------------------------------------------------------*/
   c = fac*fdyn->thsl*tau_supg*velint[0]*velint[2];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    += derxy[0][irn]*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation part of iteration force vector:
                  /
   + THETA*dt    | tau_SUPG * Uy * psi,y * Uz     d_gamma_FS
                /
  *----------------------------------------------------------------------*/
   c = fac*fdyn->thsl*tau_supg*velint[1]*velint[2];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    += derxy[1][irn]*c;
   }

}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_caliterhs_vhf_sep */

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
\param   *phiderxy   DOUBLE            deriv. of phi w.r.t. x
\param    fac        DOUBLE            integr. factor
\return void
\return void

------------------------------------------------------------------------*/
void f3_caltimerhs_vhf_sep(
                           ELEMENT                *ele,
                           DOUBLE                 *eforce,
                           INT                     ngnode,
                           DOUBLE                 *funct,
                           DOUBLE                **derxy,
                           DOUBLE                 *velint, /* at n+1 */
                           DOUBLE                 *vel2int,/* at n */
                           DOUBLE                  phiint,
                           DOUBLE                 *phiderxy,
                           DOUBLE                  fac
                          )
{
INT irn;
DOUBLE  c;
DOUBLE tau_supg,tau_dc;
FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("f3_caltimerhs_vhf_sep");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

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
   + (1-THETA)*dt  |  psi * (Uz)_old      d_gamma_FS
                  /
 *----------------------------------------------------------------------*/
c = fac*fdyn->thsr;
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    += funct[irn]*vel2int[2]*c;
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
   eforce[irn]    -= funct[irn]*vel2int[0]*phiderxy[0]*c;
} /* end loop over irow */

/*----------------------------------------------------------------------*
   Calculate time force vector:
                    /
   - (1-THETA)*dt  |  psi * (Uy)_old  phi,y_old    d_gamma_FS
                  /
 *----------------------------------------------------------------------*/
c =  fac*fdyn->thsr;
for(irn=0;irn<ngnode;irn++)
{
   eforce[irn]    -= funct[irn]*vel2int[1]*phiderxy[1]*c;
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
 +  |   tau_SUPG * Uy * psi,y * phi_old  d_gamma_FS
   /
 *----------------------------------------------------------------------*/
   c=fac*tau_supg*velint[1];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]     += derxy[1][irn]*phiint*c;
   } /* end loop over irn */


/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 + (1-THETA)*dt   |   tau_SUPG * Ux * psi,x * (Uz)_old  d_gamma_FS
                 /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]     += derxy[0][irn]*vel2int[2]*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 + (1-THETA)*dt   |   tau_SUPG * Uy * psi,y * (Uz)_old  d_gamma_FS
                 /
 *----------------------------------------------------------------------*/
   c=fac*fdyn->thsr*tau_supg*velint[1];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]     += derxy[1][irn]*vel2int[2]*c;
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
      eforce[irn]    -= derxy[0][irn]*vel2int[0]*phiderxy[0]*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 - (1-THETA)*dt   |   tau_SUPG * Uy * psi,y * (Ux)_old * phi,x_old  d_gamma_FS
                 /
 *----------------------------------------------------------------------*/
   c= fac*fdyn->thsr*tau_supg*velint[1];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    -= derxy[1][irn]*vel2int[0]*phiderxy[0]*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 - (1-THETA)*dt   |   tau_SUPG * Ux * psi,x * (Uy)_old * phi,y_old  d_gamma_FS
                 /
 *----------------------------------------------------------------------*/
   c= fac*fdyn->thsr*tau_supg*velint[0];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    -= derxy[0][irn]*vel2int[1]*phiderxy[1]*c;
   }

/*----------------------------------------------------------------------*
   Calculate stabilisation of time force vector:
                   /
 - (1-THETA)*dt   |   tau_SUPG * Uy * psi,y * (Uy)_old * phi,y_old  d_gamma_FS
                 /
 *----------------------------------------------------------------------*/
   c= fac*fdyn->thsr*tau_supg*velint[1];
   for(irn=0;irn<ngnode;irn++)
   {
      eforce[irn]    -= derxy[1][irn]*vel2int[1]*phiderxy[1]*c;
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
      eforce[irn]    -= derxy[0][irn]*phiderxy[0]*c;
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
      eforce[irn]    -= derxy[1][irn]*phiderxy[1]*c;
   }
}



/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f3_caltimerhs_vhf_sep */

#endif
#endif
