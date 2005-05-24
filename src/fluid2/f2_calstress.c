/*!----------------------------------------------------------------------
\file
\brief calculation of fluid stresses

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

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
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*!---------------------------------------------------------------------
\brief calculation of fluid stresses

<pre>                                                         genk 10/02

calculation of the stress tensor sigma, which es stored at the nodes

    sigma = -p_real*I + 2*nue * eps(u)

</pre>

\param   viscstr    INT         (i)    include viscose stresses yes/no
\param  *ele        ELMENT      (i/o)  actual element
\param **evel       DOUBLE      (-)    element velocities
\param  *epre       DOUBLE      (-)    element pressure
\param **funct      DOUBLE      (-)    shape functions
\param **deriv      DOUBLE      (-)    natural deriv. of shape funct.
\param **derxy      DOUBLE      (-)    global deriv. of sape funct.
\param **vderxy     DOUBLE      (-)    global vel. deriv
\param **xjm        DOUBLE      (-)    jacobian matrix
\param **xyze       DOUBLE      (-)    element coordinates
\param **sigmaint   DOUBLE      (-)    stresses at GAUSS point
\param  *ipos                   (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f2_calelestress(
                      INT             viscstr,
       	              ELEMENT        *ele,
		      DOUBLE        **evel,
		      DOUBLE         *epre,
		      DOUBLE         *funct,
		      DOUBLE        **deriv,
		      DOUBLE        **derxy,
		      DOUBLE        **vderxy,
		      DOUBLE        **xjm,
		      DOUBLE        **xyze,
		      DOUBLE        **sigmaint,
                      ARRAY_POSITION *ipos,
		      INT             is_relax
		    )
{
INT     i,j,node,lr,ls;      /* some counters                          */
INT     iel,nir,nis;         /* number of nodes/integr. points         */
INT     intc,icode;          /* flags                                  */
INT     actmat;              /* actual material number                 */
INT     iv;                  /* counter for GAUSS points               */
INT     pos;                 /* position, where to read from           */
DOUBLE  preint,det,val;      /* element values                         */
DOUBLE  e1,e2,r,s;
DIS_TYP typ;	             /* element displacement type              */
DOUBLE  dens,visc,twovisc;   /* material parameters                    */
DOUBLE  fpar[MAXGAUSS];      /* working array                          */
NODE   *actnode;             /* actual node                            */
FLUID_DYNAMIC  *fdyn;
FLUID_DATA     *data;
#ifdef D_FSI
NODE   *actfnode, *actanode;
GNODE  *actfgnode;
#endif

#ifdef DEBUG
dstrc_enter("f2_calelestress");
#endif

/*----------------------------------------------------- initialisation */
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity*dens; /* here we need dynamic viscosity! */
typ  = ele->distyp;

if (is_relax) pos = ipos->relax;
else          pos = ipos->velnp;

switch(viscstr)
{
case 0: /* only real pressure */
   for (i=0;i<iel;i++)
   {
      actnode=ele->node[i];
      ele->e.f2->stress_ND.a.da[i][0]=-actnode->sol_increment.a.da[pos][2]*dens;
      ele->e.f2->stress_ND.a.da[i][1]=-actnode->sol_increment.a.da[pos][2]*dens;
      ele->e.f2->stress_ND.a.da[i][2]= ZERO;
   }
break;
case 1: /* real pressure + viscose stresses */
/* sigma = -p_real*I + 2*nue * eps(u) */
/* calculate deformation velocity tensor */
/*------------------------------------------------ get integraton data  */
   switch (typ)
   {
   case quad4: case quad8: case quad9:  /* --> quad - element */
      icode = 2;
      nir = ele->e.f2->nGP[0];
      nis = ele->e.f2->nGP[1];
   break;
   case tri3: case tri6: /* --> tri - element */
      icode = 2;
      nir  = ele->e.f2->nGP[0];
      nis  = 1;
      intc = ele->e.f2->nGP[1];
   break;
   default:
      dserror("typ unknown!");
   } /* end switch(typ) */
/*-------- set element velocities, real pressure and coordinates -------*/
   for (j=0;j<iel;j++)
   {
      evel[0][j] = ele->node[j]->sol_increment.a.da[pos][0];
      evel[1][j] = ele->node[j]->sol_increment.a.da[pos][1];
      epre[j]    = ele->node[j]->sol_increment.a.da[pos][2]*dens;
   }/*end for (j=0;j<iel;j++) */
   switch (ele->e.f2->is_ale)
   {
   case 0:
      for (j=0;j<iel;j++)
      {
         xyze[0][j] = ele->node[j]->x[0];
         xyze[1][j] = ele->node[j]->x[1];
      }
   break;
#ifdef D_FSI
   case 1:
      for (j=0;j<iel;j++)
      {
         actfnode = ele->node[j];
         actfgnode = actfnode->gnode;
         actanode = actfgnode->mfcpnode[genprob.numaf];
         xyze[0][j] = ele->node[j]->x[0]+ actanode->sol_mf.a.da[1][0];
         xyze[1][j] = ele->node[j]->x[1]+ actanode->sol_mf.a.da[1][1];
      }
   break;
#endif
   default:
      dserror("elment flag is_ale out of range!");
   }
/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
   iv=0;
   twovisc=TWO*visc;
   for (lr=0;lr<nir;lr++)
   {
      for (ls=0;ls<nis;ls++)
      {
/*--------------- get values of  shape functions and their derivatives */
         switch(typ)
         {
         case quad4: case quad8: case quad9:   /* --> quad - element */
	    e1   = data->qxg[lr][nir-1];
	    e2   = data->qxg[ls][nis-1];
            f2_rec(funct,deriv,NULL,e1,e2,typ,icode);
         break;
         case tri3: case tri6:   /* --> tri - element */
	    e1   = data->txgr[lr][intc];
	    e2   = data->txgs[lr][intc];
	    f2_tri(funct,deriv,NULL,e1,e2,typ,icode);
         break;
         default:
            dserror("typ unknown!");
         } /* end switch(typ) */
/*-------------------------------------------- compute Jacobian matrix */
         f2_jaco2(xyze,deriv,xjm,&det,iel);
/*------------------------------------------- compute global derivates */
         f2_gder(derxy,deriv,xjm,det,iel);
/*--------------------- get velocity  derivatives at integration point */
         f2_vder(vderxy,derxy,evel,iel);
/*---------------------------------- get pressure at integration point */
	 preint=f2_scali(funct,epre,iel);
/*---------------
                         | Ux,x    Ux,y |
                vderxy = |              |
                         | Uy,x    Uy,y |

                             | Ux,x+Ux,x   Ux,y+Uy,x |
                eps(u) = 1/2 |		             |
                             |   symm.     Uy,y+Uy,y |

                SIGMA = -p_real*I + 2*nue * eps(u)                     */
         sigmaint[0][iv] = -preint + twovisc*vderxy[0][0];
	 sigmaint[1][iv] = -preint + twovisc*vderxy[1][1];
	 sigmaint[2][iv] = (vderxy[0][1] + vderxy[1][0])*visc;
	 iv++;
      } /* end loop over nis */
   } /* end loop over nir */
   /*------------------------------- extrapolate stresses to the nodes */
   for (i=0;i<3;i++)
   {
      switch (typ)
      {
      case quad4: case quad8: case quad9:
         f2_recex(NULL,fpar,ZERO,ZERO,&(sigmaint[i][0]),iv,0);
      break;
      case tri3: case tri6:
         f2_triex(NULL,fpar,ZERO,ZERO,&(sigmaint[i][0]),iv,0);
      break;
      default:
         dserror("typ unknown!\n");
      }
      for (node=0;node<iel;node++)
      {
         actnode = ele->node[node];
         r = f2_rsn(node,0,iel);
         s = f2_rsn(node,1,iel);
         switch (typ)
         {
         case quad4: case quad8: case quad9:
            f2_recex(&val,fpar,r,s,NULL,iv,2);
	 break;
         case tri3: case tri6:
            f2_triex(&val,fpar,r,s,NULL,iv,2);
	 break;
         default:
            dserror("typ unknown!\n");
	 }
         ele->e.f2->stress_ND.a.da[node][i]=val;
      }
   }
break;
default:
   dserror("parameter viscstr out of range!\n");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calelestress */

#endif
/*! @} (documentation module close)*/
