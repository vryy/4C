/*!----------------------------------------------------------------------
\file
\brief compute element contributions to tangent and internal force vector

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../ls/ls_prototypes.h"
#include "xfem_prototypes.h"
/*!
\addtogroup XFEM
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB        genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL      *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA               *alldyn;
extern struct _XFEM_DATA      xfem_data;


/*!---------------------------------------------------------------------
\brief compute element contributions to tangent and internal force vector

<pre>                                                            irhan 05/04

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2_xfem element is calculated

</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *hasext    INT             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    Jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadn    DOUBLE	   (-)    ele dead load (selfweight) at n
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *vel2int   DOUBLE	   (-)    vel at integration point
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void

------------------------------------------------------------------------*/
void xfem_f2_calint(
  FLUID_DATA      *data,
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE         **estif,
  DOUBLE         **emass,
  DOUBLE          *etforce,
  DOUBLE          *eiforce,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE          *epren,
  DOUBLE          *edeadn,
  DOUBLE          *edeadng,
  DOUBLE          *velint,
  DOUBLE          *vel2int,
  DOUBLE          *covint,
  DOUBLE         **vderxy,
  DOUBLE          *pderxy,
  DOUBLE         **vderxy2,
  DOUBLE         **wa1,
  DOUBLE         **wa2
  )
{
  INT       i;
  INT       iel;
  INT       nir,nis;
  INT       actmat;     /* material number of the element */
  INT       ihoel=0;    /* flag for higher order elements */
  INT       lr,ls;
  DOUBLE    fac;
  DOUBLE    facr,facs;  /* integration weights */
  DOUBLE    det;
  DOUBLE    e1, e2;
  DOUBLE    preint;     /* pressure at integration point */
  DIS_TYP   typ;	/* element type */



  DOUBLE         visc;
  DOUBLE         DENS;

  STAB_PAR_GLS  *gls;	    /* pointer to GLS stabilisation parameters */
  FLUID_DYNAMIC *fdyn;



/******************************XFEM**************************************/
/******************************XFEM**************************************/
/******************************XFEM**************************************/
  ELEMENT       *myls2;
  INT            is_elcut;
  DOUBLE         lset01[4];
  INT           *ind;
  INT           *polygonmat[2];
  DOUBLE        *polygonwgt[2];
  DOUBLE        *polygonGP[2][2];
  LS_POLY_DATA  *polydata;
  INT            ntri;
  INT            nsub,nsubtotal;
  INT            index[8];
/******************************XFEM**************************************/
/******************************XFEM**************************************/
/******************************XFEM**************************************/

#ifdef DEBUG
  dstrc_enter("xfem_f2_calint");
#endif
/*----------------------------------------------------------------------*/

  /* initialization */
  iel    = ele->numnp;
  typ    = ele->distyp;

  gls    = ele->e.f2->stabi.gls;
  fdyn = alldyn[genprob.numff].fdyn;

  if (ele->e.f2->stab_type != stab_gls)
    dserror("routine with no or wrong stabilisation called");

  /* set index */
  for (i=0; i<iel; i++)
  {
    index[i] = 2*i;
    index[i+iel] = 3*iel+2*i;
  }
  /* set associated ls2 element */
  myls2 = ele->e.f2->my_ls;
  polydata = &(myls2->e.ls2->polydata[0]);
  /* is element cut? */
  if (xfem_data.xfem_on_off==1)
  {
    is_elcut = myls2->e.ls2->is_elcut;  /* enriched formulation! */
  }
  else
  {
    is_elcut = 0;                       /* standard formulation (with smeared material properties) */
  }
  /* access to the nodal values of level set profile */
  ls2_calset1(myls2,1,lset01);
  /* check! */
  if (is_elcut==1)
  {
    /* set ind */
    ind = polydata->ind;
    /* access to local coordinates of Gauss points */
    polygonGP[0][0] = polydata->polygonGP[0][0];
    polygonGP[1][0] = polydata->polygonGP[1][0];
    polygonGP[0][1] = polydata->polygonGP[0][1];
    polygonGP[1][1] = polydata->polygonGP[1][1];
    /* access to weights of Gauss points */
    polygonwgt[0] = polydata->polygonwgt[0];
    polygonwgt[1] = polydata->polygonwgt[1];
    /* access to material index */
    polygonmat[0] = polydata->polygonmat[0];
    polygonmat[1] = polydata->polygonmat[1];
    /* set icode and ihoel */
    if (typ==quad4 || typ== quad8 || typ== quad9) ihoel = 1;
    else if (typ== tri3 || typ==tri6) ihoel=0;
    else dserror("typ not set properly!");
    /**************************************/
    /* START LOOP OVER INTEGRATION POINTS */
    /**************************************/
    /* loop over triangles */
    for (ntri=0; ntri<2; ntri++)
    {
      /* check */
      if (ind[ntri]==-1)
      {
        continue;
      }
      else if (ind[ntri]==0)
      {
        nsubtotal = 1;
      }
      else
      {
        nsubtotal = 7;
      }
      /* loop over subtriangles */
      for (nsub=0; nsub<nsubtotal; nsub++)
      {
        /* set viscosity and density */
        actmat = polygonmat[ntri][nsub]-1;
        visc   = mat[actmat].m.fluid->viscosity;
        DENS   = mat[actmat].m.fluid->density;
        /* access to the local coordinates of the GP */
        e1   = polygonGP[ntri][0][nsub];
        facr = ONE;
        e2   = polygonGP[ntri][1][nsub];
        facs = polygonwgt[ntri][nsub];
        /* evaluate shape functions and their derivatives at Gauss point */
        xfem_f2_funct(funct,deriv,deriv2,e1,e2,typ,lset01,iel,is_elcut);
        /* Jacobian matrix */
        f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
        fac = facr*facs*det;
        /* compute global first derivatives */
        xfem_f2_derxy(derxy,deriv,xjm,det,iel,funct,lset01,is_elcut);
        /* compute velocities (n+1,i) at integration point */
        xfem_f2_veli(velint,funct,evelng,iel);
        /* compute velocity (n+1,i) derivatives at integration point */
        xfem_f2_vder(vderxy,derxy,evelng,iel);
        /*
           compute "Standard Galerkin" matrices
           NOTE =>
           Standard Galerkin matrices are stored in one matrix "estif"
           Standard Galerkin mass matrix is stored in "emass"
        */
        /* compute matrix Kvv */
        xfem_f2_calkvv(
          ele,estif,velint,NULL,vderxy,funct,derxy,
          fac,visc,iel,index,DENS
          );
        /* compute matrix Kvp and Kpv */
        xfem_f2_calkvp(
          estif,funct,derxy,fac,iel,index
          );
        /* compute matrix Mvv */
        if (fdyn->nis==0)
          xfem_f2_calmvv(
            emass,funct,fac,iel,index,DENS
            );
        /*
           compute Stabilization matrices
           NOTE =>
           Stabilization matrices are all stored in one matrix "estif"
           Stabilization mass matrices are stored in one matrix "emass"
        */
        /* compute stabilization parameter during integration loop */
        if (gls->iduring!=0)
          f2_calelesize2(ele,xyze,funct,velint,visc,iel,typ);
        /* compute second global derivative */
        if (ihoel!=0)
          xfem_f2_derxy2(
            xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel,funct,
            lset01,is_elcut
            );
          /* stabilization for matrix Kvv */
          xfem_f2_calstabkvv(
            ele,estif,velint,velint,NULL,vderxy,funct,
            derxy,derxy2,fac,visc,iel,ihoel,index,DENS
            );
          /* stabilization for matrix Kvp */
          xfem_f2_calstabkvp(
            ele,estif,velint,funct,derxy,derxy2,fac,
            visc,iel,ihoel,index,DENS
            );
          /* stabilization for matrix Mvv */
          if (fdyn->nis==0)
            xfem_f2_calstabmvv(
              ele,emass,velint,funct,derxy,derxy2,
              fac,visc,iel,ihoel,index,DENS
              );
          if (gls->ipres!=0)
          {
            /* stabilization for matrix Kpv */
            xfem_f2_calstabkpv(
              ele,estif,velint,NULL,vderxy,funct,derxy,
              derxy2,fac,visc,iel,ihoel,index,DENS
              );
            /* stabilization for matrix Mpv */
            if (fdyn->nis==0)
              xfem_f2_calstabmpv(
                emass,funct,derxy,fac,iel,index,DENS
                );
          }
        /* stabilization for matrix Kpp */
        if (gls->ipres!=0)
          f2_calstabkpp(
            gls,estif,derxy,fac,iel
            );

        /* compute "external" Force Vector (b) => */
        if (*hasext != 0)
        {
          /*  compute stabilisation part of external RHS (vel dofs) at (n+1) */
          xfem_f2_calstabexfv(
            ele,eiforce,derxy,derxy2,edeadng,velint,fac,visc,
            iel,ihoel,1,index,DENS
            );
          /*  compute stabilisation part of external RHS (pre dofs) at (n+1) */
          if (gls->ipres != 0)
            xfem_f2_calstabexfp(
              &(eiforce[2*iel]),derxy,edeadng,fac,iel,1,DENS
              );
        }

        /* compute "Time" Force Vector => */
        if (fdyn->nif!=0)
        {
          if (fdyn->iprerhs>0)
          {
            /* get pressure (n) at integration point */
            preint=f2_scali(
              funct,epren,iel
              );
            /* get pressure derivatives (n) at integration point */
            f2_pder(
              pderxy,derxy,epren,iel
              );
          }

          /* velocities (n) at integration point */
          xfem_f2_veli(
            velint,funct,eveln,iel
            );
          /* velocity (n) derivatives at integration point */
          xfem_f2_vder(
            vderxy,derxy,eveln,iel
            );
          /* get 2nd velocities derivatives (n) at integration point */
          if (ihoel!=0)
            xfem_f2_vder2(
              vderxy2,derxy2,eveln,iel
              );
          /* due to historical reasons there exist two velocities */
          vel2int=velint;
          /* get convective velocities (n) at integration point */
          f2_covi(
            vderxy,velint,covint
            );
          /* calculate galerkin part of "Time-RHS" (vel-dofs)*/
          xfem_f2_calgaltfv(
            etforce,vel2int,covint,funct,derxy,vderxy,
            preint,visc,fac,iel,index,DENS
            );
          /* calculate galerkin part of "Time-RHS" (pre-dofs) */
          f2_calgaltfp(
            &(etforce[2*iel]),funct,vderxy,fac,iel
            );
            /* calculate stabilization for "Time-RHS" (vel-dofs) */
          xfem_f2_calstabtfv(
            ele,etforce,velint,vel2int,covint,derxy,
            derxy2,vderxy,vderxy2,pderxy,fac,visc,ihoel,
            iel,index,DENS
            );
          /* calculate stabilization for "Time-RHS" (pre-dofs) */
          if (gls->ipres!=0)
            xfem_f2_calstabtfp(
              &(etforce[2*iel]),derxy,vderxy2,
              vel2int,covint,pderxy,visc,fac,ihoel,iel,DENS
              );

          if (*hasext!=0)
          {
            /*  compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
            xfem_f2_calgalexfv(
              etforce,funct,edeadn,edeadng,fac,iel,index,DENS
              );
              /*  compute stabilization part of external RHS (vel dofs) at (n) */
            xfem_f2_calstabexfv(
              ele,etforce,derxy,derxy2,edeadn,velint,fac,visc,
              iel,ihoel,0,index,DENS
              );
            /*  compute stabilization part of external RHS (pre dofs) */
            if (gls->ipres != 0)
              xfem_f2_calstabexfp(
                &(etforce[2*iel]),derxy,edeadn,fac,iel,0,DENS
                );
          }
        }
      }
    }
    /********************************************************************/
    /* END LOOP OVER INTEGRATION POINTS (MODIFIED INTEGRATION SCHEME)   */
    /********************************************************************/
  }
  else
  {
    /* set viscosity and density */
    actmat = ele->mat-1;
    visc   = mat[actmat].m.fluid->viscosity;
    DENS   = mat[actmat].m.fluid->density;
    /* get integration data and check if elements are "higher order" */
    switch (typ)
    {
        case quad4: case quad8: case quad9:
          ihoel   = 1;
          /* initialize integration */
          nir = ele->e.f2->nGP[0];
          nis = ele->e.f2->nGP[1];
          break;
        case tri3: case tri6:
          if (ele->e.f2->nGP[0]!=4) dserror("nGP not set properly!");
          ihoel = 0;
          /* initialize integration */
          nir = 1;
          nis = ele->e.f2->nGP[1];
          break;
        default:
          dserror("typ unknown!");
    }
    /********************************************************************/
    /* START LOOP OVER INTEGRATION POINTS (STANDARD INTEGRATION SCHEME) */
    /********************************************************************/
    for (lr=0; lr<nir; lr++)
    {
      for (ls=0; ls<nis; ls++)
      {
        /* get values of shape functions and their derivatives */
        switch(typ)
        {
            case quad4: case quad8: case quad9:   /* --> quad - element */
              e1   = data->qxg[lr][nir-1];
              facr = data->qwgt[lr][nir-1];
              e2   = data->qxg[ls][nis-1];
              facs = data->qwgt[ls][nis-1];
              xfem_f2_funct(
                funct,deriv,deriv2,e1,e2,typ,lset01,iel,
                is_elcut
                );
              break;
            case tri3: case tri6:
              e1   = data->txgr[ls][nis-1];
              facr = ONE;
              e2   = data->txgs[ls][nis-1];
              facs = data->twgt[ls][nis-1];
              xfem_f2_funct(
                funct,deriv,deriv2,e1,e2,typ,lset01,iel,
                is_elcut
                );
              break;
            default:
              dserror("typ unknown!");
        }
        /* Jacobian matrix */
        f2_jaco(
          xyze,funct,deriv,xjm,&det,iel,ele
          );
        fac = facr*facs*det;
        /* compute global derivatives */
        xfem_f2_derxy(
          derxy,deriv,xjm,det,iel,funct,lset01,is_elcut
          );
        /* get velocities (n+g,i) at integraton point */
        xfem_f2_veli(
          velint,funct,evelng,iel
          );
        /* get velocity (n+g,i) derivatives at integration point */
        xfem_f2_vder(
          vderxy,derxy,evelng,iel
          );
        /*
           compute "Standard Galerkin" matrices =>
           NOTE =>
           Standard Galerkin matrices are stored in one matrix "estif"
           Standard Galerkin mass matrix is stored in "emass"
        */
        /* compute matrix Kvv */
        xfem_f2_calkvv(
          ele,estif,velint,NULL,vderxy,funct,derxy,
          fac,visc,iel,index,DENS
          );
        /* compute matrix Kvp and Kpv */
        xfem_f2_calkvp(
          estif,funct,derxy,fac,iel,index
          );
        /* compute matrix Mvv */
        if (fdyn->nis==0)
          xfem_f2_calmvv(
            emass,funct,fac,iel,index,DENS
            );
        /*
          compute Stabilization matrices =>
          NOTE =>
          Stabilization matrices are all stored in one matrix "estif"
          Stabilization mass matrices are stored in one matrix "emass"
        */
          /* compute stabilization parameter during integration loop */
          if (gls->iduring!=0)
            f2_calelesize2(
              ele,xyze,funct,velint,visc,iel,typ
              );
          /* compute second global derivative */
          if (ihoel!=0)
            xfem_f2_derxy2(
              xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel,funct,
              lset01,is_elcut
              );

          /* stabilization for matrix Kvv */
          xfem_f2_calstabkvv(
            ele,estif,velint,velint,NULL,vderxy,
            funct,derxy,derxy2,fac,visc,iel,
            ihoel,index,DENS
            );
          /* stabilization for matrix Kvp */
          xfem_f2_calstabkvp(
            ele,estif,velint,funct,derxy,derxy2,
            fac,visc,iel,ihoel,index,DENS
            );
          /* stabilization for matrix Mvv */
          if (fdyn->nis==0)
            xfem_f2_calstabmvv(
              ele,emass,velint,funct,derxy,derxy2,
              fac,visc,iel,ihoel,index,DENS
              );
          if (gls->ipres!=0)
          {
            /* stabilization for matrix Kpv */
            xfem_f2_calstabkpv(
              ele,estif,velint,NULL,vderxy,funct,derxy,
              derxy2,fac,visc,iel,ihoel,index,DENS
              );
            /* stabilization for matrix Mpv */
            if (fdyn->nis==0)
              xfem_f2_calstabmpv(
                emass,funct,derxy,fac,iel,index,DENS
                );
          }

          /* stabilization for matrix Kpp */
          if (gls->ipres!=0)
            f2_calstabkpp(
              gls,estif,derxy,fac,iel
              );

        /* compute "external" Force Vector (b) => */
        if (*hasext != 0)
        {
          /*  compute stabilisation part of external RHS (vel dofs) at (n+1) */
          xfem_f2_calstabexfv(
            ele,eiforce,derxy,derxy2,edeadng,velint,fac,visc,
            iel,ihoel,1,index,DENS
            );
          /*  compute stabilisation part of external RHS (pre dofs) at (n+1) */
          if (gls->ipres != 0)
            xfem_f2_calstabexfp(
              &(eiforce[2*iel]),derxy,edeadng,fac,iel,1,DENS
              );
        }

        /* compute "Time" Force Vector => */
        if (fdyn->nif!=0)
        {
          if (fdyn->iprerhs>0)
          {
            /* get pressure (n) at integration point */
            preint=f2_scali(
              funct,epren,iel
              );
            /* get pressure derivatives (n) at integration point */
            f2_pder(
              pderxy,derxy,epren,iel
              );
          }

          /* get velocities (n) at integration point */
          xfem_f2_veli(
            velint,funct,eveln,iel
            );
          /* get velocitiederivatives (n) at integration point */
          xfem_f2_vder(
            vderxy,derxy,eveln,iel
            );
          /* get 2nd velocities derivatives (n) at integration point */
          if (ihoel!=0)
            xfem_f2_vder2(
              vderxy2,derxy2,eveln,iel
              );
          /* due to historical reasons there exist two velocities */
          vel2int=velint;
          /* get convective velocities (n) at integration point */
          f2_covi(
            vderxy,velint,covint
            );
          /* calculate galerkin part of "Time-RHS" (vel-dofs) */
          xfem_f2_calgaltfv(
            etforce,vel2int,covint,funct,derxy,vderxy,
            preint,visc,fac,iel,index,DENS
            );
          /* calculate galerkin part of "Time-RHS" (pre-dofs) */
          f2_calgaltfp(
            &(etforce[2*iel]),funct,vderxy,fac,iel
            );
            /* calculate stabilization for "Time-RHS" (vel-dofs) */
          xfem_f2_calstabtfv(
            ele,etforce,velint,vel2int,covint,derxy,
            derxy2,vderxy,vderxy2,pderxy,fac,visc,ihoel,iel,index,DENS
            );
          /* calculate stabilization for "Time-RHS" (pre-dofs) */
          if (gls->ipres!=0)
            xfem_f2_calstabtfp(
              &(etforce[2*iel]),derxy,vderxy2,
              vel2int,covint,pderxy,visc,fac,ihoel,iel,DENS
              );

          if (*hasext!=0)
          {
            /*  compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
            xfem_f2_calgalexfv(
              etforce,funct,edeadn,edeadng,fac,iel,index,DENS
              );
              /*  compute stabilization part of external RHS (vel dofs) at (n) */
            xfem_f2_calstabexfv(
             ele,etforce,derxy,derxy2,edeadn,velint,fac,visc,
             iel,ihoel,0,index,DENS
             );
           /*  compute stabilization part of external RHS (pre dofs) */
            if (gls->ipres != 0)
              xfem_f2_calstabexfp(
               &(etforce[2*iel]),derxy,edeadn,fac,iel,0,DENS
               );
          }
        }
      }
    }
    /********************************************************************/
    /* END LOOP OVER INTEGRATION POINTS (STANDARD INTEGRATION SCHEME)   */
    /********************************************************************/
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of xfem_f2_calint */
/*! @} (documentation module close)*/
#endif
