#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../xfem/xfem_prototypes.h"
#include "ls_prototypes.h"



extern struct _GENPROB       genprob;
extern struct _FIELD        *field;
extern struct _MATERIAL     *mat;
extern         ALLDYNA      *alldyn;



static ARRAY               lset00_a; 
static DOUBLE             *lset00;       /* element nodal values at (n)         */
static ARRAY               lset01_a; 
static DOUBLE             *lset01;       /* element nodal values at (n+1)       */
static ARRAY               lset02_a; 
static DOUBLE             *lset02;       /* element nodal values for sign comp. */
static ARRAY               velocity_a; 
static DOUBLE            **velocity;     /* elem. nodal velocities (from fluid) */
static ARRAY               funct_a;  
static DOUBLE             *funct;        /* shape functions                     */
static ARRAY               deriv_a;  
static DOUBLE            **deriv;        /* first natural derivatives           */
static ARRAY               xyze_a;
static DOUBLE            **xyze;   
static ARRAY               xjm_a;    
static DOUBLE            **xjm;          /* Jacobian matrix                     */
static ARRAY               derxy_a;  
static DOUBLE            **derxy;        /* global 1st derivatives              */
static ARRAY               id_a;
static INT               **id;   
static DOUBLE            **estif;        /* pointer to element tangent matrix   */
static DOUBLE            **emass;        /* pointer to element mass matrix      */
static DOUBLE             *etforce;      /* pointer to Time RHS                 */
static DOUBLE             *eiforce;      /* pointer to Iteration RHS            */
static INT                 algo;	 /* treatment of velocity (exp. or imp.)*/
static DOUBLE	           thdt;	 /* (  theta)*dt	  		*/
static DOUBLE              thdt1;	 /* (1-theta)*dt 			*/
static DOUBLE              epsilon;	 /* smoothing parameter			*/
static LS_DYNAMIC         *lsdyn;        /* pointer to dynamic control struct   */
static LS2_INTG_DATA      *data;         /* integration data                    */



static INT         nlayer;
static INT         is_elcut;



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2(
  PARTITION   *actpart,
  INTRA       *actintra,
  ELEMENT     *ele,             
  ARRAY       *estif_global,   
  ARRAY       *emass_global,   
  ARRAY       *etforce_global,
  ARRAY       *eiforce_global,
  CALC_ACTION *action
  )
{
#ifdef DEBUG 
  dstrc_enter("ls2");
#endif
/*----------------------------------------------------------------------*/
  
  switch (*action)
  {
      case calc_ls_init: /* initialize */
        estif   = estif_global->a.da;
        emass   = emass_global->a.da;
        eiforce = eiforce_global->a.dv;
        etforce = etforce_global->a.dv;
        
        ls2_init();
        
        break;
      case calc_ls: /* element contributions */  	   	               
        amzero(estif_global);
        amzero(emass_global);
        amzero(eiforce_global);
        amzero(etforce_global);
        
        ls2_calc(ele);
        
        break;
      default:
        dserror("action unknown\n");
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  
  return;
} /* end of ls2 */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_init()
{
#ifdef DEBUG 
  dstrc_enter("ls2_init");
#endif
/*----------------------------------------------------------------------*/

  lset00    = amdef("lset00", &lset00_a   , MAXNOD_LS2    , ONE           , "DV");
  lset01    = amdef("lset01", &lset01_a   , MAXNOD_LS2    , ONE           , "DV");
  lset02    = amdef("lset02", &lset02_a   , MAXNOD_LS2    , ONE           , "DV");
  funct     = amdef("funct" , &funct_a    , TWO*MAXNOD_LS2, ONE           , "DV");
  deriv     = amdef("deriv" , &deriv_a    , TWO           , MAXNOD_LS2    , "DA");
  xjm       = amdef("xjm"   , &xjm_a      , TWO           , TWO           , "DA");
  xyze      = amdef("xyze"  , &xyze_a     , TWO           , MAXNOD_LS2    , "DA");
  derxy     = amdef("derxy" , &derxy_a    , TWO           , MAXNOD_LS2    , "DA");
  id        = amdef("id"    , &id_a       , TWO           , TWO           , "IA");  
  velocity  = amdef("vel"   , &velocity_a , TWO           , TWO*MAXNOD_LS2, "DA");
  
  if (genprob.numls==-1)
    dserror("\n**ERROR** there is no levelset field initialized!\n");
  
 end:
  /* set lsdyn */
  lsdyn = alldyn[genprob.numls].lsdyn;

  /* smoothing parameter for the computation of sign function */
  epsilon = lsdyn->lsdata->epsilon;

  /* set identity array */
  id[0][0] = 1;
  id[0][1] = 0;
  id[1][0] = 0;
  id[1][1] = 1;

  data = (LS2_INTG_DATA*)CCACALLOC(1,sizeof(LS2_INTG_DATA));
  /* set integration data */
  ls2_intg(data);
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_init */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_calc(
  ELEMENT* ele
  )
{
#ifdef DEBUG 
  dstrc_enter("ls2_calc");
#endif
/*----------------------------------------------------------------------*/
  
  /* set nlayer */
  nlayer = ele->e.ls2->nlayer;
  /* set is_elcut */
  is_elcut = ele->e.ls2->is_elcut;
  
  if (lsdyn->lsdata->reinitflag==1)
    ls2_calc_reinitialized(ele);
  else
  {
    if (lsdyn->lsdata->localization==1 && nlayer==0)
      ls2_calc_localized(ele);
    else
      ls2_calc_nonlocalized(ele);
  }
  
/*----------------------------------------------------------------------*/
 end: 
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_calc */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_calc_nonlocalized(
  ELEMENT* ele
  )
{
  INT         i,j,k,l;
  INT         iel;
  INT         nir,nis;
  INT         lr,ls;
  DOUBLE      fac;
  DOUBLE      facr,facs;
  DOUBLE      det;
  DOUBLE      e1,e2;
  INT         ntyp;
  DIS_TYP     typ;
  
  DOUBLE      area=0.0;
  DOUBLE      tau=0.0; /* stabilization parameter */
  DOUBLE      velGP[2]={0.0,0.0};
  DOUBLE      velnorm=0.0;
  
#ifdef DEBUG 
  dstrc_enter("ls2_calc_nonlocalized");
#endif
/*----------------------------------------------------------------------*/

  /* set integration parameters */
  thdt  = (      lsdyn->theta)*lsdyn->dt;               /* (  theta)*dt */
  thdt1 = (1.0 - lsdyn->theta)*lsdyn->dt;		/* (1-theta)*dt */
  
  /* set element parameters */
  iel = ele->numnp;
  typ = ele->distyp;
  ntyp = ele->e.ls2->ntyp;
  /*
   * NOTE =>
   * ntyp = 1 for quadrilateral element
   * ntyp = 2 for triangular element
   */

  /*
   * get nodal coordinates of the element and nodal values of the levelset 
   * function at time step (n) and (n+1)
   */
  ls2_calset(ele,xyze,lset00,lset01);

  /* set integration data */
  switch (ntyp)
  {
      case 1: /* => quadrilateral element */
        nir = ele->e.ls2->nGP[0];
        nis = ele->e.ls2->nGP[1];
        break;
      case 2: /* => triangular element */
        nir = 1;
        nis = ele->e.ls2->nGP[0];
        break;
      default:
        dserror("typ unknown!");
  }
  /* set the velocity field */
  ls2_setfluidvel(ele);
  /* compute the area */
  if (lsdyn->lsdata->isstab==1)
  {
    for (lr=0; lr<nir; lr++)
    {
      for (ls=0; ls<nis; ls++)
      {
        /* shape functions and their derivatives */
        switch(ntyp)  
        {
            case 1:
              e1   = data->xgq[lr][nir-1];
              facr = data->wgtq[lr][nir-1];
              e2   = data->xgq[ls][nis-1];
              facs = data->wgtq[ls][nis-1];
              ls2_funct(funct,deriv,e1,e2,typ);
              break;
            case 2:
              e1   = data->xgtr[ls][nis-1];
              facr = ONE;
              e2   = data->xgts[ls][nis-1];
              facs = data->wgtt[ls][nis-1];
              ls2_funct(funct,deriv,e1,e2,typ);
              break;
            default:
              dserror("typ unknown!");
        }
        /* Jacobian matrix */
        ls2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
        fac = facr*facs*det;
        area += fac;      
      }
    }

  }
  /* loop over integration points*/
  for (lr=0; lr<nir; lr++)
  {
    for (ls=0; ls<nis; ls++)
    {
      /* shape functions and their derivatives */
      switch(ntyp)  
      {
          case 1:
            e1   = data->xgq[lr][nir-1];
            facr = data->wgtq[lr][nir-1];
            e2   = data->xgq[ls][nis-1];
            facs = data->wgtq[ls][nis-1];
            ls2_funct(funct,deriv,e1,e2,typ);
            break;
          case 2:
            e1   = data->xgtr[ls][nis-1];
            facr = ONE;
            e2   = data->xgts[ls][nis-1];
            facs = data->wgtt[ls][nis-1];
            ls2_funct(funct,deriv,e1,e2,typ);
            break;
          default:
            dserror("typ unknown!");
      }
      /* Jacobian matrix */
      ls2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;
      /* global derivates */
      ls2_gder(derxy,deriv,xjm,det,iel);
      /* velocity at Gauss point */
      xfem_f2_funct1(funct,e1,e2,typ,iel,lset00,is_elcut);
      velGP[0] = 0.0;
      velGP[1] = 0.0;
      for (i=0; i<2; i++)
        for (j=0; j<TWO*iel; j++)
          velGP[i] += funct[j]*velocity[i][j];
      /* stabilization parameter */
      tau = 0.0;
      if (lsdyn->lsdata->isstab==1)
      {
        velnorm  = 0.0;
        /* norm of the velocity */
        velnorm = velGP[0]*velGP[0] + velGP[1]*velGP[1];
        velnorm = sqrt(velnorm);
        /* stabilization parameter */
        tau = 0.0;
        if (velnorm!=0.0) tau = 0.5*sqrt(area)/velnorm;
      }
      /*
       * => MASS matrix
       */
      /*
       * contribution due to Galerkin part
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          emass[i][j] += funct[i]*funct[j]*fac;     
      /*
       * contribution due to stabilization
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            emass[i][j] += tau*velGP[k]*derxy[k][i]*funct[j]*fac;
      /*
      /* => TANGENT matrix
      /*
       * contribution due to Galerkin part
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            estif[i][j] += thdt*funct[i]*velGP[k]*derxy[k][j]*fac;
      /*
       * contribution due to stabilization
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              estif[i][j] += thdt*tau*velGP[k]*velGP[l]*derxy[k][i]*derxy[l][j]*fac;
    }
  }
  /* TIME right hand side */
  for (i=0; i<iel; i++)
    for (j=0; j<iel; j++)
      etforce[i] += emass[i][j]*lset00[j] - thdt1/thdt*estif[i][j]*lset00[j];
  /* ITERATION right hand side */
  for (i=0; i<iel; i++)
    for (j=0; j<iel; j++)
      eiforce[i] += -emass[i][j]*lset01[j] -  estif[i][j]*lset01[j];
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_calc_nonlocalized */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_calc_reinitialized(
  ELEMENT* ele
  )
{
  INT         i,j,k,l;
  INT         iel;
  INT         nir,nis;
  INT         lr,ls;
  DOUBLE      fac;
  DOUBLE      facr,facs;
  DOUBLE      det;
  DOUBLE      e1,e2;
  INT         ntyp;
  DIS_TYP     typ;

  DOUBLE      area=0.0;
  DOUBLE      tau=0.0; /* stabilization parameter */
  DOUBLE      velnorm=0.0;
  DOUBLE      normal[2]={0.0,0.0};  
  DOUBLE      velGP[2]={0.0,0.0};
  DOUBLE      gradnorm=0.0;
  DOUBLE      phi=0.0;
  DOUBLE      sign=0.0;
  DOUBLE      fs[4]={0.0,0.0,0.0,0.0};
 
#ifdef DEBUG 
  dstrc_enter("ls2_calc_reinitialized");
#endif
/*----------------------------------------------------------------------*/

  /* set integration parameters */
  thdt  = (      lsdyn->theta)*lsdyn->lsdata->rdt;      /* (  theta)*dt */
  thdt1 = (1.0 - lsdyn->theta)*lsdyn->lsdata->rdt;	/* (1-theta)*dt */

  /* set element parameters */
  iel = ele->numnp;
  typ = ele->distyp;
  ntyp = ele->e.ls2->ntyp;
  /*
   * get nodal coordinates of the element and nodal values of the levelset 
   * function at time step (n) and (n+1)
   */
  ls2_calset(ele,xyze,lset00,lset01);
  /* set lset02 */
  ls2_calset1(ele,2,lset02);
  /* integration data */
  switch (ntyp)
  {
      case 1:
        nir = ele->e.ls2->nGP[0];
        nis = ele->e.ls2->nGP[1];
        break;
      case 2:
        nir = 1;
        nis = ele->e.ls2->nGP[0];
        break;
      default:
        dserror("typ unknown!");
  }
  /* compute the area */
  if (lsdyn->lsdata->isstab==1)
  {
    for (lr=0;lr<nir;lr++)
    {
      for (ls=0;ls<nis;ls++)
      {
        /* shape functions and their derivatives */
        switch(ntyp)  
        {
            case 1:
              e1   = data->xgq[lr][nir-1];
              facr = data->wgtq[lr][nir-1];
              e2   = data->xgq[ls][nis-1];
              facs = data->wgtq[ls][nis-1];
              ls2_funct(funct,deriv,e1,e2,typ);
              break;
            case 2:
              e1   = data->xgtr[ls][nis-1];
              facr = ONE;
              e2   = data->xgts[ls][nis-1];
              facs = data->wgtt[ls][nis-1];
              ls2_funct(funct,deriv,e1,e2,typ);
              break;
            default:
              dserror("typ unknown!");
        }
        /* Jacobian matrix */
        ls2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
        fac = facr*facs*det;
        area += fac;      
      }
    }
  }
  /* loop over integration points */
  for (lr=0; lr<nir; lr++)
  {
    for (ls=0; ls<nis; ls++)
    {
      /* shape functions and their derivatives */
      switch(ntyp)  
      {
          case 1:
            e1   = data->xgq[lr][nir-1];
            facr = data->wgtq[lr][nir-1];
            e2   = data->xgq[ls][nis-1];
            facs = data->wgtq[ls][nis-1];
            ls2_funct(funct,deriv,e1,e2,typ);
            break;
          case 2:
            e1   = data->xgtr[ls][nis-1];
            facr = ONE;
            e2   = data->xgts[ls][nis-1];
            facs = data->wgtt[ls][nis-1];
            ls2_funct(funct,deriv,e1,e2,typ);
            break;
          default:
            dserror("typ unknown!");
      }
      /* Jacobian matrix */
      ls2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;
      /* global derivates */
      ls2_gder(derxy,deriv,xjm,det,iel);
      /* velocity at the Gauss point */
      velGP[0] = 0.0; velGP[1] = 0.0;
      normal[0] = 0.0; normal[1] = 0.0;
      /* norm of gradient of the level set field */
      gradnorm = 0.0;
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            gradnorm += derxy[k][i]*derxy[k][j]*lset00[i]*lset00[j];
      gradnorm = sqrt(gradnorm);
      /* normal */
      for (i=0; i<2; i++)
        for (j=0; j<iel; j++)
          normal[i] += derxy[i][j]*lset00[j]/gradnorm;
      /* compute smoothed sign function */
      sign = 0;
      phi = 0;
      for (i=0; i<iel; i++)
        phi += funct[i]*lset02[i];
      sign = phi/sqrt(phi*phi+epsilon*epsilon);
      /* velocity field */
      for (i=0; i<2; i++)
        velGP[i] = sign*normal[i];
      /* stabilization parameter */
      tau = 0.0;
      if (lsdyn->lsdata->isstab==1)
      {
        velnorm  = 0.0;
        /* norm of the velocity */
        velnorm = velGP[0]*velGP[0] + velGP[1]*velGP[1];
        velnorm = sqrt(velnorm);
        /* stabilization parameter */
        tau = 0.0;
        if (velnorm!=0.0) tau = 0.5*sqrt(area)/velnorm;
      }
      /*
      /* MASS matrix */
      /*
       * contribution due to Galerkin part
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          emass[i][j] += funct[i]*funct[j]*fac;     
      /*
       * contribution due to stabilization
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            emass[i][j] += tau*velGP[k]*derxy[k][i]*funct[j]*fac;
      /*
      /* TANGENT matrix */
      /*
       * contribution due to Galerkin part
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            estif[i][j] += thdt*funct[i]*velGP[k]*derxy[k][j]*fac;
      /*
       * contribution due to stabilization
       */
      for (i=0; i<iel; i++)
        for (j=0; j<iel; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              estif[i][j] += thdt*tau*velGP[k]*velGP[l]*derxy[k][i]*derxy[l][j]*fac;      
      /*
      /* SOURCE term
       */
      /*
       * contribution due to Galerkin part
       */      
      for (i=0; i<iel; i++)
        fs[i] -= sign*funct[i]*fac;
      /*
       * contribution due to stabilization
       */
        for (i=0; i<iel; i++)
          for (j=0; j<2; j++)
            fs[i] -= tau*sign*velGP[j]*derxy[j][i]*fac;
    }
  }
  /* TIME right hand side */
  for (i=0; i<iel; i++)
    for (j=0; j<iel; j++)
      etforce[i] += emass[i][j]*lset00[j] - thdt1/thdt*estif[i][j]*lset00[j]
                                          - thdt1*fs[i];
  /* ITERATION right hand side */
  for (i=0; i<iel; i++)
    for (j=0; j<iel; j++)
      eiforce[i] += -emass[i][j]*lset01[j] - estif[i][j]*lset01[j]
                                           - thdt*fs[i];        
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_calc_reinitialized */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_calc_localized(
  ELEMENT* ele
  )
{
  INT     i,j;
  INT     iel;
  NODE   *actnode;

#ifdef DEBUG 
  dstrc_enter("ls2_calc_localized");
#endif
/*----------------------------------------------------------------------*/

  /* check whether element is active */
  if (nlayer!=0) dserror("\nelement is in active region");
  
  iel = ele->numnp;
  /* initialize */
  for (i=0; i<iel; i++)
    for (j=0; j<iel; j++)
    {
      emass[i][j] = 0.0;
      estif[i][j] = 0.0;
    }
  /* stiffness matrix */
  for (i=0; i<iel; i++)
    estif[i][i] = 1.0;
  /* time right hand side */
  for (i=0; i<iel; i++)
    etforce[i] = 0.0;
  /* iteration right hand side */
  for (i=0; i<iel; i++)
    eiforce[i] = 0.0;
  /* check whether element has active nodes */
  for (i=0; i<iel; i++)
  {
    actnode = ele->node[i];
    if (actnode->gnode->is_node_active==1 ||
        actnode->gnode->is_node_active==2) estif[i][i] = 0.0;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_calc_localized */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_setfluidvel(
  ELEMENT* ele
  )
{
#ifdef DEBUG 
  dstrc_enter("ls2_setfluidvel");
#endif
/*----------------------------------------------------------------------*/
  
  switch (lsdyn->lsdata->setvel)
  {
      case 1:
        /* by user */
        ls2_setfluidvel_byuser(ele->numnp,lsdyn->lsdata->flag_vel);
        break;
      case 2:
        /* by fluid field */
        ls2_setfluidvel_byfluid(ele);
        break;
      default:
        dserror("action unknown\n");
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_setfluidvel */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_setfluidvel_byuser(
  INT iel,
  INT flag
  )
{
  INT         i,j;
  DOUBLE      dx1,dx2,dxl;
  DOUBLE      rad,rad2;
  DOUBLE      height = 1.0;
  
#ifdef DEBUG 
  dstrc_enter("ls2_setfluidvel_byuser");
#endif
/*----------------------------------------------------------------------*/
  
  switch (flag)
  {
      case 1: /* uniform velocity in the horizontal direction */
        for (i=0; i<iel;i++ )
        {
          velocity[0][i] = 1.0;
          velocity[1][i] = 0.0;
        }
        break;
      case 2: /* uniform velocity in the vertical direction */
        for (i=0; i<iel;i++ )
        {
          velocity[0][i] = 0.0;
          velocity[1][i] = 1.0;
        }
        break;
      case 3: /* uniform velocity in the direction of diagonal */
        for (i=0; i<iel;i++ )
        {
          velocity[0][i] = 1.0;
          velocity[1][i] = 0.5;
        }
        break;
      case 4: /* uniform velocity in the direction of diagonal */
        for (i=0; i<iel;i++ )
        {
          velocity[0][i] = 1.0;
          velocity[1][i] = -0.5;
        }
        break;
      case 5: /* modify the velocity field so that we have an expanding circle */
        for (i=0; i<iel; i++)
        {
          dx1 = xyze[0][i]-lsdyn->lsdata->xc1;
          dx2 = xyze[1][i]-lsdyn->lsdata->yc1;
          dxl = sqrt(dx1*dx1+dx2*dx2);
          velocity[0][i] = dx1/dxl;
          velocity[1][i] = dx2/dxl;
        }
        break;
      case 6: /* modify the velocity field so that we have a rotating circle */
        for (i=0; i<iel; i++)
        {
          dx1 = xyze[0][i]-height/2.0;
          dx2 = xyze[1][i]-height/2.0;
          dxl = sqrt(dx1*dx1+dx2*dx2);
          velocity[0][i] = dx2/(height/2.0);
          velocity[1][i] = -dx1/(height/2.0);
        }
        break;
      case 7: /* modify the velocity field so that we have linearly varying
                 velocity field along the height */
        for (i=0; i<iel; i++)
        {
          velocity[0][i] = xyze[1][i]/height;
          velocity[1][i] = 0.0;
        }
        break;
      case 8: /* modify the velocity field so that we have quadratically varying
                 velocity field along the height */
        for (i=0; i<iel; i++)
        {
          velocity[0][i] = -xyze[1][i]*(xyze[1][i]-height)/(height/2.0)/(height/2.0);
          velocity[1][i] = 0.0;
        }
        break;
      default:
        dserror("action unknown\n");
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_setfluidfield_byuser */



/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls2_setfluidvel_byfluid(
  ELEMENT* ele
  )
{
  INT          i;
  INT          iel;
  ELEMENT     *my_fluid;
  NODE	      *actnode;
  
#ifdef DEBUG 
  dstrc_enter("ls2_setfluidvel_byfluid");
#endif
/*----------------------------------------------------------------------*/

  /* corresponding fluid element */
  my_fluid = ele->e.ls2->my_fluid;
  iel = my_fluid->numnp;
  
  for (i=0; i<iel; i++)
  {
    actnode = my_fluid->node[i];
    /*
     * NOTE =>
     * to convect the level set profile use the velocity field
     * corresponding to time t(n+1)(velocity field obtained by
     * combining the values to t(n) and t(n+1) is also possible)
     */
    velocity[0][    i] = actnode->sol_increment.a.da[3][0];
    velocity[1][    i] = actnode->sol_increment.a.da[3][1];
    velocity[0][iel+i] = actnode->sol_increment.a.da[3][3];
    velocity[1][iel+i] = actnode->sol_increment.a.da[3][4];
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls2_setfluidvel_byfluid */
#endif
