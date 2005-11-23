/*-----------------------------------------------------------------------*/
/*!
\file
\brief service routines for the fast fluid3 element

  Very detailed description.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid3/fluid3.h"

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

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of material laws                                              |
  | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


static FLUID_DYNAMIC   *fdyn;



/*-----------------------------------------------------------------------*/
/*!
  \brief set element coordinates

  The nodal coordinates of the elements are copied into the vector elecord
  that can be used in fortran.

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param elecord         DOUBLE   (o) vector containing nodal coordinates
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalelecord(
    ELEMENT     *ele[LOOPL],
    DOUBLE      *elecord,
    INT          sizevec[6]
    )
{
  INT l,j;

#ifdef DEBUG
  dstrc_enter("f3_calelecord");
#endif

  for (l=0; l<sizevec[1]; l++)
  {
    for(j=0;j<sizevec[4];j++)
    {
      elecord[                   LOOPL*l+j]=ele[j]->node[l]->x[0];
      elecord[  sizevec[0]*LOOPL+LOOPL*l+j]=ele[j]->node[l]->x[1];
      elecord[2*sizevec[0]*LOOPL+LOOPL*l+j]=ele[j]->node[l]->x[2];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}




/*-----------------------------------------------------------------------*/
/*!
  \brief set all arrays for element calculation

  get the element velocities and the pressure at different times

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param eveln          *DOUBLE   (o) vels at time n
  \param evelng         *DOUBLE   (o) vels at time n+g
  \param epren          *DOUBLE   (o) pres at time n
  \param edeadn         *DOUBLE   (o) dead load at n
  \param edeadng        *DOUBLE   (o) dead load at n+g
  \param ipos                     (i) node array positions
  \param hasext         *INT      (i) flag for external loads
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalset(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *ehist,
    DOUBLE          *evelng,
    DOUBLE          *epren,
    DOUBLE          *edeadn,
    DOUBLE          *edeadng,
    ARRAY_POSITION  *ipos,
    INT             *hasext,
    INT              sizevec[6]
    )
{
  INT i,j;
  INT    actmat  ;    /* material number of the element */
  INT    velnp, hist; /* position flags in sol_increment field */
  DOUBLE dens;        /* density */
  NODE  *actnode;     /* actual node */
  GVOL  *actgvol;



#ifdef DEBUG
  dstrc_enter("f3_calset");
#endif

  fdyn    = alldyn[genprob.numff].fdyn;
  velnp   = ipos->velnp;
  hist    = ipos->hist;

/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[ipos][i]: solution at some time level           |
 |       ipos-flags:  ipos->velnp ... solution at time (n+1)            |
 |                    ipos->veln  ... solution at time (n)              |
 |                    ipos->velnm ... solution at time (n-1)            |
 |                    ipos->gridv ... mesh velocity in actual time step |
 |                    ipos->convn ... convective velocity at time (n)   |
 |                    ipos->convnp... convective velocity at time (n+1) |
 |                    ipos->hist  ... sol. data for rhs                 |
 *---------------------------------------------------------------------*/

  for(i=0;i<sizevec[1];i++) /* loop nodes */
  {

    for(j=0;j<sizevec[4];j++)
    {

      actnode=ele[j]->node[i];

      /* set element velocities (n+gamma) */
      evelng[LOOPL*i+j]=actnode->sol_increment.a.da[velnp][0];
      evelng[sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[velnp][1];
      evelng[2*sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[velnp][2];

      ehist[                   LOOPL*i+j]=actnode->sol_increment.a.da[hist][0];
      ehist[  sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[hist][1];
      ehist[2*sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[hist][2];
      /* pressure at time (n+1) */
      epren[LOOPL*i+j]=actnode->sol_increment.a.da[velnp][3];
    } /*loop*/

  } /* end of loop over nodes */


  /* check for dead load */
  actgvol = ele[0]->g.gvol;
  if (actgvol->neum!=NULL)
  {
    actmat=ele[0]->mat-1;
    dens = mat[actmat].m.fluid->density;
    for (i=0;i<3;i++)
    {
      if (actgvol->neum->neum_onoff.a.iv[i]==0)
      {

        for(j=0;j<sizevec[4];j++)
        {

          edeadn[LOOPL*i+j]  = ZERO;
          edeadng[LOOPL*i+j] = ZERO;
        }/*loop*/

      }

      if (actgvol->neum->neum_type==neum_dead  &&
          actgvol->neum->neum_onoff.a.iv[i]!=0)
      {

        for(j=0;j<sizevec[4];j++)
        {
          actgvol = ele[j]->g.gvol;

          edeadn[LOOPL*i+j]  = actgvol->neum->neum_val.a.dv[i]*dens;
          edeadng[LOOPL*i+j] = actgvol->neum->neum_val.a.dv[i]*dens;
          (hasext[j])++;
        }/*loop*/

      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_calset */




#ifdef D_FSI

/*-----------------------------------------------------------------------*/
/*!
  \brief set all arrays for element calculation

  get the element velocities and the pressure at different times

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param eveln          *DOUBLE   (o) vels at time n
  \param evelng         *DOUBLE   (o) vels at time n+g
  \param ealecovn       *DOUBLE   (o) ALE-convective vels at time n
  \param ealecovng      *DOUBLE   (o) ALE-convective vels at time n+g
  \param egridv         *DOUBLE   (o) grid velocity
  \param epren          *DOUBLE   (o) pres at time n
  \param edeadn         *DOUBLE   (o) dead load at n
  \param edeadng        *DOUBLE   (o) dead load at n+g
  \param ipos                     (i) node array positions
  \param hasext         *INT      (i) flag for external loads
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalseta(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *ehist,
    DOUBLE          *evelng,
    DOUBLE          *ealecovng,
    DOUBLE          *egridv,
    DOUBLE          *epren,
    DOUBLE          *edeadn,
    DOUBLE          *edeadng,
    ARRAY_POSITION  *ipos,
    INT             *hasext,
    INT              sizevec[6]
    )
{
  INT i,j;
  INT    actcurve;    /* actual time curve */
  INT    velnp, hist, convnp, convn, gridv;
  DOUBLE acttimefac;  /* time factor from actual curve */
  DOUBLE acttimefacn; /* time factor at time (n) */
  DOUBLE acttime;
  NODE  *actnode;     /* actual node */
  GVOL  *actgvol;



#ifdef DEBUG
  dstrc_enter("f3fcalseta");
#endif

  fdyn    = alldyn[genprob.numff].fdyn;
  velnp   = ipos->velnp;
  hist    = ipos->hist;
  convnp  = ipos->convnp;
  convn   = ipos->convn;
  gridv   = ipos->gridv;

  /*---------------------------------------------------------------------*
    | position of the different solutions:                                |
    | node->sol_incement: solution history used for calculations          |
    |     sol_increment[ipos->flag][i]: solution at (n-1)                  |
    |     ipos->flgas: velnp  ... velocity at time (n+1)                   |
    |                 veln   ... velocity at time (n)                     |
    |                 hist   ... old solution data needed for rhs         |
    |                 convnp ... convective velocity at time (n+1)        |
    |                 convn  ... convective velocity at time (n)          |
    |                 gridv  ... grid velocity within actual time step    |
   *---------------------------------------------------------------------*/


  for(i=0;i<sizevec[1];i++) /* loop nodes */
  {

    for(j=0;j<sizevec[4];j++)
    {

      actnode=ele[j]->node[i];

      /* set element velocities (n+gamma) */
      evelng[LOOPL*i+j]=actnode->sol_increment.a.da[velnp][0];
      evelng[sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[velnp][1];
      evelng[2*sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[velnp][2];

      ealecovng[LOOPL*i+j]=actnode->sol_increment.a.da[convnp][0];
      ealecovng[sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[convnp][1];
      ealecovng[2*sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[convnp][2];

      egridv[LOOPL*i+j]=actnode->sol_increment.a.da[gridv][0];
      egridv[sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[gridv][1];
      egridv[2*sizevec[0]*LOOPL+LOOPL*i+j]=actnode->sol_increment.a.da[gridv][2];

      ehist[                   LOOPL*i+j]= actnode->sol_increment.a.da[hist][0];
      ehist[  sizevec[0]*LOOPL+LOOPL*i+j]= actnode->sol_increment.a.da[hist][1];
      ehist[2*sizevec[0]*LOOPL+LOOPL*i+j]= actnode->sol_increment.a.da[hist][2];
      /* pressure at time (n+1) */
      epren[LOOPL*i+j]=actnode->sol_increment.a.da[velnp][3];

    }/*loop*/

  } /* end of loop over nodes */


  /* check for dead load */
  actgvol = ele[0]->g.gvol;
  if (actgvol->neum!=NULL)
  {
    if (actgvol->neum->neum_type==neum_LAS)
    {
      actcurve = actgvol->neum->curve-1;
      if (actcurve<0)
        dserror("No Time curve given for neum_LAS!\n");
      acttime = fdyn->acttime;
      dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
      acttime=fdyn->acttime-fdyn->dta;
      dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;

      for(j=0;j<sizevec[4];j++)
      {
         edeadn[LOOPL*0+j] = actgvol->neum->neum_val.a.dv[0]*acttimefacn;
        edeadng[LOOPL*0+j] = actgvol->neum->neum_val.a.dv[0]*acttimefac;
         edeadn[LOOPL*1+j] = ZERO;
        edeadng[LOOPL*1+j] = ZERO;
         edeadn[LOOPL*2+j] = actgvol->neum->neum_val.a.dv[2];
        edeadng[LOOPL*2+j] = actgvol->neum->neum_val.a.dv[2];
        (hasext[j])++;
      }/*loop*/
    }

    else
    {
      actcurve = actgvol->neum->curve-1;
      if (actcurve<0)
      {
        acttimefac =ONE;
        acttimefacn=ONE;
      }
      else
      {
        acttime = fdyn->acttime;
        dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
        acttime = fdyn->acttime-fdyn->dta;
        dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
      }

      for (i=0;i<3;i++)
      {
        if (actgvol->neum->neum_onoff.a.iv[i]==0)
        {

          for(j=0;j<sizevec[4];j++)
          {

            edeadn[LOOPL*i+j]  = ZERO;
            edeadng[LOOPL*i+j] = ZERO;
          }/*loop*/

        }
        if (actgvol->neum->neum_type==neum_dead  &&
            actgvol->neum->neum_onoff.a.iv[i]!=0)
        {

          for(j=0;j<sizevec[4];j++)
          {
            actgvol = ele[j]->g.gvol;

            edeadn[LOOPL*i+j]  = actgvol->neum->neum_val.a.dv[i]*acttimefacn;
            edeadng[LOOPL*i+j] = actgvol->neum->neum_val.a.dv[i]*acttimefac;
            (hasext[j])++;
          }/*loop*/

        }
      }
    }
  }


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3fcalseta */

#endif /* ifdef D_FSI */



/*-----------------------------------------------------------------------*/
/*!
  \brief routine to calculate element size and stabilisation parameter

  Very detailed description what the function is doing why.

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param funct          *DOUBLE   (i) shape functions
  \param deriv          *DOUBLE   (i) deriv. of shape funcs
  \param deriv2         *DOUBLE   (i) 2nd deriv. of sh. funcs
  \param derxy          *DOUBLE   (i) global derivatives
  \param xjm            *DOUBLE   (i) jacobian matrix
  \param evel           *DOUBLE   (i) element velocities
  \param velint         *DOUBLE   (i) vel at int point
  \param wa1            *DOUBLE   (o) working array
  \param elecord         DOUBLE   (i) vector containing nodal coordinates
  \param tau            *DOUBLE   (o) stabilisation parameter
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalelesize(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *funct,
    DOUBLE          *deriv,
    DOUBLE          *deriv2,
    DOUBLE          *derxy,
    DOUBLE          *xjm,
    DOUBLE          *evel,
    DOUBLE          *velint,
    DOUBLE          *wa1,
    DOUBLE          *elecord,
    DOUBLE          *tau,
    INT              sizevec[6]
    )
{
  INT ieval = 0;             /* evaluation flag */
  INT l,ilen, inod;          /* simply a counter */
  INT istrnint;              /* evaluation flag */
  INT ishvol;                /* evaluation flag */
  INT actmat;                /* number of actual material */
  DOUBLE visc;               /* fluid viscosity */
  DOUBLE det[LOOPL];         /* determinant of jacobian */
  DOUBLE vol[LOOPL];         /* element volume */
  DOUBLE val[LOOPL];         /* temporary calculation value */
  DOUBLE velno;              /* velocity norm */
  DOUBLE strle[LOOPL];       /* streamlength */
  DOUBLE e1,e2,e3;           /* natural coordinates of inegration point */
  DOUBLE fac[LOOPL],facr=0.0;/* factors */
  DOUBLE facs=0.0,fact=0.0;  /* factors */
  DOUBLE velino[3];          /* normed velocity vector at integration point */
  DIS_TYP typ;
  STAB_PAR_GLS *gls;         /* pointer to GLS stabilisation parameters */
  FLUID_DYNAMIC *fdyn;
  FLUID_DATA    *data;


  /*fortran variables */
  INT typint;
  INT icode;

#ifdef DEBUG
  dstrc_enter("f3_calelesize");
#endif

  /* initialise */
  icode=2;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  typ    = ele[0]->distyp;
  gls    = ele[0]->e.f3->stabi.gls;

  if (ele[0]->e.f3->stab_type != stab_gls)
    dserror("routine with no or wrong stabilisation called");

  istrnint = gls->istrle * gls->ninths;
  ishvol  = fdyn->ishape * gls->iareavol;

  /*----------------------------------------------------------------------*
    | calculations at element center: area & streamlength                  |
    | NOTE:                                                                |
    |    volume is always calculated using only 1 integrationpoint         |
    |    --> it may be possible to save some operations here by replacing  |
    |         e1,e2,e3, facr,facs,fact with their constant values in the   |
    |         calls of f3_hex / f3_tet!!!!!!                               |
   *----------------------------------------------------------------------*/

  if (ishvol==1)
  {
    for(l=0;l<sizevec[4];l++)
    {
      vol[l]   = ZERO;
      strle[l] = ZERO;
    }/*loop*/

    /*get values of integration parameters, shape functions and their derivatives */
    switch(typ)
    {
      case hex8: case hex20: case hex27:   /* hex - element */
        e1   = data->qxg[0][0];
        facr = data->qwgt[0][0];
        e2   = data->qxg[0][0];
        facs = data->qwgt[0][0];
        e3   = data->qxg[0][0];
        fact = data->qwgt[0][0];

        if(typ==hex8) typint=8;
        else if(typ==hex20) typint=20;
        else if(typ==hex27) typint=27;

        f3fhex(funct, deriv, deriv2, &e1, &e2, &e3, &typint,
            &icode, sizevec);
        break;


      case tet4: case tet10:  /* tet - element */
        e1   = data->txgr[0][0];
        facr = data->twgt[0][0];
        e2   = data->txgs[0][0];
        facs = ONE;
        e3   = data->txgs[0][0];
        fact = ONE;

        if(typ==tet4) typint=4;
        else if(typ==tet10) typint=10;

        f3ftet(funct, deriv, deriv2, &e1, &e2, &e3, &typint,
            &icode, sizevec);
        break;

      default:
        dserror("typ unknown!");
    } /*end switch(typ) */

    ieval++;

    /*  compute jacobian matrix */
    f3fjaco(funct, deriv, xjm, det, elecord, sizevec);

    for(l=0;l<sizevec[4];l++)
    {
      fac[l]=facr*facs*fact*det[l];
      vol[l] += fac[l];
    } /*loop*/


    if (istrnint==1)    /* compute streamlength */
    {
      f3fveli(velint, funct, evel, sizevec);
      f3fgder(derxy, deriv, xjm, wa1, det, sizevec);

      ieval++;

      for(l=0;l<sizevec[4];l++)
      {

        val[l] = ZERO;
        velno=sqrt( velint[l]*velint[l] \
            + velint[LOOPL+l]*velint[LOOPL+l] \
            + velint[2*LOOPL+l]*velint[2*LOOPL+l]);
        if(velno>=EPS6)
        {
          velino[0] = velint[l]/velno;
          velino[1] = velint[LOOPL+l]/velno;
          velino[2] = velint[2*LOOPL+l]/velno;
        }
        else
        {
          velino[0] = ONE;
          velino[1] = ZERO;
          velino[2] = ZERO;
        }
        for (inod=0;inod<sizevec[1];inod++) /* loop element nodes */
        {
          val[l] += FABS(velino[0]*derxy[                   inod*LOOPL+l]
                        +velino[1]*derxy[  sizevec[0]*LOOPL+LOOPL*inod+l]
                        +velino[2]*derxy[2*sizevec[0]*LOOPL+LOOPL*inod+l]);
        } /* end of loop over elements */
        strle[l]=TWO/val[l];

      } /*loop*/

    } /* endif (istrnint==1) */


    /* set element sizes *
       loop over 3 different element sizes: vel/pre/cont  */
    for(l=0;l<sizevec[4];l++)
    {
      for(ilen=0;ilen<3;ilen++)
      {
        if (gls->ihele[ilen]==1)
          ele[l]->e.f3->hk[ilen] = pow(vol[l],(ONE/THREE));
        else if (gls->ihele[ilen]==2)
          ele[l]->e.f3->hk[ilen] = pow((SIX*vol[l]/PI),(ONE/THREE));
        else if (gls->ihele[ilen]==3)
          ele[l]->e.f3->hk[ilen] = pow((SIX*vol[l]/PI),(ONE/THREE))/sqrt(THREE);
        else if (gls->ihele[ilen]==4)
          dserror("ihele[i] = 4: calculation of element size not possible!!!");
        else if (gls->ninths==1)
          ele[l]->e.f3->hk[ilen] = strle[l];
      } /* end of loop over ilen */
    }/*loop*/
  } /* endif (ishvol==1) */


  /*----------------------------------------------------------------------*
    | calculations at element center: only streamlength                    |
    | NOTE:                                                                |
    |    volume is always calculated using only 1 integrationpoint         |
    |    --> it may be possible to save some operations here by replacing  |
    |         e1,e2,e3, facr,facs,fact with their constant values in the   |
    |         calls of f3_hex / f3_tet!!!!!!                               |
   *----------------------------------------------------------------------*/
  else if (istrnint==1 && ishvol !=1)
  {
    for(l=0;l<sizevec[4];l++)
    {
      vol[l]   = ZERO;
      strle[l] = ZERO;
    } /*loop*/

    /* get values of integration parameters, shape functions and their derivatives */
    switch(typ)
    {
      case hex8: case hex20: case hex27:   /* hex - element */
        e1   = data->qxg[0][0];
        facr = data->qwgt[0][0];
        e2   = data->qxg[0][0];
        facs = data->qwgt[0][0];
        e3   = data->qxg[0][0];
        fact = data->qwgt[0][0];

        if(typ==hex8) typint=8;
        else if(typ==hex20) typint=20;
        else if(typ==hex27) typint=27;

        f3fhex(funct, deriv, deriv2, &e1, &e2, &e3, &typint,
            &icode, sizevec);

        break;

      case tet4: case tet10:  /* tet - element */
        e1   = data->txgr[0][0];
        facr = data->twgt[0][0];
        e2   = data->txgs[0][0];
        facs = ONE;
        e3   = data->txgs[0][0];
        fact = ONE;

        if(typ==tet4) typint=4;
        else if(typ==tet10) typint=10;

        f3ftet(funct, deriv, deriv2, &e1, &e2, &e3, &typint,
            &icode, sizevec);

        break;

      default:
        dserror("typ unknown!");
    } /*end switch(typ) */

    /* compute jacobian matrix */
    f3fjaco(funct, deriv, xjm, det, elecord, sizevec);

    /* compute stream length */
    f3fveli(velint,funct,evel,sizevec);
    f3fgder(derxy, deriv, xjm, wa1, det, sizevec);

    ieval++;

    for(l=0;l<sizevec[4];l++)
    {

      val[l] = ZERO;
      velno=sqrt( velint[l]*velint[l] \
          + velint[LOOPL+l]*velint[LOOPL+l] \
          + velint[2*LOOPL+l]*velint[2*LOOPL+l]);
      if(velno>=EPS6)
      {
        velino[0] = velint[l]/velno;
        velino[1] = velint[LOOPL+l]/velno;
        velino[2] = velint[2*LOOPL+l]/velno;
      }
      else
      {
        velino[0] = ONE;
        velino[1] = ZERO;
        velino[2] = ZERO;
      }
      for (inod=0;inod<sizevec[1];inod++) /* loop element nodes */
      {
        val[l] += FABS(velino[0]*derxy[                   LOOPL*inod+l]
                      +velino[1]*derxy[  sizevec[0]*LOOPL+LOOPL*inod+l]
                      +velino[2]*derxy[2*sizevec[0]*LOOPL+LOOPL*inod+l]);
      } /* end of loop over element nodes */
      strle[l]=TWO/val[l];

      /* set element sizes *
         loop over 3 different element sizes: vel/pre/cont */
      for (ilen=0;ilen<3;ilen++)
      {
        if (gls->ihele[ilen]==5)
          ele[l]->e.f3->hk[ilen] = strle[l];
      } /* end of loop over ilen */

    }/*loop*/

  } /* endif (istrnint==1 && ishvol !=1) */



  /* calculate stabilisation parameter */
  if(gls->istapc==1 || istrnint==1)
  {
    switch(ieval) /* ival>2: vel at intpoint already available! */
    {
      case 0:
        /* get only values of integration parameters and shape functions
           no derivatives */
        switch(typ)
        {
          case hex8: case hex20: case hex27:    /* quad - element */
            e1   = data->qxg[0][0];
            facr = data->qwgt[0][0];
            e2   = data->qxg[0][0];
            facs = data->qwgt[0][0];
            e3   = data->qxg[0][0];
            fact = data->qwgt[0][0];

            if(typ==hex8) typint=8;
            else if(typ==hex20) typint=20;
            else if(typ==hex27) typint=27;

            f3fhex(funct, deriv, deriv2, &e1, &e2, &e3, &typint,
                &icode, sizevec);
            break;

          case tet4: case tet10:       /* tet - element */
            e1   = data->txgr[0][0];
            facr = data->twgt[0][0];
            e2   = data->txgs[0][0];
            facs = ONE;
            e3   = data->txgs[0][0];
            fact = ONE;

            if(typ==tet4) typint=4;
            else if(typ==tet10) typint=10;

            f3ftet(funct, deriv, deriv2, &e1, &e2, &e3, &typint,
                &icode, sizevec);

            break;

          default:
            dserror("typ unknown!");
        } /* end switch (typ) */

        f3fveli(velint,funct,evel,sizevec);
        break;


      case 1:
        f3fveli(velint,funct,evel,sizevec);
        break;


      case 2:
        break;


      default:
        dserror("wrong value for ieval");
    } /* end switch (ieval) */

    /* calculate stabilisation parameter */
    actmat=ele[0]->mat-1;
    visc = mat[actmat].m.fluid->viscosity;
    f3fcalstabpar(ele,velint,tau,visc,typint,-1,sizevec);

  } /* endif (ele->e.f3->istapc==1 || istrnint==1) */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_calelesize */



/*-----------------------------------------------------------------------*/
/*!
  \brief routine to calculate element size and stabilisation parameter

  in this routine the element size and the stabilisation parameter
  is calculated for one element during the integration loop

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param velint         *DOUBLE   (i) vel at int point
  \param derxy          *DOUBLE   (i) global derivatives
  \param tau            *DOUBLE   (o) stabilisation parameter
  \param inttyp          INT      (i) typ of integration
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalelesize2(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *velint,
    DOUBLE          *derxy,
    DOUBLE          *tau,
    DOUBLE           visc,
    INT              inttyp,
    INT              sizevec[6]
    )
{
  INT    l,ilen, inod;     /* simply a counter */
  INT    istrnint;         /* evaluation flag */
  DOUBLE strle[LOOPL];     /* stream length */
  DOUBLE val[LOOPL];       /* temporary calculation value */
  DOUBLE velno;            /* velocity norm */
  DOUBLE velino[3];        /* normed velocity vector at integration point */
  STAB_PAR_GLS *gls;       /* pointer to GLS stabilisation parameters */


#ifdef DEBUG
  dstrc_enter("f3_calelesize2");
#endif

  /* initialise */
  gls      = ele[0]->e.f3->stabi.gls;
  istrnint = gls->istrle * gls->ninths;

  if (ele[0]->e.f3->stab_type != stab_gls)
    dserror("routine with no or wrong stabilisation called");

  for(l=0;l<sizevec[4];l++)
  {
    val[l] = ZERO;
  }

  if (istrnint==2)
  {
    /* compute streamlength */
    for(l=0;l<sizevec[4];l++)
    {
      velno=sqrt( velint[l]*velint[l] \
          + velint[LOOPL+l]*velint[LOOPL+l] \
          + velint[2*LOOPL+l]*velint[2*LOOPL+l]);
      if(velno>=EPS6)
      {
        velino[0] = velint[l]/velno;
        velino[1] = velint[LOOPL+l]/velno;
        velino[2] = velint[2*LOOPL+l]/velno;
      }
      else
      {
        velino[0] = ONE;
        velino[1] = ZERO;
        velino[2] = ZERO;
      }
      for (inod=0;inod<sizevec[1];inod++)
      {
        val[l] += FABS(velino[0]*derxy[                   inod*LOOPL+l]
                      +velino[1]*derxy[  sizevec[0]*LOOPL+LOOPL*inod+l]
                      +velino[2]*derxy[2*sizevec[0]*LOOPL+LOOPL*inod+l]);
      }
      strle[l]=TWO/val[l];

      /* set element sizes *
         loop over 3 different element sizes: vel/pre/cont */
      for (ilen=0;ilen<3;ilen++)
      {
        if (gls->ihele[ilen]==5)
          ele[l]->e.f3->hk[ilen] = strle[l];
      } /* end of loop over ilen */

    } /*loop*/

  } /* endif (istrnint==2) */

  /* calculate stabilisation parameter */
  f3fcalstabpar(ele,velint,tau,visc,inttyp,1,sizevec);

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_calelesize2 */



/*-----------------------------------------------------------------------*/
/*!
  \brief routine to calculate the element dirichlet load vector

  in this routine the element load vector due to dirichlet conditions
  is calcluated. The prescribed values are taken from the node solution
  history at (n+1) 'dirich[j] = actnode->sol_increment.a.da[ipos->velnp][j]'.
  the element load vector 'dforce' is calculated by eveluating
  </pre>
  \code
  dforces[i] -= estif[i][j] * dirich[j];
  \endcode

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param dforces        *DOUBLE   (o) dirichlet force vector
  \param estif          *DOUBLE   (i) element stiffness matrix
  \param ipos                     (i) node array positions
  \param hasdirich      *INT      (i) flag if s.th. was written to dforces
  \param isrelax         INT      (i) flag for relaxation
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcaldirich(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *dforces,
    DOUBLE          *estif,
    ARRAY_POSITION  *ipos,
    INT             *hasdirich,
    INT              is_relax,
    INT              sizevec[6]
    )
{

  INT         i,j,k;
  INT         nrow;
  INT         numdf;                      /* number of fluid dofs */
  INT         nd=0;
  DOUBLE      dirich[MAXDOFPERELE];       /* dirichlet values of act. ele */
  INT         dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele */
  GNODE      *actgnode;                   /* actual GNODE */
  NODE       *actnode;                    /* actual NODE */


#ifdef DEBUG
  dstrc_enter("f3fcaldirich");
#endif

  /* check if there are any dirichlet conditions *
     for the nodes of this element */
  for(k=0;k<sizevec[4];k++)
  {

    for (i=0; i<ele[0]->numnp; i++)
    {
      actgnode = ele[k]->node[i]->gnode;
      if (actgnode->dirich==NULL)
        continue;
      else
        hasdirich[k]=1;
      break;
    } /* end loop over nodes */

    if (hasdirich[k]==0) /* --> no nodes with DBC for this element */
      continue;

    /* set number of dofs on this element */
    nd=0;
    for (i=0; i<ele[0]->numnp; i++) nd += ele[k]->node[i]->numdf;

    /* init the vectors dirich and dirich_onoff */
    for (i=0; i<nd; i++)
    {
      dirich[i] = 0.0;
      dirich_onoff[i] = 0;
    }

    /* fill vectors dirich and dirich_onoff
       dirichlet values at (n+1) were already
       written to the nodes (sol_increment[ipos->velnp][j]) */
    nrow=0;
    for (i=0; i<ele[0]->numnp; i++) /* loop nodes */
    {
      numdf    = ele[k]->node[i]->numdf;
      actnode  = ele[k]->node[i];
      actgnode = actnode->gnode;
      for (j=0; j<numdf; j++) /* loop dofs */
      {
        if (actgnode->dirich==NULL) continue;
        dirich_onoff[nrow+j] = actgnode->dirich->dirich_onoff.a.iv[j];
        if (is_relax) /* calculation for relax.-Param. reads
                         dbc from (sol_increment[ipos->relax][j]) */
          dirich[nrow+j] = actnode->sol_increment.a.da[ipos->relax][j];
        else
          dirich[nrow+j] = actnode->sol_increment.a.da[ipos->velnp][j];
      } /* end loop over dofs */
      nrow+=numdf;
    } /* end loop over nodes */
    dsassert(nrow==nd,"failure during calculation of dirich forces\n");

    /* loop rows of element matrix */
    for (i=0; i<nd; i++)
    {
      /* do nothing for supported row */
      if (dirich_onoff[i]!=0) continue;
      /* loop columns of unsupported row */
      for (j=0; j<nd; j++)
      {
        /* do nothing for unsupported columns */
        if (dirich_onoff[j]==0) continue;
        dforces[i*LOOPL+k] -= estif[i*sizevec[2]*LOOPL+LOOPL*j+k] * dirich[j];
      }/* loop j over columns */
    }/* loop i over rows */


  }/*loop -- k*/


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of f3fcaldirich*/



/*-----------------------------------------------------------------------*/
/*!
  \brief add emass and estif to estif

  Very detailed description what the function is doing why.

  \param estif          *DOUBLE   (i/o) element stiffness matrix
  \param emass          *DOUBLE   (i) element mass matrix
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fmake_estif(
    DOUBLE            *estif,
    DOUBLE            *emass,
    INT                sizevec[6]
    )

{
  INT i,j,k,rows,l;
  DOUBLE thsl;

  fdyn    = alldyn[genprob.numff].fdyn;

  k    = 4*sizevec[1];

  rows = sizevec[2];

  thsl = fdyn->thsl;

  for(i=0;i<k;i++)
    for(j=0;j<k;j++)
      for(l=0;l<sizevec[4];l++)
        estif[i*rows*LOOPL+LOOPL*j+l] *= thsl;

  if (fdyn->nis==0)
    for(i=0;i<k;i++)
      for(j=0;j<k;j+=4)
        for(l=0;l<sizevec[4];l++)
        {
          estif[i*rows*LOOPL+LOOPL*j    +l] += emass[i*rows*LOOPL+LOOPL*j    +l];
          estif[i*rows*LOOPL+LOOPL*(j+1)+l] += emass[i*rows*LOOPL+LOOPL*(j+1)+l];
          estif[i*rows*LOOPL+LOOPL*(j+2)+l] += emass[i*rows*LOOPL+LOOPL*(j+2)+l];
        }

}/* End of make_estif*/



#ifdef D_FSI

/*-----------------------------------------------------------------------*/
/*!
  \brief set element coordinates during ALE calculations

  nodal coordinates of actual element are evaluated. since nodes at the
  free surface have changing coordinates during the nonlin. iteration
  one has to treat them separately.

  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param elecord         DOUBLE   (o) vector containing nodal coordinates
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3falecord(
    ELEMENT     *ele[LOOPL],
    DOUBLE      *elecord,
    INT          sizevec[6]
    )

{

  INT    l,i;
  DOUBLE xy;
  NODE  *actfnode;    /* actual fluid node */
  NODE  *actanode;    /* actual ale node */
  GNODE *actfgnode;   /* actual fluid gnode */


#ifdef DEBUG
  dstrc_enter("f3falecord");
#endif


  /* set element coordinates */
  for (l=0; l<sizevec[1]; l++)
  {
    for(i=0;i<sizevec[4];i++)
    {
      actfnode  = ele[i]->node[l];
      actfgnode = actfnode->gnode;
      actanode  = actfgnode->mfcpnode[genprob.numaf];

      xy     = actfnode->x[0];
      elecord[                   LOOPL*l+i]=xy + actanode->sol_mf.a.da[1][0];

      xy     = actfnode->x[1];
      elecord[  sizevec[0]*LOOPL+LOOPL*l+i]=xy + actanode->sol_mf.a.da[1][1];

      xy     = actfnode->x[2];
      elecord[2*sizevec[0]*LOOPL+LOOPL*l+i]=xy + actanode->sol_mf.a.da[1][2];

    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

return;
} /* end of f3falecord */

#endif


#endif /* ifdef D_FLUID3_F */


/*! @} (documentation module close)*/


