/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale_setdirich', 'ale_caldirich',
'ale_setdirich_increment', 'ale_caldirich_increment' and
'check_ale_dirich'

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"

#ifdef FSI_NONMATCH
#include "../fluid3/fluid3_prototypes.h"
#include "../shell8/shell8.h"
#endif

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
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief sets dirichlet boundary conditions on at time t

<pre>                                                              mn 06/02
This routine reads the initial value for the dirichlet condition from
actgnode->dirich->dirich_val.a.dv[j], gets the appropriate factor from
the timecurve and writes the value for the dirichlet conditions at the
time t to actnode->sol.a.da[0][j].

</pre>
\param *actfield  FIELD          (i)  my field
\param *adyn      ALE_DYNAMIK    (i)  structure containing time information
\param  readstructpos INT        (i)  position, where to read FSI structural
                                      displacement from

\warning For (dirich_val.a.dv == 90) the boundary conditions for a special
         example (rotating hole) are calculated.
\return void
\sa calling: dyn_facfromcurve(); called by: dyn_ale_lin()

*----------------------------------------------------------------------*/
void ale_setdirich(
    FIELD        *actfield,
    INT           disnum,
    ALE_DYNAMIC  *adyn,
    INT           readstructpos
    )
{
GNODE                *actagnode;
NODE                 *actsnode;
NODE                 *actanode;
NODE                 *actfnode;
ELEMENT              *actele;
INT                   i,j;
INT                   numnp_total;
INT                   numele_total;
INT                   actcurve;
INT                   numff,numsf;
INT                   dim;
DOUBLE                funct_fac;	             /* factor from spatial function */
DOUBLE                timefac[ALENUMTIMECURVE];
DOUBLE                T;
DOUBLE                dt;
DOUBLE                acttimefac;
DOUBLE                initval;
DOUBLE                delta, deltav[MAXDOFPERNODE];
DOUBLE                cx,cy,win,wino,dd;
ARRAY_POSITION       *ipos;
#ifdef FSI_NONMATCH
INT          m;          /*counters*/
ARRAY        funct_a;
DOUBLE       *funct;
DIS_TYP      typ;
DOUBLE       r,s,t;
ELEMENT      *hostele;
#endif

#ifdef DEBUG
dstrc_enter("ale_setdirich");
#endif

numnp_total  = actfield->dis[disnum].numnp;
numele_total = actfield->dis[disnum].numele;
T            = adyn->time;
dt           = adyn->dt;
numff        = genprob.numff;
numsf        = genprob.numsf;
dim          = genprob.ndim;
ipos   = &(actfield->dis[disnum].ipos);

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
}

/*------------------------------------------------- loop over all nodes */
for (i=0;i<numnp_total;i++)
{
   actanode  = &(actfield->dis[disnum].node[i]);
   actagnode = actanode->gnode;
   if (actagnode->dirich==NULL) continue;

    switch(actagnode->dirich->dirich_type)
    {
      case dirich_none:
        for (j=0;j<actanode->numdf;j++)
        {
          /*--------- to make sure that there's no garbage we zero the sol-array */
          actanode->sol_increment.a.da[ipos->dispnp][j] = ZERO;
          if (actagnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
          actcurve = actagnode->dirich->curve.a.iv[j]-1;
          if (actcurve<0)
            acttimefac = 1.0;
          else
            acttimefac = timefac[actcurve];
          funct_fac = actagnode->d_funct[j];
          initval  = actagnode->dirich->dirich_val.a.dv[j];
          /*=====================================================================*
            |    example: rotating hole (dirich_val.a.dv == 90)                   |
            |    sonst: Normalfall:                                               |
            |    actanode->sol.a.da[0][j] = initval*acttimefac;                   |
           *=====================================================================*/
          if (FABS(initval-90.0) > EPS13)
          {
            actanode->sol_increment.a.da[ipos->dispnp][j] = initval*acttimefac*funct_fac;
          }
          else
          {
            cx = actanode->x[0]-0.2;
            cy = actanode->x[1]-0.2;
            win = (initval * acttimefac * 3.14159265359)/180.0;
            if(FABS(cx) < EPS12) wino=3.14159265359/2.0;
            else
               wino= atan(cy/cx);
            dd = sqrt(cx*cx+cy*cy);
            if(cx < 0.0) wino += 3.14159265359;
            if (j==0)
            {
              actanode->sol_increment.a.da[ipos->dispnp][j] = dd * cos(win+wino) - cx;
            }
            else
            {
              actanode->sol_increment.a.da[ipos->dispnp][j] = dd * sin(win+wino) - cy;
            }
          }
        }
        break;

#ifdef D_FSI
   case dirich_FSI: /* dirichvalues = displacements of structure -------*/
#ifdef FSI_NONMATCH
/**************************************************************************/
/*The displacements are passed from the structure nodes to the ale nodes.
 *Since the position of all FSI ale nodes in relation to their structure host
 *elements is known (GNODE->coupleptr), the shape functions are used to
 *interpolate the ale displacements from the structure displacements. The
 *sum of all contributions from struct nodes delivers the nodal ale displacements.
 *           _            _          _
 *           u_ale = Sum( N(r,s,t) * u_struct )                             */
/**************************************************************************/

      if (adyn->coupmethod == 1) /* nonconfornming: Wurde hier einfach
                                  * auf 0 gesetzt!!! 0 heisst eigentlich: conforming*/
      {
        dsassert(actanode->locsysId==0,"no locsys for ALE dirich_FSI!\n");

        /*access the coupled structure host element*/
        hostele=actagnode->coupleptr->hostele;
        typ=hostele->distyp;

        /*initialize the ale displacemants with 0*/
        for (j=0;j<actanode->numdf;j++)
          actanode->sol_increment.a.da[ipos->dispnp][j]=0;

        switch (hostele->eltyp)
        {
#ifdef D_BRICK1
         case el_brick1:

         /*allocate space for evaluation of shape functions*/
         funct=amdef("funct" ,&funct_a ,MAXNOD_F3,1,"DV");
         /*access the local coordinates of the ale node*/
         r=actagnode->coupleptr->x[0];
         s=actagnode->coupleptr->x[1];
         t=actagnode->coupleptr->x[2];

         /*evaluate shape functions at the position of the ale node*/
         f3_hex(funct,NULL,NULL,r,s,t,typ,0);

         /*loop the dofs*/
         for (j=0;j<actanode->numdf;j++)
         {
           /*loop all structure nodes of the host element*/
           for (m=0;m<hostele->numnp;m++)
           {
             actsnode=hostele->node[m];
             /*calculate the ale displacemant by interpolation via the shape functions*/
             actanode->sol_increment.a.da[ipos->dispnp][j] += funct[m] * (actsnode->sol_mf.a.da[readstructpos][j]);
           }
         }
         amdel(&funct_a);
         break;
#endif
#ifdef D_SHELL8
        case el_shell8:
        {
          funct=amdef("funct" ,&funct_a ,MAXNOD_SHELL8,1,"DV");
          /*access the local coordinates of the ale node*/
          r=actagnode->coupleptr->x[0];
          s=actagnode->coupleptr->x[1];

          /*evaluate shape functions at the position of the ale node*/
          s8_funct_deriv(funct,NULL,r,s,typ,0);

          /*loop the dofs*/
          for (j=0;j<actanode->numdf;j++)
          {
            /*loop all structure nodes of the host element*/
            for (m=0;m<hostele->numnp;m++)
            {
              actsnode=hostele->node[m];
              /*calculate the ale displacemant by interpolation via the shape functions*/
              actanode->sol_increment.a.da[ipos->dispnp][j] += funct[m] * (actsnode->sol_mf.a.da[readstructpos][j]);
            }
          }
          amdel(&funct_a);
          break;
        }
#endif
         default:
           dserror("unknown element type %d", hostele->eltyp);
        }
      }
#else
      if(adyn->coupmethod == 1) /* conforming */
      {
        dsassert(actanode->locsysId==0,"no locsys for ALE dirich_FSI!\n");
        actsnode = actagnode->mfcpnode[numsf];
        for (j=0;j<actanode->numdf;j++)
        {
          actanode->sol_increment.a.da[ipos->dispnp][j] =
	    actsnode->sol_mf.a.da[readstructpos][j];

        } /* readstructpos = 0 for 'ordinary' calculation *
                           = 6 for calculation for Relaxation parameter via
	  		       steepest descent method */
      }
#endif
      else if(adyn->coupmethod == 0) /* mortar method */
      {
        for (j=0;j<actanode->numdf;j++)
        {
	   if (readstructpos != 6)
	   {
	     /* 'ordinary' calculation */
             actanode->sol_increment.a.da[ipos->dispnp][j] =
	       actanode->sol_mf.a.da[3][j];
	   }
           else if (readstructpos == 6)
	   {
	     /* steepest descent method */
             actanode->sol_increment.a.da[ipos->dispnp][j] =
	       actanode->sol_mf.a.da[3][j] -
	       actanode->sol_mf.a.da[ipos->mf_dispn][j];
	   }
           else
	     dserror("parameter readstructpos not known in ale_setdirich()!!!\n");
        } /* readstructpos = 0 for 'ordinary' calculation *
                           = 6 for calculation for Relaxation parameter via
			     steepest descent method */
      }

      else dserror("adyn->couptyp not known in ale_setdirich()!!!");
   break;
   case dirich_freesurf: /* dirichvalues = displacement of fluid
                            free surface                                */
      if (actanode->locsysId==0)
      {
         actfnode = actagnode->mfcpnode[numff];
         for (j=0;j<actanode->numdf;j++)
         {
            dsassert(actfnode->xfs != NULL, "actfnode->xfs[j] == NULL");
            delta = actfnode->xfs[j]-actanode->x[j];
            actanode->sol_increment.a.da[ipos->dispnp][j]=delta;
         }
      }
      else /* local co-system */
           /* NOTE: given is the position of the free surface from which
                    we get the prescribed displacement for this node
                    in the XYZ co-system. So the displacement vector
                    has to be tranformed to the xyz* co-system          */
      {
         actfnode = actagnode->mfcpnode[numff];
         actele = actanode->element[0];
         for (j=0;j<actanode->numdf;j++)
            deltav[j]=actfnode->xfs[j]-actanode->x[j];
         locsys_trans_nodval(actele,&(deltav[0]),actanode->numdf,
                             actanode->locsysId-1,0);
         for (j=0;j<actanode->numdf;j++)
            actanode->sol_increment.a.da[ipos->dispnp][j]=deltav[j];

      }
   break;

#endif
   default:
      dserror("dirich type unknown!\n");
   }
}

/*------------------------ dirichlet values are applied in the xyz* co-sys
   so transform nodal values with dirichlet conditions for
   sol_increment                                                          */
locsys_trans_sol_dirich(actfield,0,1,0,1);

/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of ale_setdirich*/


/*!----------------------------------------------------------------------
  \brief sets incremental dirichlet boundary conditions on at time t

  <pre>                                                             ck 06/03
  This routine reads the initial value for the dirichlet condition from
  actgnode->dirich->dirich_val.a.dv[j], gets the appropriate factor from
  the timecurve and writes the value for the incremental dirichlet conditions
  at the time t to actnode->sol_increment.a.da[0][j].
  Routine is used for fsi-ale problems with changing ale stiffness over time.

  </pre>
  \param *actfield  FIELD          (i)  my field
  \param *sdyn      STRUCT_DYNAMIK (i)  structure containing time information
  \param  actpos    INT            (i)  actual position in solution history

  \warning For (dirich_val.a.dv == 90) the boundary conditions for a special
  example (rotating hole) are calculated.
  \return void
  \sa calling: dyn_facfromcurve();
  called by: fsi_ale_nln(), fsi_ale_2step(), fsi_ale_spring(),
  fsi_ale_laplace()

 *----------------------------------------------------------------------*/
void ale_setdirich_increment_fsi(
    FIELD              *actfield,
    INT                 disnum,
    ALE_DYNAMIC        *adyn,
    INT                 actpos)
{
GNODE                *actagnode;
NODE                 *actsnode;
NODE                 *actanode;
NODE                 *actfnode;
INT                   i,j;
INT                   numnp_total;
INT                   numele_total;
INT                   actcurve;
INT                   diff,max;
INT                   numff,numsf;
INT                   dim;
DOUBLE     funct_fac;	             /* factor from spatial function */
DOUBLE                timefac[ALENUMTIMECURVE];
DOUBLE                T;
DOUBLE                dt;
DOUBLE                acttimefac;
DOUBLE                initval;
ARRAY_POSITION       *ipos;

#ifdef DEBUG
  dstrc_enter("ale_setdirich_increment_fsi");
#endif

numnp_total  = actfield->dis[disnum].numnp;
numele_total = actfield->dis[disnum].numele;
T            = adyn->time;
dt           = adyn->dt;
numff        = genprob.numff;
numsf        = genprob.numsf;
dim          = genprob.ndim;
ipos         = &(actfield->dis[disnum].ipos);

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T,&timefac[actcurve]) ;
}

/*------------------------------------------------- loop over all nodes */
for (i=0;i<numnp_total;i++)
{
   actanode  = &(actfield->dis[disnum].node[i]);
   dsassert(actanode->locsysId==0,"locsys not implemented yet!\n");
   actagnode = actanode->gnode;
   dsassert(actanode->locsysId==0,"incremental DBCs at ALE nodes!\n");
   if (actpos >= actanode->sol.fdim)
   {
      diff = actpos - actanode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actanode->sol),actanode->sol.fdim+max+1,actanode->sol.sdim,"DA");
    }
    if (actagnode->dirich==NULL)
      continue;

    switch(actagnode->dirich->dirich_type)
    {

      case dirich_none:
        for (j=0;j<actanode->numdf;j++)
        {
          if (actagnode->dirich->dirich_onoff.a.iv[j]==0)
            continue;
          actcurve = actagnode->dirich->curve.a.iv[j]-1;
          if (actcurve<0)
            acttimefac = 1.0;
          else
            acttimefac = timefac[actcurve];
          initval  = actagnode->dirich->dirich_val.a.dv[j];
          funct_fac = actagnode->d_funct[j];
          /* actanode->sol.a.da[actpos][j] = initval*acttimefac; */

          /*=====================================================================*
            |    example: rotating hole (dirich_val.a.dv == 90)                   |
            |    sonst: Normalfall:                                               |
            |    actanode->sol.a.da[0][j] = initval*acttimefac;                   |
           *=====================================================================*/
          if (FABS(initval-90.0) > EPS13)
          {
             actanode->sol_increment.a.da[ipos->dispnp][j] = initval*funct_fac*acttimefac
               - actanode->sol_increment.a.da[ipos->dispn][j];
          }
          else
          {
            DOUBLE cx, cy, win, wino, dd;
            cx = actanode->x[0]-0.2;
            cy = actanode->x[1]-0.2;
            win = (initval * acttimefac * 3.14159265359)/180.0;
            if(FABS(cx) < EPS12) wino=3.14159265359/2.0;
            else
               wino= atan(cy/cx);
            dd = sqrt(cx*cx+cy*cy);
            if(cx < 0.0) wino += 3.14159265359;
            if (j==0)
            {
              actanode->sol_increment.a.da[ipos->dispnp][j] = dd * cos(win+wino) - cx
               - actanode->sol_increment.a.da[ipos->dispn][j];
            }
            else
            {
              actanode->sol_increment.a.da[ipos->dispnp][j] = dd * sin(win+wino) - cy
               - actanode->sol_increment.a.da[ipos->dispn][j];
            }
          }
        }
        break;

#ifdef D_FSI
   case dirich_FSI: /* dirichvalues = displacements of structure -------*/
      actsnode = actagnode->mfcpnode[numsf];
      for (j=0;j<actanode->numdf;j++)
      {
         actanode->sol_increment.a.da[ipos->dispnp][j] =
	   actsnode->sol_mf.a.da[ipos->mf_dispn][j] -
	   actanode->sol_increment.a.da[ipos->dispn][j];
      }
   break;
   case dirich_freesurf: /* dirichvalues = displacement of fluid
                            free surface                                */
      actfnode = actagnode->mfcpnode[numff];
      for (j=0;j<actanode->numdf;j++)
      {
	actanode->sol_increment.a.da[ipos->dispnp][j] =
	  actfnode->xfs[j] - actanode->x[j] -
	  actanode->sol_increment.a.da[ipos->dispn][j];
      }
   break;
#endif

      default:
        dserror("dirich type unknown!\n");
    }
  }

#ifdef DEBUG
dstrc_exit();
#endif

  return;
} /* end of ale_setdirich_increment_fsi*/




/*!----------------------------------------------------------------------
  \brief sets dirichlet boundary conditions on at time t for incremental
  calculation

  <pre>                                                             ck 12/02
  This routine reads the initial value for the dirichlet condition from
  actgnode->dirich->dirich_val.a.dv[j], gets the appropriate factor from
  the timecurve at t and t - dt and writes the value for the dirichlet
  conditions at the time t to actnode->sol_increment.a.da[0][j].

  </pre>
  \param *actfield  FIELD          (i)  my field
  \param *sdyn      STRUCT_DYNAMIK (i)  structure containing time information

  \warning There is nothing special to this routine.
  \return void
  \sa calling: dyn_facfromcurve();
  called by: dyn_ale_nln(), dyn_ale_2step(), dyn_ale_spring(),
  dyn_ale_laplace()

 *----------------------------------------------------------------------*/
void ale_setdirich_increment(
    FIELD        *actfield,
    INT           disnum,
    ALE_DYNAMIC  *adyn
    )
{
  GNODE                *actgnode;
  NODE                 *actnode;
  INT                   i,j;
  INT                   numnp_total;
  INT                   numele_total;
  INT                   actcurve;
  DOUBLE                timefac[2][ALENUMTIMECURVE];
  DOUBLE                T;
  DOUBLE                T_prev;
  DOUBLE                acttimefac, prevtimefac;
  DOUBLE     funct_fac;	             /* factor from spatial function */
  DOUBLE                initval;

  DOUBLE                cx,cy,win,wino,winp,dd;


#ifdef DEBUG
  dstrc_enter("ale_setdirich_increment");
#endif

numnp_total  = actfield->dis[disnum].numnp;
numele_total = actfield->dis[disnum].numele;
T            = adyn->time;
T_prev       = T - adyn->dt;

/*------------------------------------------ get values from time curve */
for (actcurve=0;actcurve<numcurve;actcurve++)
{
  dyn_facfromcurve(actcurve,T     ,&timefac[0][actcurve]);
  dyn_facfromcurve(actcurve,T_prev,&timefac[1][actcurve]);
}

/*------------------------------------------------- loop over all nodes */
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[disnum].node[i]);
   dsassert(actnode->locsysId==0,"locsys not implemented yet!\n");
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL)
         continue;
   for (j=0;j<actnode->numdf;j++)
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0)
        continue;
      actcurve = actgnode->dirich->curve.a.iv[j]-1;
      if (actcurve<0)
      {
        acttimefac = 1.0;
        prevtimefac = 1.0;
      }
      else
      {
        acttimefac = timefac[0][actcurve];
        prevtimefac = timefac[1][actcurve];
      }
      initval  = actgnode->dirich->dirich_val.a.dv[j];
      funct_fac = actgnode->d_funct[j];
      /*=====================================================================*
        |    example: rotating hole (dirich_val.a.dv == 90)                   |
        |    sonst: Normalfall:                                               |
        |     actnode->sol.a.da[0][j] = initval*acttimefac;                   |
       *=====================================================================*/
      if (FABS(initval-90.0) > EPS13)
        actnode->sol_increment.a.da[0][j] = initval*funct_fac*(acttimefac-prevtimefac);
      else
      {
        cx = actnode->x[0];
        cy = actnode->x[1];
        win = (initval * acttimefac * 3.14159265359)/180.0;
        winp = (initval * prevtimefac * 3.14159265359)/180.0;
        if (cx != 0.0) wino= atan(cy/cx);
        else wino = 3.14159265359/2.;
        dd = sqrt(cx*cx+cy*cy);
        if(cx < 0.0) wino += 3.14159265359;
        if (j==0)
          actnode->sol_increment.a.da[0][j] = dd * cos(win+wino)
            - dd * cos(winp+wino);
        else
          actnode->sol_increment.a.da[0][j] = dd * sin(win+wino)
            - dd * sin(winp+wino);
      }
    }
    /*=====================================================================*/
    /*	 actnode->sol_increment.a.da[actpos][j] = initval*(acttimefac-prevtimefac); */
    /*=====================================================================*/
  } /* end loop over nodes */

#ifdef DEBUG
dstrc_exit();
#endif

  return;
} /* end of ale_setdirich_increment*/



/*!----------------------------------------------------------------------
  \brief calculates the element dirichlet load vector

  <pre>                                                              mn 06/02
  This routine calculates the element dirichlet load vector, reading the
  values of the dirichlet conditions from actnode->sol.a.da[0][j]

  </pre>
  \param *actele        ELEMENT (i)  my element
  \param *fullvec        DOUBLE  (o)  the dirichlet load vector as full vector
  \param dim            INT     (i)  dimension
  \param *estif_global  ARRAY   (i)  the element stiffness matrix

  \return void
  \sa calling: ---; called by: ale_rhs()

 *----------------------------------------------------------------------*/
void ale_caldirich(
                     ELEMENT   *actele,
		     DOUBLE    *fullvec,
		     INT        dim,
                     ARRAY     *estif_global
		    )
{

INT                   i,j;
INT                   numdf;
INT                   nd=0;
INT                   ilocsys;
DOUBLE              **estif;
DOUBLE                dirich[MAXDOFPERELE];
INT                   dirich_onoff[MAXDOFPERELE];
DOUBLE                val[MAXDOFPERNODE];
DOUBLE                dforces[MAXDOFPERELE];
GNODE                *actgnode;
NODE                 *actnode;
INT                   lm[MAXDOFPERELE];

#ifdef DEBUG
dstrc_enter("ale_caldirich");
#endif

/*----------------------------------------------------------------------*/
estif  = estif_global->a.da;
/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
   dforces[i] = 0.0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                                 dirichlet values at (n) were already */
/*                            written to the nodes (sol_icrement[0][j]) */
for (i=0; i<actele->numnp; i++)
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];
   actgnode = actnode->gnode;
   ilocsys=actnode->locsysId-1;
   if (ilocsys>=0) /* local co-sys */
   {
      /*------------------------- transform values at from XYZ to xyz* */
      for (j=0;j<numdf;j++)
         val[j]=actnode->sol_increment.a.da[0][j];
      locsys_trans_nodval(actele,&(val[0]),actnode->numdf,ilocsys,0);
      for (j=0; j<numdf; j++)
      {
         lm[i*numdf+j] = actele->node[i]->dof[j];
         if (actgnode->dirich==NULL) continue;
         dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
         dirich[i*numdf+j] = val[j];
      }
   }
   else
   {
      for (j=0; j<numdf; j++)
      {
         lm[i*numdf+j] = actele->node[i]->dof[j];
         if (actgnode->dirich==NULL) continue;
         dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
         dirich[i*numdf+j] = actnode->sol_increment.a.da[0][j];
      }
   }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */
/*-------- now assemble the vector dforces to the global vector fullvec */
for (i=0; i<nd; i++)
{
   if (lm[i] >= dim) continue;
   fullvec[lm[i]] += dforces[i];
}


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

  return;
} /* end of ale_caldirich*/




/*!----------------------------------------------------------------------
  \brief calculates the element dirichlet load vector for incremental
  calculations

  <pre>                                                            ck 05/03
  This routine calculates the element dirichlet load vector, reading the
  values of the dirichlet conditions from actnode->sol_increment.a.da[0][j]

  </pre>
  \param *actele        ELEMENT (i)  my element
  \param *fullvec       DOUBLE  (o)  the dirichlet load vector as full vector
  \param dim            INT     (i)  dimension
  \param *estif_global  ARRAY   (i)  the element stiffness matrix
  \param  place         INT     (i)  where to read in sol_increment

  \return void
  \sa calling: ---; calelem()

 *----------------------------------------------------------------------*/
void ale_caldirich_increment(
    ELEMENT   *actele,
    DOUBLE    *fullvec,
    INT        dim,
    ARRAY     *estif_global,
    INT        place
    )
{

  INT                   i,j;
  INT                   numdf;
  INT                   nd=0;
  DOUBLE              **estif;
  DOUBLE                dirich[MAXDOFPERELE];
  INT                   dirich_onoff[MAXDOFPERELE];
  DOUBLE                dforces[MAXDOFPERELE];
  GNODE                *actgnode;
  NODE                 *actnode;
  INT                   lm[MAXDOFPERELE];

#ifdef DEBUG
  dstrc_enter("ale_caldirich_increment");
#endif

  estif  = estif_global->a.da;

  /* set number of dofs on this element */
  for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

  /* init the vectors dirich and dirich_onoff */
  for (i=0; i<nd; i++)
  {
    dirich[i] = 0.0;
    dirich_onoff[i] = 0;
    dforces[i] = 0.0;
  }

  /* fill vectors dirich and dirich_onoff */
  /* dirichlet values at (n) were already written to the nodes (sol_increment[place][j]) */
  for (i=0; i<actele->numnp; i++)
  {
    numdf    = actele->node[i]->numdf;
    actnode  = actele->node[i];
    actgnode = actnode->gnode;
    for (j=0; j<numdf; j++)
    {
      lm[i*numdf+j] = actele->node[i]->dof[j];
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[i*numdf+j] = actnode->sol_increment.a.da[place][j];
    }
  }

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
      dforces[i] -= estif[i][j] * dirich[j];
    }/* loop j over columns */
  }/* loop i over rows */

  /* now assemble the vector dforces to the global vector fullvec */
  for (i=0; i<nd; i++)
  {
    if (lm[i] >= dim) continue;
    fullvec[lm[i]] += dforces[i];
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale_caldirich_increment*/


/*!----------------------------------------------------------------------
  \brief checks if element has node with Dirichlet condition

  <pre>                                                             ck 05/03
  This routine checks if element has Dirichlet condition.

  </pre>

 *-----------------------------------------------------------------------*/
INT check_ale_dirich(
    ELEMENT     *actele
    )
{
  INT j, k;
  INT hasdirich=0;
  GNODE *actgnode;

#ifdef DEBUG
  dstrc_enter("check_ale_dirich");
#endif


  for (j=0; j<actele->numnp; j++)
  {
    actgnode = actele->node[j]->gnode;
    if (actgnode->dirich==NULL)
      continue;
    else
    {
      for(k=0; k<actele->node[j]->numdf; k++)
      {
        /* if (actgnode->dirich->dirich_val.a.dv[k]!=0.0)*/
        if (actgnode->dirich != NULL)
        {
          hasdirich=1;
          goto end;
        }
      }
    }
  }

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
  dstrc_exit();
#endif

return hasdirich;

}

#endif
/*! @} (documentation module close)*/
