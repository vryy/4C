/*!----------------------------------------------------------------------
\file
\brief contains the routines 'wge_stiff_de' -> calculates stiffness
       matrix K_de as sum over gaussian points,
*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*!
\addtogroup WALLGE
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief routine 'wge_stiff_de' -> calculates stiffness
       matrix K_de as sum over gaussian points

*----------------------------------------------------------------------*/
void wge_stiff_de(DOUBLE  **Kde,      /* stiffnes displ - equiv,strain */
                  DOUBLE  **bopd,     /* B-operator for displacements  */
                  DOUBLE   *E,        /* Material tangent vector       */
                  DOUBLE   *functe,   /* Anatzfunc. for equiv.strain   */
                  DOUBLE    fac,      /* integration factor            */
                  INT       numdfd,   /* element displacement DOF      */
                  INT       iele,     /* element equiv.strain DOF      */
                  INT       neps)     /* number of strain components   */
{
#ifdef D_WALLGE

INT            i,j,k,m;
DOUBLE         dum;
DOUBLE         EN[3];
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_enter("wge_stiff_de");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<iele; j++)
{
  for (k=0; k<neps; k++)
  {
   EN[k] = E[k]*functe[j]*fac;
  }
  for (i=0; i<numdfd; i++)
  {
    dum = 0.0;
    for (m=0; m<neps; m++)
    {
     dum = dum + bopd[m][i]*EN[m];
    }
    Kde[i][j] = Kde[i][j] + dum;
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_stiff_de */


/*!----------------------------------------------------------------------
\brief routine 'wge_stiff_ed' -> calculates stiffness
       matrix K_ed as sum over gaussian points

*----------------------------------------------------------------------*/
void wge_stiff_ed(DOUBLE  **Ked,      /* stiffnes equiv,strain - disp  */
                  DOUBLE  **bopd,     /* B-operator for displacements  */
                  DOUBLE   *F,        /* Material tangent vector       */
                  DOUBLE   *functe,   /* Anatzfunc. for equiv.strain   */
                  DOUBLE    fac,      /* integration factor            */
                  INT       numdfd,   /* element displacement DOF      */
                  INT       iele,     /* element equiv.strain DOF      */
                  INT       neps)     /* number of strain components   */
{
#ifdef D_WALLGE

INT            i,j,k;
DOUBLE         dum,FB;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("wge_stiff_ed");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<numdfd; j++)
{
  FB = 0.0;
  for (k=0; k<neps; k++)
  {
   FB = FB + F[k] * bopd[k][j] * (-fac);
  }
  for (i=0; i<iele; i++)
  {
    dum = functe[i] * FB;
    Ked[i][j] = Ked[i][j] + dum;
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_stiff_ed */


/*!----------------------------------------------------------------------
\brief routine 'wge_stiff_ed' -> calculates stiffness
       matrix K_ed as sum over gaussian points

*----------------------------------------------------------------------*/
void wge_stiff_ee(DOUBLE  **Kee,      /* stiffnes equiv,strain - disp  */
                  DOUBLE  **bope,     /* B-operator for displacements  */
                  DOUBLE    crad,     /* influenceradius nonloc.strain */
                  DOUBLE   *functe,   /* Anatzfunc. for equiv.strain   */
                  DOUBLE    fac,      /* integration factor            */
                  INT       iele)      /* element equiv.strain DOF      */
{
#ifdef D_WALLGE

INT            i,j,k;
DOUBLE         dum1,dum2;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("wge_stiff_ee");
#endif
/*----------------------------------------------------------------------*/
for (j=0; j<iele; j++)
{
  for (i=0; i<iele; i++)
  {
    dum1 = 0.0;
    for (k=0; k<2; k++)
    {
      dum1 = dum1 + bope[k][i]*bope[k][j]*fac*crad;
    }
    dum2 = functe[i]*functe[j]*fac;
    Kee[i][j] = Kee[i][j] + dum1 + dum2;
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_stiff_ee */


/*!----------------------------------------------------------------------
\brief routine 'wge_fintd' -> calculates internal forces for displacement
       DOF's as sum over gauss points

*----------------------------------------------------------------------*/
void wge_fintd(DOUBLE    *stress,      /* stresses                      */
               DOUBLE     fac,         /* integration factor            */
               DOUBLE   **bopd,       /* B-operator for displacements  */
               INT        numdfd,     /* element-displacement DOF      */
               INT        neps,       /* number of strain components   */
               DOUBLE    *fintd)      /* int. force(displacementDOF)   */
{

INT            i,k;
DOUBLE         dum;

#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("wge_fintd");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<numdfd; i++)
{
  dum = 0.0;
  for (k=0; k<neps; k++)
  {
    dum = dum + bopd[k][i]*stress[k]*fac;
  }
  fintd[i] = fintd[i] + dum;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_fintd */


/*!----------------------------------------------------------------------
\brief routine 'wge_finte' -> calculates internal forces for equiv.strain
       DOF's as sum over gauss points

*----------------------------------------------------------------------*/
void wge_finte(DOUBLE     eps_vl,       /* local equiv. strain         */
               DOUBLE     eps_vnl,      /* nonlocal equiv. strain      */
               DOUBLE    *grad_eps_vnl, /* grad nonlocal equiv. strain */
               DOUBLE     crad,         /* material parameter          */
               DOUBLE     fac,          /* integration factor          */
               DOUBLE   **bope,         /* B-operator for equiv.strain */
               DOUBLE    *functe,       /* Ansatz-fun for equiv.strain */
               INT        iele,         /* eleDOF for equiv.strain     */
               DOUBLE    *finte)        /* int. force(equiv.strainDOF) */
{

INT            i,k;
DOUBLE         dum1,dum2;

#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("wge_finte");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<iele; i++)
{
  dum1 = 0.0;
  for (k=0; k<2; k++)
  {
    dum1 = dum1 + bope[k][i]*grad_eps_vnl[k]*fac*crad;
  }
  dum2 = functe[i]*fac*(eps_vnl - eps_vl);
  finte[i] = finte[i] + dum1 + dum2;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_finte */


/*!----------------------------------------------------------------------
\brief routine 'wge_permstiff' -> resort of stiffness parts Kdd,Kde,Ked,Kee
       into the element stiffness matrix "estif"

*----------------------------------------------------------------------*/
void wge_permstiff(DOUBLE    **Kdd,     /* 1.stiffness part            */
                   DOUBLE    **Kde,     /* 2.stiffness part            */
                   DOUBLE    **Ked,     /* 3.stiffness part            */
                   DOUBLE    **Kee,     /* 4.stiffness part            */
                   INT         iele,    /* node-number of equiv.strain */
                   INT         ield,    /* node-number of displacement */
                   INT         numdf,   /* total DOF                   */
                   DOUBLE    **estif)   /* "mixed" element stiffness   */
{

INT            i,j;
INT            nodei,nodej;
INT            nodestarti,nodestartj;
INT            dofi,dofj;

#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("wge_permstiff");
#endif
/*-------------------- (upper left part of estif) node 1-4 * node1-4 ---*/
for (nodei=0; nodei<iele; nodei++)
{
  for (nodej=0; nodej<iele; nodej++)
  {
    dofi       = 2 * nodei;
    dofj       = 2 * nodej;
    nodestarti = 3 * nodei;
    nodestartj = 3 * nodej;
    /*-----------------------------------------------------*/
    estif[nodestarti][nodestartj]     = Kdd[dofi][dofj];
    estif[nodestarti+1][nodestartj]   = Kdd[dofi+1][dofj];
    estif[nodestarti][nodestartj+1]   = Kdd[dofi][dofj+1];
    estif[nodestarti+1][nodestartj+1] = Kdd[dofi+1][dofj+1];
    /*-----------------------------------------------------*/
    estif[nodestarti][nodestartj+2]   = Kde[dofi][nodej];
    estif[nodestarti+1][nodestartj+2] = Kde[dofi+1][nodej];
    /*-----------------------------------------------------*/
    estif[nodestarti+2][nodestartj]   = Ked[nodei][dofj];
    estif[nodestarti+2][nodestartj+1] = Ked[nodei][dofj+1];
    /*-----------------------------------------------------*/
    estif[nodestarti+2][nodestartj+2] = Kee[nodei][nodej];
  }
}
/*-------------------------------------------- (rest part of estif) ---*/
if(ield>iele)
{
/*----------------- (upper right part of estif) node1-4 * node4-8/9 ---*/
  for (nodei=0; nodei<iele; nodei++)
  {
    for (j=12; j<numdf; j++)
    {
      dofi       = 2 * nodei;
      dofj       = j-4;
      nodestarti = 3 * nodei;
      /*-----------------------------------------------------*/
      estif[nodestarti][j]     = Kdd[dofi][dofj];
      estif[nodestarti+1][j]   = Kdd[dofi+1][dofj];
      /*-----------------------------------------------------*/
      estif[nodestarti+2][j]   = Ked[nodei][dofj];
    }
  }
/*------------------- (down left part of estif) node4-8/9 * node1-4 ---*/
  for (i=12; i<numdf; i++)
  {
    for (nodej=0; nodej<iele; nodej++)
    {
      dofi       = i-4;
      dofj       = 2 * nodej;
      nodestartj = 3 * nodej;
      /*-----------------------------------------------------*/
      estif[i][nodestartj]     = Kdd[dofi][dofj];
      estif[i][nodestartj+1]   = Kdd[dofi][dofj+1];
      /*-----------------------------------------------------*/
      estif[i][nodestartj+2]   = Kde[dofi][nodej];
    }
  }
/*---------------- (down right part of estif) node4-8/9 * node4-8/9 ---*/
  for (i=12; i<numdf; i++)
  {
    for (j=12; j<numdf; j++)
    {
      dofi       = i-4;
      dofj       = j-4;
      /*-----------------------------------------------------*/
      estif[i][j] = Kdd[dofi][dofj];
    }
  }
} /*-----endif-*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_permstiff */



/*!----------------------------------------------------------------------
\brief routine 'wge_permforce' -> resort of internal force parts fintd,finte
       into the element internal force "force"

*----------------------------------------------------------------------*/
void wge_permforce(DOUBLE     *fintd,   /* 1.part int. force           */
                   DOUBLE     *finte,   /* 2.part int. force           */
                   INT         iele,    /* num.of equiv.strain nodes   */
                   INT         ield,    /* num.of displacement nodes   */
                   INT         numdf,   /* total element DOF           */
                   DOUBLE     *force)   /* "mixed" element int. force  */
{

INT            i;
INT            nodei;
INT            nodestarti;
INT            dofi;

#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("wge_permforce");
#endif
/*----------------------------------- (upper part of force) node 1-4 ---*/
for (nodei=0; nodei<iele; nodei++)
{
    dofi       = 2 * nodei;
    nodestarti = 3 * nodei;
    /*---------------------------------------*/
    force[nodestarti]   = fintd[dofi];
    force[nodestarti+1] = fintd[dofi+1];
    /*---------------------------------------*/
    force[nodestarti+2] = finte[nodei];
}
/*--------------------------------------------- (rest part of force) ---*/
if(ield>iele)
{
  for (i=12; i<numdf; i++)
  {
    dofi = i - 4;
    force[i] = fintd[dofi];
  }
} /*-----endif-*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return;
} /* end of wge_permforce */





/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
