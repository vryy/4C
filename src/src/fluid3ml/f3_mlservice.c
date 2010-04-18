/*!----------------------------------------------------------------------
\file
\brief service routines for multilevel fluid3 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef FLUID3_ML
#include "../headers/standardtypes.h"
#include "fluid3ml_prototypes.h"
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

static FLUID_DYNAMIC *fdyn;

static INT PREDOF = 3;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief set all arrays for large-scale element for multi-level fluid3

<pre>                                                       gravem 07/03

In this routine, the element velocities, pressure and external loads
at time steps (n) and (n+1) are set.

</pre>
\param   *ele      ELEMENT	   (i)    actual element
\param  **ehist    DOUBLE	   (o)    ele history data
\param  **evel     DOUBLE	   (o)    ele vels at time step n+1
\param   *epre     DOUBLE	   (o)    ele pres at time step n+1
\param   *edeadn   DOUBLE          (o)    ele dead load at time step n
\param   *edead    DOUBLE          (o)    ele dead load at time step n+1
\param   *ipos                     (i)    node array positions
\param   *hasext   INT             (o)    flag for external loads
\return void

------------------------------------------------------------------------*/
void f3_lsset(ELEMENT	      *ele,
              DOUBLE	     **ehist,
	      DOUBLE	     **evel,
	      DOUBLE	      *epre,
	      DOUBLE	      *edeadn,
	      DOUBLE	      *edead,
              ARRAY_POSITION *ipos,
	      INT	      *hasext)
{
INT    i;           /* simply a counter                                 */
INT    actmat  ;    /* material number of the element                   */
DOUBLE dens;        /* density                                          */
NODE  *actnode;     /* actual node                                      */
GVOL  *actgvol;     /* actual g-volume                                  */

#ifdef DEBUG
dstrc_enter("f3_lsset");
#endif


/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[ipos.velnm][i]: solution at (n-1)               |
 |	 sol_increment[ipos.veln][i]:  solution at (n)                 |
 |	 sol_increment[ipos.velnp][i]: solution at (n+1)               |
 *---------------------------------------------------------------------*/
fdyn = alldyn[genprob.numff].fdyn;

for(i=0;i<ele->numnp;i++) /* loop nodes of large-scale element */
{
  actnode=ele->node[i];
/*------------------------------------ set element velocities at (n+1) */
  evel[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
  evel[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];
  evel[2][i]=actnode->sol_increment.a.da[ipos->velnp][2];
/*------------------------------------- set element pressures at (n+1) */
  epre[i]   =actnode->sol_increment.a.da[ipos->velnp][PREDOF];
/*---------------------------------------- set element history data ---*/
  ehist[0][i]=actnode->sol_increment.a.da[ipos->hist][0];
  ehist[1][i]=actnode->sol_increment.a.da[ipos->hist][1];
  ehist[2][i]=actnode->sol_increment.a.da[ipos->hist][2];
} /* end of loop over nodes of large-scale element */

/*------------------------------------------------ check for dead load */
actgvol = ele->g.gvol;
if (actgvol->neum!=NULL)
{
   actmat=ele->mat-1;
   dens = mat[actmat].m.fluid->density;
   for (i=0;i<3;i++)
   {
      if (actgvol->neum->neum_onoff.a.iv[i]==0)
      {
         edeadn[i] = ZERO;
	 edead[i]  = ZERO;
      }
      if (actgvol->neum->neum_type==neum_dead  &&
          actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadn[i] = actgvol->neum->neum_val.a.dv[i]*dens;
	 edead[i]  = actgvol->neum->neum_val.a.dv[i]*dens;
	 (*hasext)++;
      }
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_lsset */

/*!---------------------------------------------------------------------
\brief set all arrays for (sub-)submesh element for multi-level fluid3

<pre>                                                       gravem 07/03

In this routine, the (sub-)submesh element coordinates, topology array
and location matrix are set.

</pre>
\param   *smesh    FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *smlme    INT   	   (o)    (sub-)submesh location matrix
\param   *smitope  INT   	   (o)    (sub-)submesh topology array
\param   *smxyze   DOUBLE   	   (o)    (sub-)submesh coordinates
\param   *smxyzep  DOUBLE   	   (o)    (s)sm parent domain coordinates
\param    iele     INT             (i)    actual submesh element number
\param    flag     INT             (i)    flag: submesh or sub-submesh?
\return void

------------------------------------------------------------------------*/
void f3_smset(FLUID_ML_SMESH  *smesh,
	      ELEMENT	      *ele,
              INT 	      *smlme,
	      INT	      *smitope,
	      DOUBLE	     **smxyze,
	      DOUBLE	     **smxyzep,
	      INT	       iele,
	      INT	       flag)
{
INT i;       /* simply a counter                                       */

#ifdef DEBUG
dstrc_enter("f3_smset");
#endif

for (i=0; i<smesh->numen; i++)/* loop nodes of (sub-)submesh element */
{
  smitope[i]    = smesh->ien.a.ia[iele][i];
  smlme[i]      = smesh->id.a.iv[smitope[i]];
  smxyzep[0][i] = smesh->xyzpd.a.da[0][smitope[i]];
  smxyzep[1][i] = smesh->xyzpd.a.da[1][smitope[i]];
  smxyzep[2][i] = smesh->xyzpd.a.da[2][smitope[i]];
  switch (flag)
  {
/*-------------------------------------------------------- sub-submesh */
  case 1:
    smxyze[0][i] = ele->e.f3->xyzssm.a.da[0][smitope[i]];
    smxyze[1][i] = ele->e.f3->xyzssm.a.da[1][smitope[i]];
    smxyze[2][i] = ele->e.f3->xyzssm.a.da[2][smitope[i]];
  break;
/*------------------------------------------------------------ submesh */
  default:
    smxyze[0][i] = ele->e.f3->xyzsm.a.da[0][smitope[i]];
    smxyze[1][i] = ele->e.f3->xyzsm.a.da[1][smitope[i]];
    smxyze[2][i] = ele->e.f3->xyzsm.a.da[2][smitope[i]];
  break;
  }
}
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_smset */

/*!---------------------------------------------------------------------
\brief set all bubble functions for submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh element velocity, pressure and rhs bubble
functions are set.

</pre>
\param   *mlvar    FLUID_DYN_ML    (i)
\param   *submesh  FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *smlme    INT   	   (i)    submesh location matrix
\param  **evbub    DOUBLE   	   (o)    submesh element velocity bubble
\param  **epbub    DOUBLE   	   (o)    submesh element pressure bubble
\param  **efbub    DOUBLE   	   (o)    submesh element rhs bubble
\param    flag     INT             (i)    flag: time step (n) or (n+1)?
\return void

------------------------------------------------------------------------*/
void f3_bubset(FLUID_DYN_ML    *mlvar,
               FLUID_ML_SMESH  *submesh,
	       ELEMENT	       *ele,
               INT	       *smlme,
	       DOUBLE         **evbub,
	       DOUBLE         **epbub,
	       DOUBLE         **efbub,
	       INT	        flag)
{
INT i,j,ipbub,ifbub;/* simply some counters                            */
INT km;             /* value in location matrix                        */

#ifdef DEBUG
dstrc_enter("f3_bubset");
#endif

/*---- calculate velocity bubble functions at nodes of submesh element */
for (i=0;i<mlvar->nvbub;i++)
{
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) evbub[j][i]=ZERO;
    else if (km>=submesh->numeq) continue;
    else
    {
/*---------------------------------------------------- time step (n+1) */
      if (flag==0) evbub[j][i]= ele->e.f3->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */
      else         evbub[j][i]= ele->e.f3->solsmn.a.da[km][i];
    }
  }
}

/*---- calculate pressure bubble functions at nodes of submesh element */
for (i=mlvar->nvbub;i<mlvar->nelbub-3;i++)
{
  ipbub=i-mlvar->nvbub;
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) epbub[j][ipbub]=ZERO;
    else if (km>submesh->numeq) continue;
    else
    {
/*---------------------------------------------------- time step (n+1) */
      if (flag==0) epbub[j][ipbub]= ele->e.f3->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */
      else         epbub[j][ipbub]= ele->e.f3->solsmn.a.da[km][i];
    }
  }
}

/*--------- calculate rhs bubble functions at nodes of submesh element */
for (i=mlvar->nelbub-3;i<mlvar->nelbub;i++)
{
  ifbub=i-(mlvar->nvbub+mlvar->npbub);
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) efbub[j][ifbub]=ZERO;
    else if (km>submesh->numeq) continue;
    else
    {
/*---------------------------------------------------- time step (n+1) */
      if (flag==0) efbub[j][ifbub]= ele->e.f3->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */
      else         efbub[j][ifbub]= ele->e.f3->solsmn.a.da[km][i];
    }
  }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_bubset */

/*!---------------------------------------------------------------------
\brief set all arrays for integration of sub-submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the sub-submesh element arrays for the elementwise
integration of a normalized bubble function are set.

</pre>
\param   *ssmesh   FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *sslme    INT   	   (o)    sub-submesh location matrix
\param   *ssitope  INT   	   (o)    sub-submesh topology array
\param   *ssxyze   DOUBLE   	   (o)    sub-submesh coordinates
\param   *ebub     DOUBLE   	   (o)    sub-submesh element bubble
\param    iele     INT             (i)    actual sub-submesh ele. number
\return void

------------------------------------------------------------------------*/
void f3_ssset(FLUID_ML_SMESH  *ssmesh,
	      ELEMENT	      *ele,
              INT 	      *sslme,
	      INT	      *ssitope,
	      DOUBLE	     **ssxyze,
	      DOUBLE          *ebub,
	      INT	       iele)
{
INT i;              /* simply a counter                                */
INT km;             /* value in location matrix                        */

#ifdef DEBUG
dstrc_enter("f3_ssset");
#endif

for (i=0; i<ssmesh->numen; i++)/* loop nodes of sub-submesh element */
{
/*----------------------------------------- sub-submesh element arrays */
  ssitope[i]   = ssmesh->ien.a.ia[iele][i];
  sslme[i]     = ssmesh->id.a.iv[ssitope[i]];
  ssxyze[0][i] = ele->e.f3->xyzssm.a.da[0][ssitope[i]];
  ssxyze[1][i] = ele->e.f3->xyzssm.a.da[1][ssitope[i]];
  ssxyze[2][i] = ele->e.f3->xyzssm.a.da[2][ssitope[i]];

/*- calculate normalized bubble funct. at nodes of sub-submesh element */
  km=sslme[i];
  if (km==-1) ebub[i]=ZERO;
  else  ebub[i]= ssmesh->rhs.a.dv[km];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_ssset */

/*!---------------------------------------------------------------------
\brief copy submesh solution to element array for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh solution is copied to the respective
element array.

</pre>
\param   *smrhs    DOUBLE          (i)    submesh solution array
\param   *ele      ELEMENT	   (o)    actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void

------------------------------------------------------------------------*/
void f3_smcopy(DOUBLE  **smrhs,
               ELEMENT  *ele,
               INT       numeq,
               INT       numrhs)
{
INT  i,j;

#ifdef DEBUG
dstrc_enter("f3_smcopy");
#endif

for (i=0;i<numrhs;i++)
{
  for (j=0;j<numeq;j++)
  {
    ele->e.f3->solsm.a.da[j][i] = smrhs[j][i];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_smcopy */

/*!---------------------------------------------------------------------
\brief copy submesh element solution at (n+1) to (n) for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh element solution array at time step (n+1)
is copied to the respective array at time step (n).

</pre>
\param   *ele      ELEMENT	   (i/o)  actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void

------------------------------------------------------------------------*/
void f3_smcopy2(ELEMENT  *ele,
                INT       numeq,
                INT       numrhs)
{
INT  i,j;

#ifdef DEBUG
dstrc_enter("f3_smcopy2");
#endif

for (i=0;i<numrhs;i++)
{
  for (j=0;j<numeq;j++)
  {
    ele->e.f3->solsmn.a.da[j][i] = ele->e.f3->solsm.a.da[j][i];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_smcopy2 */

/*!---------------------------------------------------------------------
\brief routine to calculate small-scale pressure at int. p. for fluid3

<pre>                                                       gravem 07/03

</pre>
\param   *smpreint  DOUBLE     (o)   small-scale pressure at int. point
\param  **pbubint   DOUBLE     (i)   pressure bubble functions at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    INT        (i)   number of large-scale ele. nodes
\return void

------------------------------------------------------------------------*/
void f3_smprei(DOUBLE  *smpreint,
               DOUBLE **pbubint,
	       DOUBLE  *epre,
	       INT      iel)
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_smprei");
#endif

for (i=0;i<3;i++) /* loop spatial directions i */
{
   smpreint[i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpreint[i] += pbubint[i][j]*epre[j];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_smprei */

/*!---------------------------------------------------------------------
\brief routine to calculate s-s pressure deriv. at int. p. for fluid3

<pre>                                                       gravem 07/03

In this routine, the derivatives of the small-scale pressure w.r.t x/y/z
are calculated.

</pre>
\param   *smpderxy  DOUBLE     (o)   s-s pressure deriv. at int. point
\param  **pbubderxy DOUBLE     (i)   pre. bubble fun. deriv. at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    INT        (i)   number of large-scale ele. nodes
\return void

------------------------------------------------------------------------*/
void f3_smpder(DOUBLE  **smpderxy,
               DOUBLE ***pbubderxy,
	       DOUBLE   *epre,
	       INT       iel)
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_smpder");
#endif

for (i=0;i<3;i++) /* loop spatial directions i */
{
   smpderxy[0][i]=ZERO;
   smpderxy[1][i]=ZERO;
   smpderxy[2][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpderxy[0][i] += pbubderxy[i][0][j]*epre[j];
      smpderxy[1][i] += pbubderxy[i][1][j]*epre[j];
      smpderxy[2][i] += pbubderxy[i][2][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_smpder */

/*!---------------------------------------------------------------------
\brief routine to calc. 2nd s-s pressure deriv. at int. p. for fluid3

<pre>                                                       gravem 07/03

In this routine, the 2nd derivatives of the small-scale pressure w.r.t
x/y/z are calculated.

</pre>
\param  **smpderxy2  DOUBLE     (o)   2nd s-s pre. deriv. at int. point
\param ***pbubderxy2 DOUBLE     (i)   2nd pre. bub. fun. deriv. at i. p.
\param   *epre       DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	     INT        (i)   number of large-scale ele. nodes
\return void

------------------------------------------------------------------------*/
void f3_smpder2(DOUBLE  **smpderxy2,
                DOUBLE ***pbubderxy2,
	        DOUBLE   *epre,
	        INT       iel)
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f3_smpder2");
#endif

for (i=0;i<6;i++) /* loop spatial directions i */
{
   smpderxy2[0][i]=ZERO;
   smpderxy2[1][i]=ZERO;
   smpderxy2[2][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpderxy2[0][i] += pbubderxy2[i][0][j]*epre[j];
      smpderxy2[1][i] += pbubderxy2[i][1][j]*epre[j];
      smpderxy2[2][i] += pbubderxy2[i][2][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_smpder2 */

void f3_mlpermestif(DOUBLE         **estif,
		    DOUBLE	 **emass,
		    DOUBLE	 **tmp,
		    INT		   iel)
{
INT i,j,icol,irow;          /* simply some counters  	        	*/
INT nvdof;                  /* number of vel dofs 			*/
INT npdof;                  /* number of pre dofs 			*/
INT totdof;                 /* total number of dofs			*/
DOUBLE thsl,thpl;	    /* factor for LHS (THETA*DT)		*/

#ifdef DEBUG
dstrc_enter("f3_mlpermestif");
#endif

fdyn = alldyn[genprob.numff].fdyn;

nvdof  = NUM_F3_VELDOF*iel;
npdof  = iel;
totdof = (NUM_F3_VELDOF+1)*iel;
thsl   = fdyn->thsl;
thpl   = fdyn->thpl;

/*--------------------------------------------- copy estif to tmp-array *
                and multitply stiffness matrix with respective THETA*DT */
for (i=0;i<totdof;i++)
{
/*--------------------------------------------------------- Kvv and Kpv */
   for (j=0;j<nvdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl;
   } /* end of loop over j */
/*--------------------------------------------------------- Kvp and Kpp */
   for (j=nvdof;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thpl;
   } /* end of loop over j */
} /* end of loop over i */
/*------------------------------- add mass matrix for instationary case */
if (fdyn->nis==0)
{
   for (i=0;i<totdof;i++)
   {
      for (j=0;j<nvdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (fdyn->nis==0) */

/*--------------------------------------------------------- compute Kvv */
irow = 0;
for (i=0;i<nvdof;i+=3)
{
   icol = 0;
   for (j=0;j<nvdof;j+=3)
   {
      estif[irow][icol]     = tmp[i][j];
      estif[irow+1][icol]   = tmp[i+1][j];
      estif[irow+2][icol]   = tmp[i+2][j];
      estif[irow][icol+1]   = tmp[i][j+1];
      estif[irow+1][icol+1] = tmp[i+1][j+1];
      estif[irow+2][icol+1] = tmp[i+2][j+1];
      estif[irow][icol+2]   = tmp[i][j+2];
      estif[irow+1][icol+2] = tmp[i+1][j+2];
      estif[irow+2][icol+2] = tmp[i+2][j+2];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*--------------------------------------------------------- compute Kvp */
irow = 0;
for (i=0;i<nvdof;i+=3)
{
   icol = 3;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow+1][icol] = tmp[i+1][j];
      estif[irow+2][icol] = tmp[i+2][j];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*--------------------------------------------------------- compute Kpv */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   icol = 0;
   for (j=0;j<nvdof;j+=3)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow][icol+1] = tmp[i][j+1];
      estif[irow][icol+2] = tmp[i][j+2];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*--------------------------------------------------------- compute Kpp */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   icol = 3;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol] = tmp[i][j];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_mlpermestif */

void f3_mljaco(DOUBLE     *funct,
               DOUBLE    **deriv,
               DOUBLE    **xjm,
               DOUBLE     *det,
               ELEMENT    *ele,
               INT         iel)

{
INT i,j,l;
DOUBLE dum;

#ifdef DEBUG
dstrc_enter("f3_mljaco");
#endif

/*-------------------------------- determine jacobian at point r,s,t ---*/
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      dum=0.0;
      for (l=0; l<iel; l++)
      {
         dum += deriv[i][l]*ele->node[l]->x[j];
      }
      xjm[i][j]=dum;
   } /* end loop j */
} /* end loop i */
/*------------------------------------------ determinant of jacobian ---*/
*det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
       xjm[0][1]*xjm[1][2]*xjm[2][0]+
       xjm[0][2]*xjm[1][0]*xjm[2][1]-
       xjm[0][2]*xjm[1][1]*xjm[2][0]-
       xjm[0][0]*xjm[1][2]*xjm[2][1]-
       xjm[0][1]*xjm[1][0]*xjm[2][2];

if(*det<ZERO)
{
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   dserror("NEGATIVE JACOBIAN DETERMINANT\n");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_mljaco */

void f3_mljaco3(DOUBLE    **xyze,
                DOUBLE	 *funct,
                DOUBLE	**deriv,
                DOUBLE	**xjm,
                DOUBLE	 *det,
                INT	  iel,
                ELEMENT	 *ele)
{
INT i,j,l;        /* just some counters                                 */
DOUBLE dum;       /* dummy variable                                     */

#ifdef DEBUG
dstrc_enter("f3_mljaco3");
#endif

/*-------------------------------- determine jacobian at point r,s,t ---*/
for (i=0; i<3; i++)
{
   for (j=0; j<3; j++)
   {
      dum=0.0;
      for (l=0; l<iel; l++)
      {
         dum += deriv[i][l]*xyze[j][l];
      }
      xjm[i][j]=dum;
   } /* end loop j */
} /* end loop i */
/*------------------------------------------ determinant of jacobian ---*/
*det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
       xjm[0][1]*xjm[1][2]*xjm[2][0]+
       xjm[0][2]*xjm[1][0]*xjm[2][1]-
       xjm[0][2]*xjm[1][1]*xjm[2][0]-
       xjm[0][0]*xjm[1][2]*xjm[2][1]-
       xjm[0][1]*xjm[1][0]*xjm[2][2];


if(*det<ZERO)
{
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   dserror("NEGATIVE JACOBIAN DETERMINANT\n");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_mljaco3 */

void f3_mlgcoor2(DOUBLE     *funct,
                 DOUBLE    **xyze,
	         INT         iel,
	         DOUBLE     *gcoor)
{
INT i;

#ifdef DEBUG
dstrc_enter("f3_mlgcoor2");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;
gcoor[2]=ZERO;

for(i=0;i<iel;i++) /* loop all nodes of the element */
{
   gcoor[0] += funct[i] * xyze[0][i];
   gcoor[1] += funct[i] * xyze[1][i];
   gcoor[2] += funct[i] * xyze[2][i];
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_mlgcoor2 */

void f3_mlgder2(ELEMENT     *ele,
	        DOUBLE	 **xjm,
	        DOUBLE	 **xder2,
	        DOUBLE	 **derxy,
	        DOUBLE	 **derxy2,
                DOUBLE	 **deriv2,
	        INT	   iel)
{
INT i,j;
DOUBLE r0,r1,r2,r3,r4,r5;
INT irc=0;
DOUBLE bb[36],ba[36];

#ifdef DEBUG
dstrc_enter("f3_mlgder2");
#endif

/*--------------------------- calculate elements of jacobian_bar matrix */
bb[0] = xjm[0][0]*xjm[0][0];
bb[1] = xjm[1][0]*xjm[1][0];
bb[2] = xjm[2][0]*xjm[2][0];
bb[3] = xjm[0][0]*xjm[1][0];
bb[4] = xjm[0][0]*xjm[2][0];
bb[5] = xjm[1][0]*xjm[2][0];

bb[6] = xjm[0][1]*xjm[0][1];
bb[7] = xjm[1][1]*xjm[1][1];
bb[8] = xjm[2][1]*xjm[2][1];
bb[9] = xjm[0][1]*xjm[1][1];
bb[10] = xjm[0][1]*xjm[2][1];
bb[11] = xjm[1][1]*xjm[2][1];

bb[12] = xjm[0][2]*xjm[0][2];
bb[13] = xjm[1][2]*xjm[1][2];
bb[14] = xjm[2][2]*xjm[2][2];
bb[15] = xjm[0][2]*xjm[1][2];
bb[16] = xjm[0][2]*xjm[2][2];
bb[17] = xjm[1][2]*xjm[2][2];

bb[18] = TWO*xjm[0][0]*xjm[0][1];
bb[19] = TWO*xjm[1][0]*xjm[1][1];
bb[20] = TWO*xjm[2][0]*xjm[2][1];
bb[21] = xjm[0][0]*xjm[1][1]+xjm[1][0]*xjm[0][1];
bb[22] = xjm[0][0]*xjm[2][1]+xjm[2][0]*xjm[0][1];
bb[23] = xjm[1][0]*xjm[2][1]+xjm[2][0]*xjm[1][1];

bb[24] = TWO*xjm[0][0]*xjm[0][2];
bb[25] = TWO*xjm[1][0]*xjm[1][2];
bb[26] = TWO*xjm[2][0]*xjm[2][2];
bb[27] = xjm[0][0]*xjm[1][2]+xjm[1][0]*xjm[0][2];
bb[28] = xjm[0][0]*xjm[2][2]+xjm[2][0]*xjm[0][2];
bb[29] = xjm[1][0]*xjm[2][2]+xjm[2][0]*xjm[1][2];

bb[30] = TWO*xjm[0][1]*xjm[0][2];
bb[31] = TWO*xjm[1][1]*xjm[1][2];
bb[32] = TWO*xjm[2][1]*xjm[2][2];
bb[33] = xjm[0][1]*xjm[1][2]+xjm[1][1]*xjm[0][2];
bb[34] = xjm[0][1]*xjm[2][2]+xjm[2][1]*xjm[0][2];
bb[35] = xjm[1][1]*xjm[2][2]+xjm[2][1]*xjm[1][2];

/*-------------------------------------- inverse of jacobian_bar matrix */
c1inv6(&bb,&ba,&irc);
if (irc!=0) dserror("matrix inversion 6x6 failed!");

/*----------------------------------------------------------- initialise*/
for (i=0;i<3;i++)
{
   for (j=0;j<6;j++) xder2[j][i]=ZERO;
}
for (i=0;i<iel;i++)
{
   for (j=0;j<6;j++) derxy2[j][i]=ZERO;
}

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++)
{
   xder2[0][0] += deriv2[0][i] * ele->node[i]->x[0];
   xder2[1][0] += deriv2[1][i] * ele->node[i]->x[0];
   xder2[2][0] += deriv2[2][i] * ele->node[i]->x[0];
   xder2[3][0] += deriv2[3][i] * ele->node[i]->x[0];
   xder2[4][0] += deriv2[4][i] * ele->node[i]->x[0];
   xder2[5][0] += deriv2[5][i] * ele->node[i]->x[0];

   xder2[0][1] += deriv2[0][i] * ele->node[i]->x[1];
   xder2[1][1] += deriv2[1][i] * ele->node[i]->x[1];
   xder2[2][1] += deriv2[2][i] * ele->node[i]->x[1];
   xder2[3][1] += deriv2[3][i] * ele->node[i]->x[1];
   xder2[4][1] += deriv2[4][i] * ele->node[i]->x[1];
   xder2[5][1] += deriv2[5][i] * ele->node[i]->x[1];

   xder2[0][2] += deriv2[0][i] * ele->node[i]->x[2];
   xder2[1][2] += deriv2[1][i] * ele->node[i]->x[2];
   xder2[2][2] += deriv2[2][i] * ele->node[i]->x[2];
   xder2[3][2] += deriv2[3][i] * ele->node[i]->x[2];
   xder2[4][2] += deriv2[4][i] * ele->node[i]->x[2];
   xder2[5][2] += deriv2[5][i] * ele->node[i]->x[2];
} /* end of loop over i */

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++)
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i] \
                     - xder2[0][2]*derxy[2][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i] \
                     - xder2[1][2]*derxy[2][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i] \
                     - xder2[2][2]*derxy[2][i];
   r3 = deriv2[3][i] - xder2[3][0]*derxy[0][i] - xder2[3][1]*derxy[1][i] \
                     - xder2[3][2]*derxy[2][i];
   r4 = deriv2[4][i] - xder2[4][0]*derxy[0][i] - xder2[4][1]*derxy[1][i] \
                     - xder2[4][2]*derxy[2][i];
   r5 = deriv2[5][i] - xder2[5][0]*derxy[0][i] - xder2[5][1]*derxy[1][i] \
                     - xder2[5][2]*derxy[2][i];

   derxy2[0][i] += ba[0]*r0 + ba[6]*r1 + ba[12]*r2 \
                +  ba[18]*r3 + ba[24]*r4 + ba[30]*r5;
   derxy2[1][i] += ba[1]*r0 + ba[7]*r1 + ba[13]*r2 \
                +  ba[19]*r3 + ba[25]*r4 + ba[31]*r5;
   derxy2[2][i] += ba[2]*r0 + ba[8]*r1 + ba[14]*r2 \
                +  ba[20]*r3 + ba[26]*r4 + ba[32]*r5;
   derxy2[3][i] += ba[3]*r0 + ba[9]*r1 + ba[15]*r2 \
                +  ba[21]*r3 + ba[27]*r4 + ba[33]*r5;
   derxy2[4][i] += ba[4]*r0 + ba[10]*r1 + ba[16]*r2 \
                +  ba[22]*r3 + ba[28]*r4 + ba[34]*r5;
   derxy2[5][i] += ba[5]*r0 + ba[11]*r1 + ba[17]*r2 \
                +  ba[23]*r3 + ba[29]*r4 + ba[35]*r5;
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_mlgder2 */

void f3_mlcogder2(DOUBLE     **xyze,
	          DOUBLE     **xjm,
	          DOUBLE     **xder2,
	          DOUBLE     **derxy,
	          DOUBLE     **derxy2,
                  DOUBLE     **deriv2,
	          INT          iel)
{
INT i,j;
DOUBLE r0,r1,r2,r3,r4,r5;
INT irc=0;
DOUBLE bb[36],ba[36];

#ifdef DEBUG
dstrc_enter("f3_mlcogder2");
#endif

/*--------------------------- calculate elements of jacobian_bar matrix */
bb[0] = xjm[0][0]*xjm[0][0];
bb[1] = xjm[1][0]*xjm[1][0];
bb[2] = xjm[2][0]*xjm[2][0];
bb[3] = xjm[0][0]*xjm[1][0];
bb[4] = xjm[0][0]*xjm[2][0];
bb[5] = xjm[1][0]*xjm[2][0];

bb[6] = xjm[0][1]*xjm[0][1];
bb[7] = xjm[1][1]*xjm[1][1];
bb[8] = xjm[2][1]*xjm[2][1];
bb[9] = xjm[0][1]*xjm[1][1];
bb[10] = xjm[0][1]*xjm[2][1];
bb[11] = xjm[1][1]*xjm[2][1];

bb[12] = xjm[0][2]*xjm[0][2];
bb[13] = xjm[1][2]*xjm[1][2];
bb[14] = xjm[2][2]*xjm[2][2];
bb[15] = xjm[0][2]*xjm[1][2];
bb[16] = xjm[0][2]*xjm[2][2];
bb[17] = xjm[1][2]*xjm[2][2];

bb[18] = TWO*xjm[0][0]*xjm[0][1];
bb[19] = TWO*xjm[1][0]*xjm[1][1];
bb[20] = TWO*xjm[2][0]*xjm[2][1];
bb[21] = xjm[0][0]*xjm[1][1]+xjm[1][0]*xjm[0][1];
bb[22] = xjm[0][0]*xjm[2][1]+xjm[2][0]*xjm[0][1];
bb[23] = xjm[1][0]*xjm[2][1]+xjm[2][0]*xjm[1][1];

bb[24] = TWO*xjm[0][0]*xjm[0][2];
bb[25] = TWO*xjm[1][0]*xjm[1][2];
bb[26] = TWO*xjm[2][0]*xjm[2][2];
bb[27] = xjm[0][0]*xjm[1][2]+xjm[1][0]*xjm[0][2];
bb[28] = xjm[0][0]*xjm[2][2]+xjm[2][0]*xjm[0][2];
bb[29] = xjm[1][0]*xjm[2][2]+xjm[2][0]*xjm[1][2];

bb[30] = TWO*xjm[0][1]*xjm[0][2];
bb[31] = TWO*xjm[1][1]*xjm[1][2];
bb[32] = TWO*xjm[2][1]*xjm[2][2];
bb[33] = xjm[0][1]*xjm[1][2]+xjm[1][1]*xjm[0][2];
bb[34] = xjm[0][1]*xjm[2][2]+xjm[2][1]*xjm[0][2];
bb[35] = xjm[1][1]*xjm[2][2]+xjm[2][1]*xjm[1][2];

/*-------------------------------------- inverse of jacobian_bar matrix */
c1inv6(&bb,&ba,&irc);
if (irc!=0) dserror("matrix inversion 6x6 failed!");

/*----------------------------------------------------------- initialise*/
for (i=0;i<3;i++)
{
   for (j=0;j<6;j++) xder2[j][i]=ZERO;
}
for (i=0;i<iel;i++)
{
   for (j=0;j<6;j++) derxy2[j][i]=ZERO;
}

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++)
{
   xder2[0][0] += deriv2[0][i] * xyze[0][i];
   xder2[1][0] += deriv2[1][i] * xyze[0][i];
   xder2[2][0] += deriv2[2][i] * xyze[0][i];
   xder2[3][0] += deriv2[3][i] * xyze[0][i];
   xder2[4][0] += deriv2[4][i] * xyze[0][i];
   xder2[5][0] += deriv2[5][i] * xyze[0][i];

   xder2[0][1] += deriv2[0][i] * xyze[1][i];
   xder2[1][1] += deriv2[1][i] * xyze[1][i];
   xder2[2][1] += deriv2[2][i] * xyze[1][i];
   xder2[3][1] += deriv2[3][i] * xyze[1][i];
   xder2[4][1] += deriv2[4][i] * xyze[1][i];
   xder2[5][1] += deriv2[5][i] * xyze[1][i];

   xder2[0][2] += deriv2[0][i] * xyze[2][i];
   xder2[1][2] += deriv2[1][i] * xyze[2][i];
   xder2[2][2] += deriv2[2][i] * xyze[2][i];
   xder2[3][2] += deriv2[3][i] * xyze[2][i];
   xder2[4][2] += deriv2[4][i] * xyze[2][i];
   xder2[5][2] += deriv2[5][i] * xyze[2][i];
} /* end of loop over i */

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++)
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i] \
                     - xder2[0][2]*derxy[2][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i] \
                     - xder2[1][2]*derxy[2][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i] \
                     - xder2[2][2]*derxy[2][i];
   r3 = deriv2[3][i] - xder2[3][0]*derxy[0][i] - xder2[3][1]*derxy[1][i] \
                     - xder2[3][2]*derxy[2][i];
   r4 = deriv2[4][i] - xder2[4][0]*derxy[0][i] - xder2[4][1]*derxy[1][i] \
                     - xder2[4][2]*derxy[2][i];
   r5 = deriv2[5][i] - xder2[5][0]*derxy[0][i] - xder2[5][1]*derxy[1][i] \
                     - xder2[5][2]*derxy[2][i];

   derxy2[0][i] += ba[0]*r0 + ba[6]*r1 + ba[12]*r2 \
                +  ba[18]*r3 + ba[24]*r4 + ba[30]*r5;
   derxy2[1][i] += ba[1]*r0 + ba[7]*r1 + ba[13]*r2 \
                +  ba[19]*r3 + ba[25]*r4 + ba[31]*r5;
   derxy2[2][i] += ba[2]*r0 + ba[8]*r1 + ba[14]*r2 \
                +  ba[20]*r3 + ba[26]*r4 + ba[32]*r5;
   derxy2[3][i] += ba[3]*r0 + ba[9]*r1 + ba[15]*r2 \
                +  ba[21]*r3 + ba[27]*r4 + ba[33]*r5;
   derxy2[4][i] += ba[4]*r0 + ba[10]*r1 + ba[16]*r2 \
                +  ba[22]*r3 + ba[28]*r4 + ba[34]*r5;
   derxy2[5][i] += ba[5]*r0 + ba[11]*r1 + ba[17]*r2 \
                +  ba[23]*r3 + ba[29]*r4 + ba[35]*r5;
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_mlcogder2 */

#endif
#endif
