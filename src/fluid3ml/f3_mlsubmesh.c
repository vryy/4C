/*!----------------------------------------------------------------------
\file
\brief submesh creation for multilevel fluid3 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID3_ML
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid3ml_prototypes.h"
#include "../fluid3/fluid3.h"

static int nxsm,nysm,nzsm,nxss,nyss,nzss;

/*!---------------------------------------------------------------------
\brief creation of (sub-)submesh on parent domain for fluid3

<pre>                                                       gravem 07/03

In this routine, the creation of submesh or sub-submesh on the parent
domain is performed, i.e. nodal coordinates on the parent domain,
id-array and ien-array are established.

</pre>
\param  *smesh      FLUID_ML_SMESH   (i/o)
\param   xele       INT              (i)    number of elements in x-dir.
\param   yele       INT              (i)    number of elements in y-dir.
\param   zele       INT              (i)    number of elements in z-dir.
\param   order      INT              (i)    polyn. interpolation order
\param   flag       INT              (i)    flag: submesh or sub-submesh?
\return void

------------------------------------------------------------------------*/
void f3_pdsubmesh(FLUID_ML_SMESH *smesh,
                  INT             xele,
                  INT             yele,
                  INT             zele,
                  INT             order,
		  INT             flag)
{
DOUBLE hpdx,hpdy,hpdz;/* element length in coord. dir. on parent domain */
INT    numeq;         /* number of equations                            */
INT    numnp;         /* number of nodal points                         */
INT    ix,iy,iz;      /* counters in coordinate directions              */
INT    nx,ny,nz;      /* number of nodes in coordinate directions       */
INT    iel;           /* element counter                                */
INT    nnglo;         /* counter for global node number                 */
INT    nnhor;         /* counter for number of nod. in a horizontal row */
INT    nrhor;         /* counter for number of horizontal rows          */
INT    nrpla;         /* counter for number of rows in a plane          */
INT    nnpla;         /* number of nodes in a plane                     */
INT    npact;         /* actual plane                                   */

#ifdef DEBUG
dstrc_enter("f3_pdsubmesh");
#endif

/*----------------------------------------------------------------------
  calculate coordinates of nodes on parent domain and evaluate id-array
------------------------------------------------------------------------*/
/*------------------------------------ divide large-scale element edges */
hpdx  = ONE/(xele*order);
hpdy  = ONE/(yele*order);
hpdz  = ONE/(zele*order);
switch (flag)
{
/*--------------------------------------------------------- sub-submesh */
case 1:
  nxss = order*xele+1;
  nyss = order*yele+1;
  nzss = order*zele+1;
  nx = nxss;
  ny = nyss;
  nz = nzss;
break;
/*------------------------------------------------------------- submesh */
default:
  nxsm = order*xele+1;
  nysm = order*yele+1;
  nzsm = order*zele+1;
  nx = nxsm;
  ny = nysm;
  nz = nzsm;
break;
}
/*--------------------------------------- set node and equation counter */
numeq = 0;
numnp = 0;

for (iz=0; iz<nz; iz++)
{
  for (iy=0; iy<ny; iy++)
  {
    for (ix=0; ix<nx; ix++)
    {
/*-------------------------------- compute coordinates on parent domain */
      smesh->xyzpd.a.da[0][numnp] = -ONE+TWO*ix*hpdx;
      smesh->xyzpd.a.da[1][numnp] = -ONE+TWO*iy*hpdy;
      smesh->xyzpd.a.da[2][numnp] = -ONE+TWO*iz*hpdz;
/*--------------------------------------------------- evaluate id-array */
      if (ix==0 || ix==nx-1 || iy==0 || iy==ny-1 || iz==0 || iz==nz-1)
/*------------------------------------------------------- boundary node */
        smesh->id.a.iv[numnp] = -1;
      else
      {
/*--------------------------------------------------- non-boundary node */
        smesh->id.a.iv[numnp] = numeq;
        numeq++;
      }
      numnp++;
    }
  }
}

/*----------------------------------------------------------------------
   evaluate ien-array (only linear elements so far)
----------------------------------------------------------------------- */
/*----------------------------------------------------- linear elements */
if (order==1)
{
/*------------------------------------------ set row and plane counters */
  nnhor = 1;
  nrhor = 1;
  npact = 1;
  nrpla = 1;
  nnpla = (xele+1)*(yele+1);
  for (iel=0; iel<smesh->numele; iel++)
  {
    if (nnhor==nx)
    {
/*-------------------------------------------- start new horizontal row */
      nnhor = 1;
      nrhor++;
      nrpla++;
    }
    if (nrpla==ny)
    {
/*----------------------------------------------------- start new plane */
      nrpla = 1;
      npact++;
    }
/*---------------------------------------- evaluate global node counter */
    nnglo = iel+nrhor-1+nx*(npact-1);
/*----------------------------------------- assign element node numbers */
    smesh->ien.a.ia[iel][0] = nnglo;
    smesh->ien.a.ia[iel][1] = nnglo+1;
    smesh->ien.a.ia[iel][2] = nnglo+xele+2;
    smesh->ien.a.ia[iel][3] = nnglo+xele+1;
    smesh->ien.a.ia[iel][4] = nnpla+nnglo;
    smesh->ien.a.ia[iel][5] = nnpla+nnglo+1;
    smesh->ien.a.ia[iel][6] = nnpla+nnglo+yele+2;
    smesh->ien.a.ia[iel][7] = nnpla+nnglo+yele+1;
    nnhor++;
  }
}
else dserror("no other than linear submeshes in 3D!\n");
/*--------------------------------------------------------------------- */
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_pdsubmesh */

/*!---------------------------------------------------------------------
\brief creation of (sub-)submesh on individual element for fluid3

<pre>                                                       gravem 07/03

In this routine, the creation of submesh or sub-submesh on an individual
element is performed, i.e. the nodal coordinates on this particular
element are established.

</pre>
\param  *ele        ELEMENT          (i/o)
\param  *smesh      FLUID_ML_SMESH   (i)
\param   flag       INT              (i)    flag: submesh or sub-submesh?
\return void

------------------------------------------------------------------------*/
void f3_elesubmesh(ELEMENT        *ele,
                   FLUID_ML_SMESH *smesh,
		   INT             flag)
{
INT       k;          /* just a counter                  		*/
INT       nelbub;     /* number of bubble functions        		*/
INT       numnp;      /* number of nodal points 			*/
INT       ix,iy,iz;   /* counters in coordinate directions      	*/
INT       nx,ny,nz;   /* number of nodes in coordinate directions       */
INT       iel;        /* number of large-scale element nodes            */
DOUBLE    r[3];       /* coordinates on parent domain                   */
DOUBLE    funct[8];   /* large-scale element shape functions            */
DOUBLE    deriv_b[3*8];/* large-scale element shape function derivatives */
DOUBLE   *deriv[3];
DOUBLE    deriv2_b[6*8];/* l-s element shape function 2nd derivatives    */
DOUBLE   *deriv2[6];
DIS_TYP   typ;	      /* element type                                   */

#ifdef DEBUG
dstrc_enter("f3_elesubmesh");
#endif

/* Hack */
deriv[0] = &deriv_b[0];
deriv[1] = &deriv_b[8];
deriv[2] = &deriv_b[16];

deriv2[0] = &deriv2_b[0];
deriv2[1] = &deriv2_b[8];
deriv2[2] = &deriv2_b[16];
deriv2[3] = &deriv2_b[24];
deriv2[4] = &deriv2_b[32];
deriv2[5] = &deriv2_b[40];

/*------------------------------------------------------- define arrays */
switch (flag)
{
/*--------------------------------------------------------- sub-submesh */
case 1:
  (DOUBLE**)amdef("xyzssm",&(ele->e.f3->xyzssm),3,smesh->numnp,"DA");
  amzero(&(ele->e.f3->xyzssm));
break;
/*------------------------------------------------------------- submesh */
default:
  nelbub = 4*ele->numnp + 3;
  (DOUBLE**)amdef("xyzsm",&(ele->e.f3->xyzsm),3,smesh->numnp,"DA");
  (DOUBLE**)amdef("solsm" ,&(ele->e.f3->solsm) ,smesh->numeq,nelbub,"DA");
  (DOUBLE**)amdef("solsmn",&(ele->e.f3->solsmn),smesh->numeq,nelbub,"DA");
  amzero(&(ele->e.f3->xyzsm));
  amzero(&(ele->e.f3->solsm));
  amzero(&(ele->e.f3->solsmn));
break;
}

/*----------------------------------------------------------------------
  calculate coordinates of nodes on this particular large-scale element
------------------------------------------------------------------------*/
switch (flag)
{
/*--------------------------------------------------------- sub-submesh */
case 1:
  nx = nxss;
  ny = nyss;
  nz = nzss;
break;
/*------------------------------------------------------------- submesh */
default:
  nx = nxsm;
  ny = nysm;
  nz = nzsm;
break;
}
/*--- set large-scale element node number/type and submesh node counter */
iel  = ele->numnp;
typ  = ele->distyp;
numnp = 0;

for (iz=0; iz<nz; iz++)
{
  for (iy=0; iy<ny; iy++)
  {
    for (ix=0; ix<nx; ix++)
    {
/*----------------- coordinates of actual submesh node on parent domain */
      r[0] = smesh->xyzpd.a.da[0][numnp];
      r[1] = smesh->xyzpd.a.da[1][numnp];
      r[2] = smesh->xyzpd.a.da[2][numnp];
/*---------------------- evaluate shape functions for this submesh node */
      f3_hex(funct,deriv,deriv2,r[0],r[1],r[2],typ,1);
/*--------------------------------- compute and store nodal coordinates */
      switch (flag)
      {
/*--------------------------------------------------------- sub-submesh */
      case 1:
        for (k=0; k<iel; k++) /* loop all nodes of the l-s element */
        {
          ele->e.f3->xyzssm.a.da[0][numnp] += funct[k]*ele->node[k]->x[0];
          ele->e.f3->xyzssm.a.da[1][numnp] += funct[k]*ele->node[k]->x[1];
          ele->e.f3->xyzssm.a.da[2][numnp] += funct[k]*ele->node[k]->x[2];
        } /* end loop over iel */
      break;
/*------------------------------------------------------------- submesh */
      default:
        for (k=0; k<iel; k++) /* loop all nodes of the l-s element */
        {
          ele->e.f3->xyzsm.a.da[0][numnp] += funct[k]*ele->node[k]->x[0];
          ele->e.f3->xyzsm.a.da[1][numnp] += funct[k]*ele->node[k]->x[1];
          ele->e.f3->xyzsm.a.da[2][numnp] += funct[k]*ele->node[k]->x[2];
        } /* end loop over iel */
      break;
      }
      numnp++;
    }
  }
}

/*------------- set flag: this element submesh has already been created */
ele->e.f3->smisal = 1;

/*--------------------------------------------------------------------- */
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_elesubmesh */

#endif
