/*!----------------------------------------------------------------------
\file
\brief submesh creation for multilevel fluid2 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"

static INT nhorsm,nversm,nhorss,nverss;

/*!---------------------------------------------------------------------
\brief creation of (sub-)submesh on parent domain for fluid2

<pre>                                                       gravem 07/03

In this routine, the creation of submesh or sub-submesh on the parent
domain is performed, i.e. nodal coordinates on the parent domain,
id-array and ien-array are established.

</pre>
\param  *smesh      FLUID_ML_SMESH   (i/o)
\param   xele       INT              (i)    number of elements in x-dir.
\param   yele       INT              (i)    number of elements in y-dir.
\param   order      INT              (i)    polyn. interpolation order
\param   flag       INT              (i)    flag: submesh or sub-submesh?
\return void

------------------------------------------------------------------------*/
void f2_pdsubmesh(FLUID_ML_SMESH *smesh,
                  INT             xele,
                  INT             yele,
                  INT             order,
		  INT             flag)
{
DOUBLE hpdx,hpdy;  /* element length in coord. dir. on parent domain   */
INT    numeq;      /* number of equations                              */
INT    numnp;      /* number of nodal points                           */
INT    ihor,iver;  /* counters in horizontal/vertical direction        */
INT    nhor,nver;  /* number of nodes in horizontal/vertical direction */
INT    iel;        /* element counter                                  */
INT    nnglo;      /* counter for global node number                   */
INT    nnhor;      /* counter for number of nodes in a horizontal row  */
INT    nrhor;      /* counter for number of horizontal rows            */

#ifdef DEBUG
dstrc_enter("f2_pdsubmesh");
#endif

/*----------------------------------------------------------------------
  calculate coordinates of nodes on parent domain and evaluate id-array
------------------------------------------------------------------------*/
/*------------------------------------ divide large-scale element edges */
hpdx  = ONE/(xele*order);
hpdy  = ONE/(yele*order);
switch (flag)
{
/*--------------------------------------------------------- sub-submesh */
case 1:
  nhorss = order*xele+1;
  nverss = order*yele+1;
  nhor = nhorss;
  nver = nverss;
break;
/*------------------------------------------------------------- submesh */
default:
  nhorsm = order*xele+1;
  nversm = order*yele+1;
  nhor = nhorsm;
  nver = nversm;
break;
}
/*--------------------------------------- set node and equation counter */
numeq = 0;
numnp = 0;

for (iver=0; iver<nver; iver++)
{
  for (ihor=0; ihor<nhor; ihor++)
  {
/*-------------------------------- compute coordinates on parent domain */
    smesh->xyzpd.a.da[0][numnp] = ONE-TWO*ihor*hpdx;
    smesh->xyzpd.a.da[1][numnp] = ONE-TWO*iver*hpdy;
/*--------------------------------------------------- evaluate id-array */
    if (iver==0 || iver==nver-1 || ihor==0 || ihor==nhor-1)
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

/*----------------------------------------------------------------------
   evaluate ien-array
------------------------------------------------------------------------*/
/*----------------------------------------------------- linear elements */
if (order==1)
{
/*---------------------------------------------------- set row counters */
  nnhor = 1;
  nrhor = 1;
  for (iel=0; iel<smesh->numele; iel++) /* loop over submesh elements */
  {
    if (nnhor==nhor)
    {
/*-------------------------------------------- start new horizontal row */
      nnhor = 1;
      nrhor++;
    }
/*---------------------------------------- evaluate global node counter */
    nnglo = iel+nrhor-1;
/*----------------------------------------- assign element node numbers */
    smesh->ien.a.ia[iel][0] = nnglo;
    smesh->ien.a.ia[iel][1] = nnglo+1;
    smesh->ien.a.ia[iel][2] = nnglo+xele+2;
    smesh->ien.a.ia[iel][3] = nnglo+xele+1;
    nnhor++;
  }
}
/*------------------------------------------------- quadratic elements */
else if (order==2)
{
/*------------------------------------ set row and global node counter */
  nnhor = 1;
  nnglo = 0;
  for (iel=0; iel<smesh->numele; iel++) /* loop over submesh elements */
  {
    if (nnhor==nhor)
    {
/*-------------------------------------------- start new horizontal row */
      nnhor  = 1;
      nnglo += 2*xele+2;
    }
/*----------------------------------------- assign element node numbers */
    smesh->ien.a.ia[iel][0] = nnglo;
    smesh->ien.a.ia[iel][1] = nnglo+2;
    smesh->ien.a.ia[iel][2] = nnglo+4*xele+4;
    smesh->ien.a.ia[iel][3] = nnglo+4*xele+2;
    smesh->ien.a.ia[iel][4] = nnglo+1;
    smesh->ien.a.ia[iel][5] = nnglo+2*xele+3;
    smesh->ien.a.ia[iel][6] = nnglo+4*xele+3;
    smesh->ien.a.ia[iel][7] = nnglo+2*xele+1;
    smesh->ien.a.ia[iel][8] = nnglo+2*xele+2;
    nnhor += 2;
    nnglo += 2;
  }
}

/*--------------------------------------------------------------------- */
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_pdsubmesh */

/*!---------------------------------------------------------------------
\brief creation of (sub-)submesh on individual element for fluid2

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
void f2_elesubmesh(ELEMENT        *ele,
                   FLUID_ML_SMESH *smesh,
		   INT             flag)
{
INT       k;          /* just a counter                  		*/
INT       nelbub;     /* number of bubble functions        		*/
INT       numnp;      /* number of nodal points 			*/
INT       ihor,iver;  /* counters in horizontal/vertical direction	*/
INT       nhor,nver;  /* number of nodes in horizontal/vertical dir.    */
INT       iel;        /* number of large-scale element nodes            */
DOUBLE    r[2];       /* coordinates on parent domain                   */
DOUBLE    funct[9];   /* large-scale element shape functions            */
DOUBLE    deriv_b[2*9];/* large-scale element shape function derivatives */
DOUBLE   *deriv[2];
DOUBLE    deriv2_b[3*9];/* l-s element shape function 2nd derivatives    */
DOUBLE   *deriv2[3];
DIS_TYP   typ;	      /* element type                                   */

#ifdef DEBUG
dstrc_enter("f2_elesubmesh");
#endif

/* Hack */
deriv[0] = &deriv_b[0];
deriv[1] = &deriv_b[9];

deriv2[0] = &deriv2_b[0];
deriv2[1] = &deriv2_b[9];
deriv2[2] = &deriv2_b[18];

/*------------------------------------------------------- define arrays */
switch (flag)
{
/*--------------------------------------------------------- sub-submesh */
case 1:
  (DOUBLE**)amdef("xyzssm",&(ele->e.f2->xyzssm),2,smesh->numnp,"DA");
  amzero(&(ele->e.f2->xyzssm));
break;
/*------------------------------------------------------------- submesh */
default:
  nelbub = 3*ele->numnp + 2;
  (DOUBLE**)amdef("xyzsm" ,&(ele->e.f2->xyzsm) ,2,smesh->numnp     ,"DA");
  (DOUBLE**)amdef("solsm" ,&(ele->e.f2->solsm) ,smesh->numeq,nelbub,"DA");
  (DOUBLE**)amdef("solsmn",&(ele->e.f2->solsmn),smesh->numeq,nelbub,"DA");
  amzero(&(ele->e.f2->xyzsm));
  amzero(&(ele->e.f2->solsm));
  amzero(&(ele->e.f2->solsmn));
break;
}

/*----------------------------------------------------------------------
  calculate coordinates of nodes on this particular large-scale element
------------------------------------------------------------------------*/
switch (flag)
{
/*--------------------------------------------------------- sub-submesh */
case 1:
  nhor = nhorss;
  nver = nverss;
break;
/*------------------------------------------------------------- submesh */
default:
  nhor = nhorsm;
  nver = nversm;
break;
}
/*--- set large-scale element node number/type and submesh node counter */
iel  = ele->numnp;
typ  = ele->distyp;
numnp = 0;

for (iver=0; iver<nver; iver++)
{
  for (ihor=0; ihor<nhor; ihor++)
  {
/*----------------- coordinates of actual submesh node on parent domain */
    r[0] = smesh->xyzpd.a.da[0][numnp];
    r[1] = smesh->xyzpd.a.da[1][numnp];
/*---------------------- evaluate shape functions for this submesh node */
    f2_rec(funct,deriv,deriv2,r[0],r[1],typ,1);
/*--------------------------------- compute and store nodal coordinates */
    switch (flag)
    {
/*--------------------------------------------------------- sub-submesh */
    case 1:
      for (k=0; k<iel; k++) /* loop all nodes of the l-s element */
      {
        ele->e.f2->xyzssm.a.da[0][numnp] += funct[k]*ele->node[k]->x[0];
        ele->e.f2->xyzssm.a.da[1][numnp] += funct[k]*ele->node[k]->x[1];
      } /* end loop over iel */
    break;
/*------------------------------------------------------------- submesh */
    default:
      for (k=0; k<iel; k++) /* loop all nodes of the l-s element */
      {
        ele->e.f2->xyzsm.a.da[0][numnp] += funct[k]*ele->node[k]->x[0];
        ele->e.f2->xyzsm.a.da[1][numnp] += funct[k]*ele->node[k]->x[1];
      } /* end loop over iel */
    break;
    }
    numnp++;
  }
}

/*------------- set flag: this element submesh has already been created */
ele->e.f2->smisal = 1;

/*--------------------------------------------------------------------- */
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_elesubmesh */

#endif
