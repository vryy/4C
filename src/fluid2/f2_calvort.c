/*!----------------------------------------------------------------------
\file
\brief calculation of vorticity

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

static DOUBLE Q12=(ONE/TWO);

/*!---------------------------------------------------------------------
\brief routine to calculate the vorticity

<pre>                                                         genk 12/02
   calculation of vorticity (for visualisation):
      vort = 1/2 (Ux,y + Uy,x)
   when extrapolating the vorticity from the gauss points to the nodes
   the results are averaged with the number of elements the actual node
   belongs to.
   An alternative would be to average with the element areas;
   smaller elements = better values!
</pre>

\param  *ele	            ELEMENT	     (i)   actual element
\param   init	            INT              (i)   init flag
\return void

------------------------------------------------------------------------*/
void f2_calvort(
	        ELEMENT        *ele,
       	        INT             init
               )
{
static ARRAY     evel_a;     /* element velocities at (n)                 */
static DOUBLE  **evel;
static ARRAY     funct_a;    /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;    /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     xjm_a;      /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     vderxy_a;   /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     derxy_a;    /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vort_a;     /* vorticity                                 */
static DOUBLE   *vort;
static ARRAY     xyze_a;     /* element coordinates                       */
static DOUBLE  **xyze;

INT j,icol,kk;
INT iv=0;		     /* index for gauss points                    */
INT ivmax;		     /* index for maximum number of gauss points  */
INT icode=2;                 /* flag for eveluation of shape functions    */
INT iel=0;		     /* number of nodes                           */
INT lr,ls;		     /* index to loop over integration points     */
INT node;		     /* index to loop over nodal points           */
INT nir,nis;                 /* number of gauss points in r&s directions  */
INT intc;
INT ncols;                   /* number of columns                         */
INT numele;                  /* number of elements to the actual element  */
DOUBLE det;                  /* determinant                               */
DOUBLE e1,e2;                /* coordinates of the current gauss points   */
DOUBLE r;                    /* local coord. of gauss points in r-dir.    */
DOUBLE s;                    /* local coord. of gauss points in s-dir     */
DOUBLE f;                    /* vorticity value at the nodes              */
DOUBLE fpar[MAXGAUSS];
NODE *actnode;               /* actual node                               */
DIS_TYP typ;                 /* element type                              */
static FLUID_DATA     *data;

#ifdef DEBUG
dstrc_enter("f2_calvort");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   evel    = amdef("evel"   ,&evel_a  ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   funct   = amdef("funct"  ,&funct_a  ,MAXNOD_F2,1,"DV");
   deriv   = amdef("deriv"  ,&deriv_a  ,2,MAXNOD_F2,"DA");
   xjm     = amdef("xjm"    ,&xjm_a    ,2,2        ,"DA");
   vderxy  = amdef("vderxy" ,&vderxy_a ,2,2        ,"DA");
   derxy   = amdef("derxy"  ,&derxy_a  ,2,MAXNOD_F2,"DA");
   vort    = amdef("vort"   ,&vort_a   ,MAXGAUSS,1 ,"DV");
   xyze    = amdef("xyze"   ,&xyze_a   ,2,MAXNOD_F2,"DA");

   fdyn    = alldyn[genprob.numff].fdyn;
   data    = fdyn->data;
   goto end;
} /* endif (init==1) */

/*------------------------------------------------- calculate vorticity */
/*	                vort = 1/2*(Ux,y - Uy,x)                        */
/*----------------------------------------------------------------------*/

/*-------------------------------------------- loop over each time step */
ncols = fdyn->ncols;
for (icol=0;icol<ncols;icol++)
{

/*------------------------------------------------------ initialization */
iel = ele->numnp;
typ = ele->distyp;

/*------------------------------------------------- get integraton data */
switch (typ)
{
case quad4: case quad8: case quad9: /* -----> quad element */
   icode = 2;
	/* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case tri3: case tri6: /* triangular element */
      icode = 2;
      nir  = ele->e.f2->nGP[0];
      nis  = 1;
      intc = ele->e.f2->nGP[1];
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */


/*---------------------- set element velocities, and coordinates -------*/
for (j=0;j<iel;j++)
{
   evel[0][j] = ele->node[j]->sol.a.da[icol][0];
   evel[1][j] = ele->node[j]->sol.a.da[icol][1];
   xyze[0][j] = ele->node[j]->x[0];
   xyze[1][j] = ele->node[j]->x[1];
}/*end for (j=0;j<iel;j++) */

/*---------------------------------- start loop over integration points */
iv = 0;
for (lr=0;lr<nir;lr++)
{
   for (ls=0;ls<nis;ls++)
      {
/*---------------- get values of  shape functions and their derivatives */
      switch(typ)
      {
      case quad4: case quad8: case quad9:  /* quad element*/
         e1 = data->qxg[lr][nir-1];
         e2 = data->qxg[ls][nis-1];
         f2_rec(funct,deriv,NULL,e1,e2,typ,icode);
      break;
      case tri3: case tri6:   /* --> tri - element */
	 e1   = data->txgr[lr][intc];
	 e2   = data->txgs[lr][intc];
	 f2_tri(funct,deriv,NULL,e1,e2,typ,icode);
      break;
      default:
         dserror("typ unknown!");
      }/*end of switch(typ) */

/*--------------------------------------------- compute Jacobian matrix */
      f2_jaco2(xyze,deriv,xjm,&det,iel);
/*--------------------------------- compute global derivative ----------*/
      f2_gder(derxy,deriv,xjm,det,iel);
/*-------- compute velocity derivatives at integration points ----------*/
      f2_vder(vderxy,derxy,evel,iel);
/*----------------------- calculate vorticity --------------------------*/
      vort[iv] = Q12*(vderxy[0][1]-vderxy[1][0]);
      iv++;
   }/* end of for (ls=0;ls<nis;ls++) */
}/* end of for (lr=0;lr<nir;lr++) */

ivmax = iv;

/*------------------------------start loop over nodal points------------*/
for (node=0;node<iel;node++)
{
   actnode = ele->node[node];                          /*--actual node--*/
   numele  = actnode->numele; /*--number of elements to the actual node-*/
   if (actnode->sol.sdim < 4)
   {
       amredef(&(actnode->sol),actnode->sol.fdim,4,"DA");
       for (kk=0;kk<ncols;kk++)
       actnode->sol.a.da[kk][3] = ZERO;
   }
/*------ --------get local coordinates of nodes-------------------------*/
   r = f2_rsn(node,0,iel);
   s = f2_rsn(node,1,iel);
/*---------------------------------- extrapolate vorticity to the nodes */
   switch (typ)
   {
   case quad4: case quad8: case quad9:
      f2_recex(&f,fpar,r,s,vort,ivmax,1);
   break;
   case tri3: case tri6:
      f2_triex(&f,fpar,r,s,vort,ivmax,1);
   break;
   default:
      dserror("typ unknown!\n");
   }
/*------------------- store the vorticity value in the solution history-*/
   actnode->sol.a.da[icol][3] += f/numele;

}/* end of loop for (node=0;node<iel;node++) */
}/* end of loop over columns of solution history */


end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calvort */

#endif
/*! @} (documentation module close)*/
