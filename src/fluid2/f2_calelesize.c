/*!----------------------------------------------------------------------
\file
\brief Calculate stabilisation parameter

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
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
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

/*!---------------------------------------------------------------------
\brief routine to calculate element size and stabilisation parameter

<pre>                                                         genk 04/02

   ele->e.f2->stabi.gls->iadvec: adevction stab.
      0 = no
      1 = yes
   ele->e.f2->stabi.gls->ipres: pressure stab.
      0 = no
      1 = yes
   ele->e.f2->stabi.gls->ivisc: diffusion stab.
      0 = no
      1 = GLS-
      2 = GLS+
   ele->e.f2->stabi.gls->icont: continuity stab.
      0 = no
      1 = yes
   ele->e.f2->stabi.gls->istapa: version of stab. parameter
      35 = diss wall instationary
      36 = diss wall stationanary
   ele->e.f2->stabi.gls->norm_P: p-norm
      p = 1<=p<=oo
      0 = max.-norm (p=oo)
   ele->e.f2->stabi.gls->mk: higher order elements control flag
      0 = mk fixed (--> (bi)linear: 1/3, biquadr.: 1/12)
      1 = min(1/3,2*C)
     -1 = mk=1/3  (--> element order via approx. nodal-distance)
   ele->e.f2->stabi.gls->ihele[]:
      x/y/z = length-def. for velocity/pressure/continuity stab
      0 = don't compute
      1 = sqrt(area)
      2 = area equivalent diameter
      3 = diameter/sqrt(2)
      4 = sqrt(2)*area/diagonal (rectangle) 4*area/s (triangle)
      5 = streamlength (element length in flow direction
   ele->e.f2->stabi.gls->ninths: number of integration points for streamlength
      1 = at center of element
      2 = at every INT pt used for element.-stab.-matrices
   ele->e.f2->stabi.gls->istapc: flag for stabilisation parameter calculation
      1 = at center of element
      2 = at every integration point
   ele->e.f2->stabi.gls->clamb \
   ele->e.f2->c1               |_>> stabilisation constants (input)
   ele->e.f2->c2               |
   ele->e.f2->c3              /
   ele->e.f2->stabi.gls->istrle: has streamlength to be computed
   ele->e.f2->stabi.gls->iareavol: calculation of area length
   ele->e.f2->stabi.gls->iduring: calculation during INT.-pt.loop
   ele->e.f2->stabi.gls->itau[0]: flag for tau_mu calc. (-1: before, 1:during)
   ele->e.f2->stabi.gls->itau[1]: flag for tau_mp calc. (-1: before, 1:during)
   ele->e.f2->stabi.gls->itau[2]: flag for tau_c calc. (-1: before, 1:during)
   ele->e.f2->hk[i]: "element sizes" (vel / pre / cont)
   ele->e.f2->stabi.gls->idiaxy: has diagonals to be computed
   fdyn->tau[0]: stability parameter momentum / velocity (tau_mu)
   fdyn->tau[1]: stability parameter momentum / pressure (tau_mp)
   fdyn->tau[2]: stability parameter continuity (tau_c)
</pre>
\param  *ele     ELEMENT	       (i)   actual element
\param  *eleke   ELEMENT	       (i)   actual element (only for turbulence)
\param **xzye    DOUBLE                (-)   nodal coordinates
\param  *funct   DOUBLE 	       (-)   shape functions
\param **deriv   DOUBLE 	       (-)   deriv. of shape funcs
\param **deriv2  DOUBLE 	       (-)   2nd deriv. of sh. funcs
\param **xjm     DOUBLE 	       (-)   jacobian matrix
\param **derxy   DOUBLE 	       (-)   global derivatives
\param **vderxy  DOUBLE 	       (-)   global derivatives of velocity
\param **evel    DOUBLE 	       (i)   element velocities
\param  *velint  DOUBLE 	       (-)   vel. at integr. point
\param  *velint  DOUBLE 	       (-)   vel. at integr. point
\param  *eddyint  DOUBLE 	       (-)   eddy-visc. at integr. point (only for turbulence)
\param  *visc     DOUBLE 	       (-)   viscosity
\return void

------------------------------------------------------------------------*/
void f2_calelesize(
                     ELEMENT         *ele,
                     ELEMENT         *eleke,
                     DOUBLE         **xyze,
                     DOUBLE          *funct,
                     DOUBLE         **deriv,
                     DOUBLE         **deriv2,
                     DOUBLE         **xjm,
                     DOUBLE         **derxy,
                     DOUBLE         **vderxy,
                     DOUBLE         **evel,
                     DOUBLE          *ephin,
                     DOUBLE          *ephing,
                     DOUBLE          *eddy,
                     DOUBLE          *visc,
                     INT              cpele
                     )
{

INT     i,ilen;         /* simply a counter	        		*/
INT     ieval = 0;	/* evaluation flag			        */
INT     igc   = 0;	/* evaluation flag			        */
INT     istrnint;       /* evaluation flag			        */
INT     isharea;        /* evaluation flag			        */
INT     actmat;         /* number of actual material		        */
INT     iel;            /* number of nodes of actual element            */
DOUBLE  area;           /* element area                                 */
DOUBLE  det;            /* determinant of jacobian                      */
DOUBLE  strle;          /* streamlength                                 */
DOUBLE  e1,e2;          /* natural coordinates of inegration point      */
DOUBLE  fac,facr,facs;  /* factors                                      */
DOUBLE  dia,dia1,dia2;  /* values used for calculation of element size  */
DOUBLE  dx,dy;          /* values used for calculation of element size  */
DOUBLE  gcoor[2];       /* global coordinates                           */
DOUBLE  velint[2];
DIS_TYP typ;
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/
FLUID_DYNAMIC *fdyn;
FLUID_DATA      *data;
#ifdef D_FLUID2TU
DOUBLE  eddyint;        /* eddy-viscosity                               */
#endif


#ifdef DEBUG
dstrc_enter("f2_calelesize");
#endif

/*---------------------------------------------------------- initialise */
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
iel    = ele->numnp;
typ    = ele->distyp;
gls    = ele->e.f2->stabi.gls;

dsassert(ele->e.f2->stab_type == stab_gls,
         "routine with no or wrong stabilisation called");

istrnint = gls->istrle * gls->ninths;
isharea  = fdyn->ishape * gls->iareavol;

/*----------------------------------------------------------------------*
 | calculations at element center: area & streamlength                  |
 | NOTE:                                                                |
 |    area is always calculated using only 1 integrationpoint           |
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,facr,facs with their constant values in the calls of   |
 |         f2_rec / f2_tri!!!!!!                                        |
 *----------------------------------------------------------------------*/

if (isharea==1)
{
   area  = ZERO;
   strle = ZERO;
/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(typ)
   {
   case quad4: case quad8: case quad9:    /* --> quad - element */
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      f2_rec(funct,deriv,deriv2,e1,e2,typ,2);
   break;
   case tri3: case tri6:       /* --> tri - element */
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      f2_tri(funct,deriv,deriv2,e1,e2,typ,2);
   break;
   default:
      dserror("typ unknown!\n");
   } /*end switch(typ) */
   ieval++;
/* -------------------------------------------- compute jacobian matrix */
   f2_jaco(xyze,deriv,xjm,&det,iel,ele);
   fac=facr*facs*det;
   area += fac;
   fdyn->totarea += area;
   if (istrnint==1)    /* compute streamlength */
   {
      f2_veci(velint,funct,evel,iel);
      ieval++;
      f2_gcoor(xyze,funct,iel,gcoor);
      igc++;
      f2_calstrlen(&strle,xyze,velint,ele,gcoor,typ);
   } /* enidf (istrnint==1) */
   if (gls->idiaxy==1)    /* compute diagonal based diameter */
   {
      switch(typ)
      {
      case quad4: case quad8: case quad9:
         dx = xyze[0][0] - xyze[0][2];
	 dy = xyze[1][0] - xyze[1][2];
	 dia1 = sqrt(dx*dx+dy*dy);
	 dx = xyze[0][1] - xyze[0][3];
	 dy = xyze[1][1] - xyze[1][3];
	 dia2 = sqrt(dx*dx+dy*dy);
/*------ dia=sqrt(2)*area/(1/2*(dia1+dia2))=sqrt(8)*area/(dia1+dia2) ---*/
	 dia = sqrt(EIGHT)*area/(dia1+dia2);
      break;
      case tri3: case tri6:    /* get global coordinate of element center */
         if (igc==0)
	    f2_gcoor(xyze,funct,iel,gcoor);
	 dia = ZERO;
	 for (i=0;i<3;i++)
	 {
	    dx = gcoor[0] - xyze[0][i];
	    dy = gcoor[1] - xyze[1][i];
	    dia += dx*dx + dy*dy;
	 }
	 dia = FOUR*area/sqrt(THREE*dia);
      break;
      default:
          dserror("typ unknown!\n");
      } /* end switch(typ) */
   } /* endif (ele->e.f2->idiaxy==1) */
/*--------------------------------------------------- set element sizes *
  ----loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for(ilen=0;ilen<3;ilen++)
   {
      if (gls->ihele[ilen]==1)
         ele->e.f2->hk[ilen] = sqrt(area);
      else if (gls->ihele[ilen]==2)
         ele->e.f2->hk[ilen] = TWO*sqrt(area/PI);
      else if (gls->ihele[ilen]==3)
         ele->e.f2->hk[ilen] = sqrt(TWO*area/PI);
      else if (gls->ihele[ilen]==4)
         ele->e.f2->hk[ilen] = dia;
      else if (gls->ninths==1)
         ele->e.f2->hk[ilen] = strle;
   } /* end loop over ilen */
} /* endif (isharea==1) */

/*----------------------------------------------------------------------*
 | calculations at element center: only streamlength                    |
 |    --> it may be possible to save some operations here by replacing  |
 |         e1,e2,facr,facs with their constant values in the calls of   |
 |         f2_rec / f2_tri!!!!!!                                        |
 *----------------------------------------------------------------------*/
else if (istrnint==1 && isharea !=1)
{
   area  = ZERO;
   strle = ZERO;
/*------ get values of integration parameters, shape functions and their
         derivatives ---------------------------------------------------*/
   switch(typ)
   {
   case quad4: case quad8: case quad9:    /* --> quad - element */
      e1   = data->qxg[0][0];
      facr = data->qwgt[0][0];
      e2   = data->qxg[0][0];
      facs = data->qwgt[0][0];
      f2_rec(funct,deriv,deriv2,e1,e2,typ,2);
   break;
   case tri3: case tri6:       /* --> tri - element */
      e1   = data->txgr[0][0];
      facr = data->twgt[0][0];
      e2   = data->txgs[0][0];
      facs = ONE;
      f2_tri(funct,deriv,deriv2,e1,e2,typ,2);
   break;
   default:
      dserror("typ unknown!\n");
   } /* end switch(typ) */
   ieval++;
/* ------------------------------------------- compute jacobian matrix */
   f2_jaco(xyze,deriv,xjm,&det,iel,ele);
/*----------------------------------------------- compute streamlength */
   f2_veci(velint,funct,evel,iel);
   ieval++;
   f2_gcoor(xyze,funct,iel,gcoor);
   igc++;
   f2_calstrlen(&strle,xyze,velint,ele,gcoor,typ);
/*--------------------------------------------------- set element sizes *
      loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for (ilen=0;ilen<3;ilen++)
   {
      if (gls->ihele[ilen]==5)
         ele->e.f2->hk[ilen] = strle;
   } /* end loop over ilen */
} /* endif (istrnint==1 && isharea !=1) */

/*----------------------------------------------------------------------*
  calculate stabilisation parameter
 *----------------------------------------------------------------------*/
if(gls->istapc==1 || istrnint==1)
{
   switch(ieval) /* ival>2: vel at intpoint already available! */
   {
   case 0:
/*------ get only values of integration parameters and shape functions
        no derivatives -------------------------------------------------*/
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
         e1   = data->qxg[0][0];
         facr = data->qwgt[0][0];
         e2   = data->qxg[0][0];
         facs = data->qwgt[0][0];
         f2_rec(funct,deriv,deriv2,e1,e2,typ,1);
      break;
      case tri3: case tri6:       /* --> tri - element */
         e1   = data->txgr[0][0];
         facr = data->twgt[0][0];
         e2   = data->txgs[0][0];
         facs = ONE;
         f2_tri(funct,deriv,deriv2,e1,e2,typ,1);
      break;
      default:
         dserror("typ unknown!\n");
      } /* end switch(typ) */
      f2_veci(velint,funct,evel,iel);
   break;
   case 1:
      f2_veci(velint,funct,evel,iel);
   break;
   case 2:
   break;
   default:
      dserror("wrong value for ieval\n");
   } /* end swtich(ieval) */
/*----------------------------------- calculate stabilisation parameter */
   actmat=ele->mat-1;
   (*visc) = mat[actmat].m.fluid->viscosity;

#ifdef D_FLUID2TU
  if (ele->e.f2->turbu == 1)
  {
   /*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
   /*---------------------------------- compute global derivates for vel. */
      f2_vder(vderxy,derxy,evel,iel);

     (*visc) += f2_calvisc(ele,vderxy);
  }

  if (ele->e.f2->turbu == 2 || ele->e.f2->turbu == 3)
  {
   f2_eddyirans(eleke,&eddyint,funct,eddy,iel);
   (*visc) += eddyint;
  }
#endif

   f2_calstabpar(ele,velint,(*visc),iel,typ,-1);
} /* endif (ele->e.f2->istapc==1 || istrnint==1) */

/*--------------------------------------------- copy stabpar to element */
if (cpele==1)
{
   for (i=0;i<3;i++)
      ele->e.f2->tau_old.a.dv[i]=fdyn->tau[i];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calelesize */

/*!---------------------------------------------------------------------
\brief routine to calculate element size and stabilisation parameter

<pre>                                                         genk 04/02

in this routine the element size and the stabilisation parameter
is calculated for one element during the integration loop

</pre>
\param  *ele     ELEMENT	        (i)    actual element
\param  *xyze    DOUBLE                 (-)    nodal coordinates
\param  *funct   DOUBLE 		(-)    natural shape funcs
\param  *velint  DOUBLE 		(-)    vel at intpoint
\param   visc    DOUBLE 		(i)    fluid viscosity
\param   iel     INT		        (i)    act. num. of ele nodes
\param   typ     DIS_TYP                (i)    element type
\return void
\sa f2_calelesize()

------------------------------------------------------------------------*/
void f2_calelesize2(
                     ELEMENT         *ele,
                     DOUBLE         **xyze,
                     DOUBLE          *funct,
                     DOUBLE          *velint,
                     DOUBLE           visc,
                     INT              iel,
                     DIS_TYP          typ
                    )
{
INT    ilen;       /* simply a counter                                  */
INT    istrnint;   /* evaluation flag                                   */
DOUBLE strle;      /* stream length                                     */
DOUBLE gcoor[2];   /* global coordinates                                */
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG
dstrc_enter("f2_calelesize2");
#endif

/*---------------------------------------------------------- initialise */
gls    = ele->e.f2->stabi.gls;
istrnint = gls->istrle * gls->ninths;

dsassert(ele->e.f2->stab_type == stab_gls,
   "routine with no or wrong stabilisation called");

if (istrnint==2)
{
/*------------------------------------------------ compute streamlength */
   f2_gcoor(xyze,funct,iel,gcoor);
   f2_calstrlen(&strle,xyze,velint,ele,gcoor,typ);
/*--------------------------------------------------- set element sizes *
      loop over 3 different element sizes: vel/pre/cont  ---------------*/
   for (ilen=0;ilen<3;ilen++)
   {
      if (gls->ihele[ilen]==5)
         ele->e.f2->hk[ilen] = strle;
   } /* end loop over ilen */
} /* endif (istrnint==2) */

/*----------------------------------- calculate stabilisation parameter */
f2_calstabpar(ele,velint,visc,iel,typ,1);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calelesize2 */

/*!---------------------------------------------------------------------
\brief routine to calculate streamlength

<pre>                                                         genk 04/02

in this routine the the streamlength, used for calculation of
stabilisation parameter is calculated. for higher order element this
is only an approximation, since the boundaries are assumed to be
straight.

</pre>
\param  *strle     DOUBLE   (o)    streamlength
\param  *velint    DOUBLE   (i)    velocities at integr. point
\param **xyze      DOUBLE   (i)    nodal coordinates
\param  *ele       ELEMENT  (i)    actual element
\param  *gcoor     DOUBLE   (i)    global coord. of INT. point
\param   typ       DIS_TYP  (i)    flag for element type
\return void
\sa f2_calelesize()

------------------------------------------------------------------------*/
void f2_calstrlen(
                  DOUBLE   *strle,
                  DOUBLE  **xyze,
                  DOUBLE   *velint,
                  ELEMENT  *ele,
                  DOUBLE   *gcoor,
                  DIS_TYP   typ
                 )
{
char string[100]="Couldn't find two cutting points for element             \n";
INT     nodcut=-1;
INT     nodmax;
INT     inod;
DOUBLE cutp[2][2];
DOUBLE dl,dx,dy,dxh,dyh;
DOUBLE dsub,dval;

#ifdef DEBUG
dstrc_enter("f2_calstrlen");
#endif

dval = FABS(velint[0])+FABS(velint[1]);
if (dval == ZERO)  /* no flow at this point - take some arbitr. measure for streamlength */
{
   dx = xyze[0][2] - xyze[0][0];
   dy = xyze[1][2] - xyze[1][0];
   goto calc2;
} /* enidf (dval == ZERO) */

/*----------------------------------------------------------------------*
   streamlength is calculated via cutting points of velocity vector
   with straight boundaries
*/
switch(typ)
{
case quad4: case quad8: case quad9:
   /* max number of nodes for quad: 4 --> C-numbering nodmax = 4-1 */
   nodmax = 3;
break;
case tri3: case tri6:
   /* max number of nodes for tri: 3 --> C-numbering nodmax = 3-1 */
   nodmax = 2;
break;
default:
   dserror("typ unknown!\n");
} /* end switch(typ) */
 /*------------------------------------------------- get cutting points */
for (inod=0;inod<nodmax;inod++)
{
   dxh = xyze[0][inod+1] - xyze[0][inod];
   dyh = xyze[1][inod+1] - xyze[1][inod];
   dsub = dxh*velint[1]-dyh*velint[0];
   if (dsub==ZERO)  /* check for parallel vectors */
      continue;
   dl = ((xyze[1][inod]-gcoor[1])*velint[0] -	\
	 (xyze[0][inod]-gcoor[0])*velint[1])/dsub;
   if (dl>=ZERO && dl<=ONE)
   {
      nodcut++;
      cutp[0][nodcut]=xyze[0][inod]+dl*dxh;
      cutp[1][nodcut]=xyze[1][inod]+dl*dyh;
      if (nodcut==1)
	 goto calc1;
   } /* endif (dl>=ZERO && dl<=ONE) */
} /* end loop over inod */
/*------------------------------------------------- the last boundary */
dxh = xyze[0][0]-xyze[0][nodmax];
dyh = xyze[1][0]-xyze[1][nodmax];
dsub = dxh*velint[1] - dyh*velint[0];
if (dsub==ZERO)
   dserror("Couldn't find two cutting points!\n");
dl = ((xyze[1][nodmax]-gcoor[1])*velint[0] -	\
      (xyze[0][nodmax]-gcoor[0])*velint[1])/dsub;
if (dl>=ZERO && dl <= ONE)
{
   nodcut++;
   cutp[0][nodcut]=xyze[0][nodmax]+dl*dxh;
   cutp[1][nodcut]=xyze[1][nodmax]+dl*dyh;
   if(nodcut==1)
      goto calc1;
} /* endif  (dl>=ZERO && dl <= ONE) */

sprintf(&string[45] ,"%-d" ,ele->Id);
dserror(string);

calc1:
dx = cutp[0][1]-cutp[0][0];
dy = cutp[1][1]-cutp[1][0];

calc2:
*strle = sqrt(dx*dx+dy*dy);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calstrlen */
#endif
/*! @} (documentation module close)*/
