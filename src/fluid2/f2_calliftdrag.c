/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/

#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*---------------------------------------------------------------------*/

/*!---------------------------------------------------------------------
\brief calculation of Lift and drag total loads

<pre>                                         genk & Mustafa Omercik 6/03

in this function element Total forces regarding to the nodal\
stress values calculated


*************************************************************************
</pre>

------------------------------------------------------------------------*/

void f2_calliftdrag(	ELEMENT       *ele,
			FLUID_DATA    *data,
			CONTAINER     *container)
{
INT		lr,ls;		/* INTegration directions		*/
INT		i,j,jj,k;	/* some loopers				*/
INT		nir;		/* number of GP's in r direction 	*/
INT		ngline;		/* number of geometrylines at element 	*/   
INT		line;		/* looper over lines			*/
INT		ngnode;		/* number of geometry-nodes on g-line	*/
INT		ngr,ngs;	/* number of GP for line-integration	*/
INT		iegnod[MAXNOD_F2]; /* (o)   edge nodes*/
const INT	numdf  = 3;	/* dof per node				*/

DOUBLE		vnx,vny;	/* comp, of normal vector at int point	*/ 
DOUBLE		forceline[2];/* lineload value in x and y direct.(inp)	*/
DOUBLE		sigmaint[3]; /* fluid stresses at integration point	*/
DOUBLE		nsigma[3][MAXNOD_F2]; /* nodal fluid stresses		*/
DOUBLE		xyzl[2][MAXNOD_F2]; /* nodal coordinates		*/
DOUBLE		e1;		/* GP-koordinates in r-s-system		*/
DOUBLE		facr;		/* integration factor  GP-info		*/
DOUBLE		det;		/* det of jacobian matrix		*/

GLINE		*gline[4];	/* geometrylines of the element		*/
GLINE		*actgline;
DIS_TYP		typ;		/* element type				*/
GNODE		gnode;
NODE		*actnode;

static DOUBLE	*funct;		/* ansatz-functions			*/
static ARRAY	funct_a;
static DOUBLE	**deriv;	/* derivatives of ansatz-funct		*/
static ARRAY	deriv_a;
static DOUBLE	**xjm;		/* jacobian matrix			*/
static ARRAY	xjm_a;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("f2_calinta");
#endif

/*--------------------------------------------------- initialization ---*/
 funct   = amdef("funct"  ,&funct_a,MAXNOD_F2,1            ,"DV");       
 deriv   = amdef("deriv"  ,&deriv_a,2            ,MAXNOD_F2,"DA");       
 xjm     = amdef("xjm_a"  ,&xjm_a  ,numdf     ,numdf       ,"DA");

amzero(&funct_a);
amzero(&deriv_a);
amzero(&xjm_a);

/*-------------------------------------------- init the gaussian points */
f2_intg(data,1);

nir     = ele->e.f2->nGP[0];
nir = IMAX(nir,2);

ngline=ele->g.gsurf->ngline;
typ=ele->distyp;

/*-------------------------------------------------- loop over Gline ---*/
for (line=0; line<ngline; line++) 
{
   actgline = ele->g.gsurf->gline[line];
   /*- confirm that the gline is the right one that lies on the body ---*/
   if (actgline->dline==NULL) continue;
   if (actgline->dline->liftdrag==0) continue;

   /*------------- interpolate the stress tensor to the Gauss points ---*/
   /*------------------------------------------------ get edge nodes ---*/
   f2_iedg(iegnod,ele,line,0);
   
   /*--------------------------------------------- read nodal values ---*/
   ngnode = actgline->ngnode;
   for (j=0; j<ngnode; j++)
   {
      for (i=0; i<3; i++)		/* stresses			*/
	 nsigma[i][j] = ele->e.f2->stress_ND.a.da[iegnod[j]][i];

      actnode=ele->node[iegnod[j]];	/* coordinates			*/
      xyzl[0][j] = actnode->x[0];
      xyzl[1][j] = actnode->x[1];
   }

   /*======================= integration loop ==========================*/
   /*------------- over liftdrag edge of actual element ----------------*/

   for (lr=0; lr<nir; lr++)
   {
      /*--------------------------------- Gauss point at weight at it ---*/
      e1   = data->qxg[lr][nir-1];
      facr = data->qwgt[lr][nir-1];
      /*-------- get values of  shape functions and their derivatives ---*/	
      f2_degrectri(funct,deriv,e1,typ,1);

      /*-------------------- interpolate stress tensor to Gauss point ---*/
      for (i=0;i<3;i++) sigmaint[i] = ZERO;
     
      for (i=0;i<3;i++)
       	 for (j=0;j<ngnode;j++)  sigmaint[i] += funct[j] * nsigma[i][j];

      /*-------------------------- compute normal vector at gauss point *
       |	       see Mok et al in Engineering Computations (1999) |
       *----------------------------------------------------------------*/
      vnx=ZERO;
      vny=ZERO;
      for(i=0;i<ngnode;i++) vnx+=deriv[0][i]*xyzl[1][i];
      for(i=0;i<ngnode;i++) vny-=deriv[0][i]*xyzl[0][i];
      
      /*-------------------------- compute stress vector at gauss point *
       |  force = sigma * n  (Cauchy's law)                             |
       *----------------------------------------------------------------*/
      for (i=0;i<2;i++) forceline[i] = ZERO;
      forceline[0] = sigmaint[0]*vnx + sigmaint[2]*vny;
      forceline[1] = sigmaint[2]*vnx + sigmaint[1]*vny;
      
      /*-------------- add load vector component to element load vector */
      /* jacobian determinant cancels with length of normal  vector	*/
      for (j=0; j<ngnode; j++)
      {
         /*-- TOTAL Lift and Drag forces in global coordinate system ---*/
         container->liftdrag[0] -= funct[j] * forceline[0] * facr;
         container->liftdrag[1] += funct[j] * forceline[1] * facr;
      }/*----------------------------------------- end loop over ngnode */
   }/*==================================== end integration loop over lr */
}/* end loop line over lines to this element */

/*---------------------------------------------------------- tidy up ---*/
amdel(&funct_a);
amdel(&deriv_a);
amdel(&xjm_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calinta */



void f2_plot_liftdrag(DOUBLE time, DOUBLE *liftdrag)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("f2_plot_liftdrag");
#endif

fprintf(allfiles.gnu,"%10.5f  %8.7f  %8.7f\n", 
        time, liftdrag[0], liftdrag[1]);
fflush(allfiles.gnu);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_plot_liftdrag */

#endif
/*! @} (documentation module close)*/         
