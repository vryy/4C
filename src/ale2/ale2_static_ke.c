/*!----------------------------------------------------------------------
\file
\brief contains classic routines which integrate the stiffness for the 2d
ale element

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
#include "ale2.h"

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element

<pre>                                                              mn 06/02
This routine integrates the linear stiffness for the 2d ale element

</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE2_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void
\sa calling: ale2_intg(), ale2_funct_deriv(), ale2_jaco(), ale2_bop(),
             ale2_mat_linel(), ale2_keku(), ale2_hourglass();
             called by: ale2()

*----------------------------------------------------------------------*/
void ale2_static_ke(ELEMENT   *ele,
                   ALE2_DATA  *data,
                   MATERIAL   *mat,
                   ARRAY      *estif_global,
                   INT         init)
{
INT                 i, j;         /* counters */
INT                 nir=0,nis=0;  /* num GP in r/s/t direction */
INT                 lr, ls;       /* loopers over GP */
INT                 iel;          /* numnp to this element */
INT                 nd;
const INT           numdf =2;
const INT           numeps=3;

DOUBLE              fac;
DOUBLE              e1,e2=0;      /*GP-coords*/
DOUBLE              facr,facs=0;  /* weights at GP */

static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    xyz_a;    /* actual element coordiantes */
static DOUBLE **xyz;
static DOUBLE **estif;    /* element stiffness matrix ke */

DOUBLE det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_static_ke");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");
deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");
D         = amdef("D"      ,&D_a   ,6,6              ,"DA");
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");
xyz       = amdef("xyz"    ,&xyz_a ,4    ,numdf      ,"DA");

bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");
goto end;
}
/*------------------------------------------- integration parameters ---*/
ale2_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
switch (ele->distyp)
{
    case quad4:
    case quad8:
    case quad9:
        nir     = ele->e.ale2->nGP[0];
        nis     = ele->e.ale2->nGP[1];
        break;
    case tri3:
    case tri6:
        nir = ele->e.ale2->nGP[0];
        nis = 1;
        break;
    default:
        dserror("unknown number of gaussian points in ale2_intg");
        break;
}

iel     = ele->numnp;
nd      = numdf * iel;
/*--------------------------------------- actual element coordinates ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
         xyz[i][j] = ele->node[i]->x[j];
   }
}
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     /*============================= gaussian point and weight at it ===*/
      switch (ele->distyp)
      {
          case quad4:
          case quad8:
          case quad9:
              e2   = data->xgps[ls];
              facs = data->wgts[ls];
              break;
          case tri3:
          case tri6:
              e2   = data->xgps[lr];
              facs = ONE;
              break;
          default:
              dserror("unknown number of gaussian points in ale2_intg");
              break;
      }
     /*---------------------- shape functions and their derivatives */
     ale2_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale2_jaco (deriv,xjm,&det,xyz,iel);
     /*----------------------------- use jacobian determinant or not ---*/
     if(ele->e.ale2->jacobi==1)
       fac = facr * facs * det;
     else
       fac = facr * facs;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale2_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     ale2_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale2_keku(estif,bop,D,fac,nd,numeps);
     /*---------------- hourglass stabalization  stiffness matrix ke ---*/
     if(ele->distyp==quad4 && nir == 1 && nis == 1)
       ale2_hourglass(ele,estif);
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */

/*----------------------------------------------------- local co-system */
if(ele->locsys==locsys_yes)
   locsys_trans(ele,estif,NULL,NULL,NULL);

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_static_ke */
/*----------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element including
additional stiffening

<pre>                                                              ck 01/03
This routine integrates the linear stiffness for the 2d ale element
The stiffness is influenced by the element distortion via minimal Jacobian
determinant

</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE2_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration
\param  quality       INT       (i)   flag quality == 0 : no quality monitoring
                                           quality == 1 : aspect ratio
                                           quality == 2 : corner angle
					   quality == 3 : min det F

\warning There is nothing special to this routine
\return void
\sa calling: ale2_intg(), ale2_funct_deriv(), ale2_jaco(), ale2_bop(),
             ale2_mat_linel(), ale2_keku(), ale2_hourglass(),
	     ale2_min_jaco(), write_element_quality();
             called by: ale2()

*----------------------------------------------------------------------*/
void ale2_static_ke_stiff(ELEMENT     *ele,
                          ALE2_DATA   *data,
                          MATERIAL    *mat,
                          ARRAY       *estif_global,
                          INT          init,
	        	  INT          quality)
{
INT                 i,j;             /* some loopers */
INT                 nir=0,nis=0;     /* num GP in r/s/t direction */
INT                 lr, ls;          /* loopers over GP */
INT                 iel;             /* numnp to this element */
INT                 nd;

const INT           numdf  = 2;
const INT           numeps = 3;

DOUBLE              fac;
DOUBLE              e1,e2=0;         /* GP-coords */
DOUBLE              facr,facs=0;     /* weights at GP */

DOUBLE              min_detF;         /* minimal Jacobian determinant */

static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    xyz_a;    /* actual element coordiantes */
static DOUBLE **xyz;
static ARRAY    fint_a;   /* internal force vector from prestress */
static DOUBLE  *fint;
static DOUBLE **estif;    /* element stiffness matrix ke */

DOUBLE det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_static_ke_stiff");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");
deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");
D         = amdef("D"      ,&D_a   ,6,6              ,"DA");
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");
xyz       = amdef("xyz"    ,&xyz_a ,4    ,numdf      ,"DA");
fint      = amdef("fint"   ,&fint_a,4*numdf,1        ,"DV");
bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");
goto end;
}
/*------------------------------------------- integration parameters ---*/
ale2_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
switch (ele->distyp)
{
    case quad4:
    case quad8:
    case quad9:
        nir     = ele->e.ale2->nGP[0];
        nis     = ele->e.ale2->nGP[1];
        break;
    case tri3:
    case tri6:
        nir = ele->e.ale2->nGP[0];
        nis = 1;
        break;
    default:
        dserror("unknown number of gaussian points in ale2_intg");
        break;
}

iel     = ele->numnp;
nd      = numdf * iel;
/*--------------------------------------- actual element coordinates ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
      xyz[i][j] = ele->node[i]->x[j] + ele->node[i]->sol_increment.a.da[1][j];
   }
}
/*================================================== element quality ===*/
/*------------------------------------------------look for min(det F)---*/
ale2_min_jaco(ele->distyp,xyz,&min_detF);
/*----------------------------------- write element quality meassure ---*/
write_element_quality(ele,quality,xyz,min_detF);
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     /*============================= gaussian point and weight at it ===*/
      switch (ele->distyp)
      {
          case quad4:
          case quad8:
          case quad9:
              e2   = data->xgps[ls];
              facs = data->wgts[ls];
              break;
          case tri3:
          case tri6:
              e2   = data->xgps[lr];
              facs = ONE;
              break;
          default:
              dserror("unknown number of gaussian points in ale2_intg");
              break;
      }
     /*-------------------------- shape functions and their derivatives */
     ale2_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale2_jaco (deriv,xjm,&det,xyz,iel);
     /*------------------------ evaluate factor including stiffening ---*/
     fac = facr * facs * det/min_detF/min_detF;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale2_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     ale2_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale2_keku(estif,bop,D,fac,nd,numeps);
     /*---------------- hourglass stabalization  stiffness matrix ke ---*/
/*     if(nir == 1 && nis == 1)
       ale2_hourglass(ele,estif); */
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_static_ke_stiff */
/*----------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element with linear
prestress

<pre>                                                              ck 01/03
This routine integrates the linear stiffness for the 2d ale element and
calculates element node forces from prestress evaluated from element
distortion

</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE2_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param  init          INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration
\param *rhs           DOUBLE    (o)   element contribution to right hand side
\param  quality       INT       (i)   flag quality == 0 : no quality monitoring
                                           quality == 1 : aspect ratio
                                           quality == 2 : corner angle
					   quality == 3 : min det F

\warning There is nothing special to this routine
\return void
\sa calling: ale2_intg(), ale2_funct_deriv(), ale2_jaco(), ale2_bop(),
             ale2_mat_linel(), ale2_keku(), ale2_hourglass(),
	     ale2_min_jaco(), write_element_quality();
             called by: ale2()

*----------------------------------------------------------------------*/
void ale2_static_ke_prestress(ELEMENT    *ele,
                              ALE2_DATA  *data,
                              MATERIAL   *mat,
                              ARRAY      *estif_global,
                              INT         init,
		              DOUBLE     *rhs,
       			      INT         total_dim,
		              INT         quality)
{
INT                 i,j;              /* some loopers */
INT                 nir=0,nis=0;      /* num GP in r/s/t direction */
INT                 lr, ls;           /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;

const INT           numdf  = 2;
const INT           numeps = 3;

DOUBLE              fac;
DOUBLE              e1,e2=0;          /* GP-coords */
DOUBLE              facr,facs=0;      /* weights at GP */
INT                 lm[8];

DOUBLE              a,b,c;            /* geometric parameters */
DOUBLE              el_area;          /* element area */
DOUBLE              P[2],Q[2];        /* two side mid points */
DOUBLE              mid[2];           /* element middle point */
DOUBLE              square[4][2];     /* nodes of square with same area*/
DOUBLE              min_detF;         /* minimal Jacobian determinant */

static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    xyz_a;    /* actual element coordiantes */
static DOUBLE **xyz;
static ARRAY    fint_a;   /* internal force vector from prestress */
static DOUBLE  *fint;
static DOUBLE **estif;    /* element stiffness matrix ke */

DOUBLE det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_static_ke_prestress");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
funct     = amdef("funct" ,&funct_a,MAXNOD_BRICK1,1 ,"DV");
deriv     = amdef("deriv" ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");
D         = amdef("D"     ,&D_a   ,6,6  	    ,"DA");
xjm       = amdef("xjm"   ,&xjm_a ,numdf,numdf      ,"DA");
xyz       = amdef("xyz"   ,&xyz_a ,4	,numdf      ,"DA");
fint      = amdef("fint"  ,&fint_a,4*numdf,1	    ,"DV");
bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");
goto end;
}
/*------------------------------------------- integration parameters ---*/
ale2_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
switch (ele->distyp)
{
    case quad4:
    case quad8:
    case quad9:
        nir     = ele->e.ale2->nGP[0];
        nis     = ele->e.ale2->nGP[1];
        break;
    case tri3:
    case tri6:
        nir = ele->e.ale2->nGP[0];
        nis = 1;
        break;
    default:
        dserror("unknown number of gaussian points in ale2_intg");
        break;
}

iel     = ele->numnp;
nd      = numdf * iel;
/*--------------------------------------- actual element coordinates ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
      xyz[i][j] = ele->node[i]->x[j] + ele->node[i]->sol_increment.a.da[1][j];
   }
}
/*------------------------------------------------look for min(det F)---*/
ale2_min_jaco(ele->distyp,xyz,&min_detF);
/*==================================================== some geometry ===*/
/*----------------------------------------------------------------------*/
el_area = ale2_el_area(xyz);
/*------------------------------------ evaluate side of equal square ---*/
a = sqrt(el_area);
/*------------------------------------------------------- mid points ---*/
P[0] = 0.5 * (xyz[1][0] + xyz[2][0]);
P[1] = 0.5 * (xyz[1][1] + xyz[2][1]); /* mid of point 1 and 2 */
Q[0] = 0.5 * (xyz[0][0] + xyz[3][0]);
Q[1] = 0.5 * (xyz[0][1] + xyz[3][1]); /* mid of point 0 and 3 */
mid[0] = 0.5 * (P[0] + Q[0]);
mid[1] = 0.5 * (P[1] + Q[1]);         /* element mid point */
/*---------------------------------------------------- shift Q and P ---*/
b = sqrt((Q[0]-P[0])*(Q[0]-P[0])+(Q[1]-P[1])*(Q[1]-P[1])); /* dist P-Q */
Q[0] = mid[0] + (Q[0]-P[0])/b * 0.5*a;
Q[1] = mid[1] + (Q[1]-P[1])/b * 0.5*a;
P[0] = 2.0*mid[0] - Q[0];
P[1] = 2.0*mid[1] - Q[1];
/*------        --------------------------------- replacement square ---*/
b = Q[0]-mid[0];
c = Q[1]-mid[1];
square[0][0] = Q[0] - c;
square[0][1] = Q[1] + b;
square[1][0] = square[0][0] - 2.0 * b;
square[1][1] = square[0][1] - 2.0 * c;
square[2][0] = square[1][0] + 2.0 * c;
square[2][1] = square[1][1] - 2.0 * b;
square[3][0] = square[0][0] + 2.0 * c;
square[3][1] = square[0][1] - 2.0 * b;
/*----------------------------------- write element quality meassure ---*/
write_element_quality(ele,quality,xyz,min_detF);
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {

     /*============================= gaussian point and weight at it ===*/
      switch (ele->distyp)
      {
          case quad4:
          case quad8:
          case quad9:
              e2   = data->xgps[ls];
              facs = data->wgts[ls];
              break;
          case tri3:
          case tri6:
              e2   = data->xgps[lr];
              facs = ONE;
              break;
          default:
              dserror("unknown number of gaussian points in ale2_intg");
              break;
      }
     /*-------------------------- shape functions and their derivatives */
     ale2_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale2_jaco (deriv,xjm,&det,xyz,iel);
     /*----------------------------- use jacobian determinant or not ---*/
     if(ele->e.ale2->jacobi==1)   fac = facr * facs * det;
     else         fac = facr * facs;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale2_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     ale2_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale2_keku(estif,bop,D,fac,nd,numeps);
     /*---------------- hourglass stabalization  stiffness matrix ke ---*/
     if(nir == 1 && nis == 1)
       ale2_hourglass(ele,estif);
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */

/*========================== evaluate internal forces from prestress ===*/
for (i=0; i<nd; i++) fint[i] = 0;

for (i=0; i<iel; i++)
{
   for (j=0; j<iel; j++)
   {
      fint[2*i]   += estif[2*i][2*j]     * (square[j][0]-xyz[j][0]);
      fint[2*i]   += estif[2*i][2*j+1]   * (square[j][1]-xyz[j][1]);
      fint[2*i+1] += estif[2*i+1][2*j]   * (square[j][0]-xyz[j][0]);
      fint[2*i+1] += estif[2*i+1][2*j+1] * (square[j][1]-xyz[j][1]);
   }
   for (j=0; j<numdf; j++)
   {
      lm[i*numdf+j] = ele->node[i]->dof[j];
   }
}
/*------------------------------------ assemble in global rhs-vector ---*/
for (i=0; i<nd; i++)
{
   if (lm[i] >= total_dim) continue;
   rhs[lm[i]] += fint[i];
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_static_ke_prestress */
/*----------------------------------------------------------------------*/




/*!----------------------------------------------------------------------
\brief  integration of stiffness ke for ALE element depending on
displacement of test step

<pre>                                                              ck 01/03
This routine integrates the stiffness for the 2d ale element depending on
the results of a test step [see Chiandussi et al. 2000]

</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE2_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration
\param quality        INT       (i)   flag quality == 0 : no quality monitoring
                                           quality == 1 : aspect ratio
                                           quality == 2 : corner angle
					   quality == 3 : min det F
\param *min_stiff     DOUBLE    (o)   minimal stiffness
\param *max_stiff     DOUBLE    (o)   maximal stiffness
\param *min           DOUBLE    (i)   minimal stiffness of previous step
\param *max           DOUBLE    (i)   maximal stiffness of previous step

\warning There is nothing special to this routine
\return void
\sa calling: ale2_min_jaco(), ale2_el_area(), ale2_intg(),
             ale2_funct_deriv(), ale2_jaco(), ale2_bop(),
	     ale2_mat_linel(), ale2_keku(), ale2_hourglass(),
	     ale2_min_jaco(), write_element_quality();
             called by: ale2()

*----------------------------------------------------------------------*/
void ale2_static_ke_step2(ELEMENT    *ele,
                          ALE2_DATA  *data,
                          MATERIAL   *mat,
                          ARRAY      *estif_global,
                          INT	      init,
                          INT	      quality,
			  DOUBLE     *min_stiff,
			  DOUBLE     *max_stiff,
			  DOUBLE     *min,
			  DOUBLE     *max)
{
INT                 i,j;            /* some loopers */
INT                 nir=0,nis=0;    /* num GP in r/s/t direction */
INT                 lr, ls;         /* loopers over GP */
INT                 iel;            /* numnp to this element */
INT                 nd;

const INT           numdf =2;
const INT           numeps=3;

DOUBLE              fac;
DOUBLE              e1,e2=0;          /* GP-coords */
DOUBLE              facr,facs=0;      /* weights at GP */
DOUBLE              min_detF;         /* minimal Jacobian determinant */
DOUBLE              el_area;          /* elemental area */
DOUBLE              stiff;            /* stiffness factor */
DOUBLE              pv;               /* Possions ratio */
DOUBLE              m,r;

static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    xyz_a;    /* actual element coordiantes */
static DOUBLE **xyz;
static ARRAY    uxyz_a;   /* actual displacements */
static DOUBLE **uxyz;
static ARRAY    strain_a; /* strains */
static DOUBLE  *strain;
static DOUBLE **estif;    /* element stiffness matrix ke */

DOUBLE det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_static_ke_step2");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");
deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");
D         = amdef("D"      ,&D_a   ,6,6              ,"DA");
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");
xyz       = amdef("xyz"    ,&xyz_a ,4    ,numdf      ,"DA");
uxyz      = amdef("uxyz"   ,&uxyz_a ,4   ,numdf      ,"DA");

bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");
strain    = amdef("strain",&strain_a ,numeps,1,"DV");
goto end;
}
/*------------------------------------------- integration parameters ---*/
ale2_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
pv        = mat->m.stvenant->possionratio;
/*------------------------------------------- integration parameters ---*/
switch (ele->distyp)
{
    case quad4:
    case quad8:
    case quad9:
        nir     = ele->e.ale2->nGP[0];
        nis     = ele->e.ale2->nGP[1];
        break;
    case tri3:
    case tri6:
        nir = ele->e.ale2->nGP[0];
        nis = 1;
        break;
    default:
        dserror("unknown number of gaussian points in ale2_intg");
        break;
}

iel     = ele->numnp;
nd      = numdf * iel;
/*--------------------------------------- actual element coordinates ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
      /* real solution so far */
      xyz[i][j]  = ele->node[i]->x[j] + ele->node[i]->sol_increment.a.da[1][j];
      /* displacement of 'trial' step */
      uxyz[i][j] = ele->node[i]->sol_increment.a.da[0][j];
   }
}
/*----------------------------------------------------------------------*/
   ale2_min_jaco(ele->distyp,xyz,&min_detF);
   el_area = ale2_el_area(xyz);
/*----------------------------------- write element quality meassure ---*/
write_element_quality(ele,quality,xyz,min_detF);
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     for (i=0; i<numeps; i++)
     {
        strain[i] = 0.0;
     }
     /*============================= gaussian point and weight at it ===*/
      switch (ele->distyp)
      {
          case quad4:
          case quad8:
          case quad9:
              e2   = data->xgps[ls];
              facs = data->wgts[ls];
              break;
          case tri3:
          case tri6:
              e2   = data->xgps[lr];
              facs = ONE;
              break;
          default:
              dserror("unknown number of gaussian points in ale2_intg");
              break;
      }
     /*-------------------------- shape functions and their derivatives */
     ale2_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale2_jaco(deriv,xjm,&det,xyz,iel);
     /*----------------------------- use jacobian determinant or not ---*/
     if(ele->e.ale2->jacobi==1)
       fac = facr * facs * det;
     else
       fac = facr * facs;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale2_bop(bop,deriv,xjm,det,iel);
     /*-------------------------------------------- evaluate strains ---*/
     for (i=0; i<iel; i++)
     {
        strain[0] += bop[0][2*i]*uxyz[i][0] + bop[0][2*i+1]*uxyz[i][1];
        strain[1] += bop[1][2*i]*uxyz[i][0] + bop[1][2*i+1]*uxyz[i][1];
        strain[2] += bop[2][2*i]*uxyz[i][0] + bop[2][2*i+1]*uxyz[i][1];
     }
     /* projection on mean strains via Mohr's circle */
     m = (strain[0]+strain[1])/2.0;
     r = sqrt( (strain[0]-strain[1])*(strain[0]-strain[1])
                + 4.0*strain[2]*strain[2] );
     strain[0] = m + r;  /* epsilon_11 */
     strain[1] = m - r;  /* epsilon_22 */
     /*------------------------- in the first step where the rhs = 0 ---*/
     if (strain[0] == 0.0 && strain[1] == 0.0)
        stiff = 1.0;
     /*--------------- else for ordinary steps, where strains appear ---*/
     else
     {
     /* equations refere to paper of Chiandussi et al. 2000 */
     /* elemental strain energy density criterion: eq. (15) */
/*        stiff = (1.0-pv)/2.0/(1.0-2.0*pv)/(1+pv) * ( strain[0]*strain[0]
              + strain[1]*strain[1] + 2.0*pv/(1.0-pv)*strain[0]*strain[1]); */
        /* elemental distortion energy density criterion: eq. (17) */
/*        stiff = (strain[0]-strain[1])*(strain[0]-strain[1])/12.0/(1.0 + pv); */
        /* square norm of element principal strain criterion: eq. (13) */
        stiff = (strain[0]*strain[0] + strain[1]*strain[1]) / 2.0;
        *min_stiff = (stiff < *min_stiff) ? stiff:*min_stiff;
        *max_stiff = (stiff > *max_stiff) ? stiff:*max_stiff;
        /* --- scale stiffening factor --- */
        stiff = ((stiff - *min)/(*max - *min) * 99.0) + 1.0;
        /* --- use stiffening factor --- */
     }
     fac = fac*stiff;
     /*------------------------------------------- call material law ---*/
     ale2_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale2_keku(estif,bop,D,fac,nd,numeps);
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_static_ke_step2 */
/*----------------------------------------------------------------------*/




/*!----------------------------------------------------------------------
\brief  linear stiffness ke for spring ALE element

<pre>                                                           ck 06/03
This routine

</pre>
\param *ele           ELEMENT   (i)   my element
\param *estif_global  ARRAY     (o)   the stifness matrix
\param  quality       INT       (i)   flag quality == 0 : no quality monitoring
                                           quality == 1 : aspect ratio
                                           quality == 2 : corner angle
					   quality == 3 : min det F
\param  init          INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void
\sa calling: ale2_min_jaco(), edge_geometry(), write_element_quality()
             called by: ale2_static_ke()

*----------------------------------------------------------------------*/
void ale2_static_ke_spring(ELEMENT   *ele,
                           ARRAY     *estif_global,
			   INT        quality,
			   INT        init)
{
INT                 i, j;
INT                 iel;             /* numnp to this element */
INT                 nd;
INT                 numdf=2;         /* degree of freedom per node */
INT                 node_i, node_j;  /* end nodes of actual spring */

DOUBLE              min_detF;
DOUBLE              length;          /* length of actual edge*/
DOUBLE              sin, cos;        /* direction of actual edge */
DOUBLE              factor;

static DOUBLE **estif;    /* element stiffness matrix ke */

static ARRAY    xyz_a;    /* actual element coordiantes */
static DOUBLE **xyz;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_static_ke_spring");
#endif
/*--------------------------------------------------- initialisation ---*/
if (init==1)
{
   xyz = amdef("xyz", &xyz_a, 4, numdf, "DA");
   ale2_tors_spring_quad4(NULL,NULL,1);
   ale2_tors_spring_tri3(NULL,NULL,1);
   goto end;
}
/*------------------------------------------- zero element stiffness ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*----------------------------------------------------------------------*/
iel     = ele->numnp;
nd      = numdf * iel;
/*--------------------------------------- actual element coordinates ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
         xyz[i][j] = ele->node[i]->x[j] + ele->node[i]->sol_increment.a.da[1][j];
   }
}
/* check 'jacobian determinant' just to determine degenerated elements -*/
ale2_min_jaco(ele->distyp,xyz,&min_detF);
/*----------------------------------- write element quality meassure ---*/
write_element_quality(ele,quality,xyz,min_detF);
/*----------------------- lineal springs from all nodes to all nodes ---*/
/*----------------- loop over all edges and diagonals of the element ---*/
for (node_i=0; node_i<iel; node_i++)
{
   for (node_j=node_i+1; node_j<iel; node_j++)
   {
      edge_geometry(node_i,node_j,xyz,&length,&sin,&cos);
      factor = 1.0 / length;
      /*-------------------------- put values in 'element stiffness' ---*/
      estif[node_i*2][node_i*2]     += cos*cos * factor;
      estif[node_i*2+1][node_i*2+1] += sin*sin * factor;
      estif[node_i*2][node_i*2+1]   += sin*cos * factor;
      estif[node_i*2+1][node_i*2]   += sin*cos * factor;

      estif[node_j*2][node_j*2]     += cos*cos * factor;
      estif[node_j*2+1][node_j*2+1] += sin*sin * factor;
      estif[node_j*2][node_j*2+1]   += sin*cos * factor;
      estif[node_j*2+1][node_j*2]   += sin*cos * factor;

      estif[node_i*2][node_j*2]     -= cos*cos * factor;
      estif[node_i*2+1][node_j*2+1] -= sin*sin * factor;
      estif[node_i*2][node_j*2+1]   -= sin*cos * factor;
      estif[node_i*2+1][node_j*2]   -= sin*cos * factor;

      estif[node_j*2][node_i*2]     -= cos*cos * factor;
      estif[node_j*2+1][node_i*2+1] -= sin*sin * factor;
      estif[node_j*2][node_i*2+1]   -= sin*cos * factor;
      estif[node_j*2+1][node_i*2]   -= sin*cos * factor;
   }
}
/*--------------------------------------- build in torsional springs ---*/
switch (ele->distyp)
{
    case quad4:
       ale2_tors_spring_quad4(estif,xyz,0);
       break;
    case quad8:
    case quad9:
       dserror("ale spring dynamic for this distyp not yet implemented");
       break;
    case tri3:
       ale2_tors_spring_tri3(estif,xyz,0);
       break;
    case tri6:
       dserror("ale spring dynamic for this distyp not yet implemented");
       break;
    default:
       dserror("unknown distyp in ale spring dynamic");
       break;
}
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_static_ke_spring */
/*----------------------------------------------------------------------*/



/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element with
Laplacian smoothing

<pre>                                                              ck 01/03
This routine integrates the linear 'stiffness' for the 2d ale element
stiffness is determined by laplacian smoothing of the velocity (or displacement
increment)
In difference to the mentioned paper the diffusion coefficient does not
depend upon the distance to the moving surface but rather on the shape of
the element which is represented by the minimal Jacobian determinant.

</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE2_DATA (i)   structure containing gaussian point and weight
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration
\param  quality       INT       (i)   flag quality == 0 : no quality monitoring
                                           quality == 1 : aspect ratio
                                           quality == 2 : corner angle
					   quality == 3 : min det F

\warning There is nothing special to this routine
\return void
\sa calling: ale2_intg(), ale2_funct_deriv(), ale2_jaco(), ale2_bop(),
             ale2_mat_linel(), ale2_keku(), ale2_hourglass(),
	     ale2_min_jaco(), write_element_quality();
             called by: ale2()

*----------------------------------------------------------------------*/
void ale2_static_ke_laplace(ELEMENT     *ele,
                            ALE2_DATA   *data,
                            ARRAY       *estif_global,
                            INT          init,
		            INT          quality)
{
INT                 i,j;             /* some loopers */
INT                 nir=0,nis=0;     /* num GP in r/s/t direction */
INT                 lr, ls;          /* loopers over GP */
INT                 iel;             /* numnp to this element */
INT                 nd;

const INT           numdf  = 2;

DOUBLE              fac;
DOUBLE              e1,e2=0;          /* GP-coords */
DOUBLE              facr,facs=0;      /* weights at GP */

DOUBLE              min_detF;         /* minimal Jacobian determinant */

DOUBLE              k_diff;

static ARRAY    D_a;        /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;    /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;    /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;      /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    deriv_xy_a; /* global derivatives */
static DOUBLE **deriv_xy;
static ARRAY    xyz_a;      /* actual element coordiantes */
static DOUBLE **xyz;
static ARRAY    fint_a;     /* internal force vector from prestress */
static DOUBLE  *fint;
static DOUBLE **estif;      /* element stiffness matrix ke */

DOUBLE det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_static_ke_laplace");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");
deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");
D         = amdef("D"      ,&D_a   ,6,6              ,"DA");
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");
xyz       = amdef("xyz"    ,&xyz_a ,4    ,numdf      ,"DA");
fint      = amdef("fint"   ,&fint_a,4*numdf,1        ,"DV");
deriv_xy  = amdef("deriv_xy",&deriv_xy_a ,numdf,(MAXNOD_BRICK1),"DA");
goto end;
}
/*------------------------------------------- integration parameters ---*/
ale2_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
switch (ele->distyp)
{
    case quad4:
    case quad8:
    case quad9:
        nir     = ele->e.ale2->nGP[0];
        nis     = ele->e.ale2->nGP[1];
        break;
    case tri3:
    case tri6:
        nir = ele->e.ale2->nGP[0];
        nis = 1;
        break;
    default:
        dserror("unknown number of gaussian points in ale2_intg");
        break;
}

iel     = ele->numnp;
nd      = numdf * iel;
/*--------------------------------------- actual element coordinates ---*/
for (i=0; i<iel; i++)
{
   for (j=0; j<numdf; j++)
   {
      xyz[i][j] = ele->node[i]->x[j] + ele->node[i]->sol_increment.a.da[1][j];
   }
}
/*================================================== element quality ===*/
/*------------------------------------------------look for min(det F)---*/
ale2_min_jaco(ele->distyp,xyz,&min_detF);
/*----------------------------------- write element quality meassure ---*/
write_element_quality(ele,quality,xyz,min_detF);
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     /*============================= gaussian point and weight at it ===*/
      switch (ele->distyp)
      {
          case quad4:
          case quad8:
          case quad9:
              e2   = data->xgps[ls];
              facs = data->wgts[ls];
              break;
          case tri3:
          case tri6:
              e2   = data->xgps[lr];
              facs = ONE;
              break;
          default:
              dserror("unknown number of gaussian points in ale2_intg");
              break;
      }
     /*-------------------------- shape functions and their derivatives */
     ale2_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale2_jaco (deriv,xjm,&det,xyz,iel);
     /*--------------------------------------------- evaluate factor ---*/
     fac = facr * facs * det;
     /*-------------------------------- calculate global derivatives ---*/
     amzero(&deriv_xy_a);
     ale2_deriv_xy(deriv_xy,deriv,xjm,det,iel);
     /*------------------------- diffusivity depends on displacement ---*/
     k_diff = 1.0/min_detF/min_detF;
     /*------------------------------- sort it into stiffness matrix ---*/
     for (i=0; i<iel; i++)
     {
        for (j=0; j<iel; j++)
	{
	   estif[i*2][j*2]     += ( deriv_xy[0][i] * deriv_xy[0][j]
	                          + deriv_xy[1][i] * deriv_xy[1][j] )*fac*k_diff;
	   estif[i*2+1][j*2+1] += ( deriv_xy[0][i] * deriv_xy[0][j]
	                          + deriv_xy[1][i] * deriv_xy[1][j] )*fac*k_diff;
	}
     }
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_static_ke_laplace */
/*----------------------------------------------------------------------*/


/*! @} (documentation module close)*/
#endif
