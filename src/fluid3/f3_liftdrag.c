/*!---------------------------------------------------------------------
  \file
  \brief

  <pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
  </pre>

---------------------------------------------------------------------*/

/*!
  \addtogroup FLUID3
  *//*! @{ (documentation module open)*/

#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3_prototypes.h"
/*!---------------------------------------------------------------------
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
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*!---------------------------------------------------------------------
  \brief lift, drag and moment calculation for 3d fluid

  <pre>                                                           mn 03/04
  In this function the lift and drag forces and moments are calcutaed for
  one 3d fluid element.

  </pre>

  \param  *ele        ELMENT      (i/o)  actual element
  \param  *container       CONTAINER    (i)    container

  \return void

  ----------------------------------------------------------------------*/
void f3_liftdrag(
    ELEMENT       *ele,
    CONTAINER     *container
    )
{


  INT              i,j,cc;              /* some loopers and counters */
  INT              nir,nis,nit;         /* num GP in r/s/t direction */
  INT              iel;                 /* numnp to this element */
  static DOUBLE    funct[MAXNOD_F3];
  static DOUBLE    deriv[3][MAXNOD_F3];
  static DOUBLE    xjm[3][3];
  DOUBLE           xyze[3][MAXNOD_F3];   /* element-node coordinates */
  DOUBLE           vn[3];
  GSURF	          *actgsurf;
  DOUBLE           center[3];	         /* center of body in fluid */
  DOUBLE           nsigma[12][MAXNOD_F3]; /* nodal fluid stresses */
  DIS_TYP          typ;		         /* element type */
  DOUBLE           sigmaint[12];          /* stresses at int point */
  DOUBLE           forceline[6];         /* traction vector */
  DOUBLE           xyz[3];
  INT              ld_id;
  INT              ngsurf;
  INT              ngnode;
  INT              surf;
  INT              ngr, ngs, ngt;
  DOUBLE           xgp[3], ygp[3], zgp[3];
  DOUBLE           wgx[3], wgy[3], wgz[3];
  INT  gpr, gps, gpt;
  DOUBLE e1, e2, e3, facr, facs, fact, fac;
  DOUBLE det0;
  FLUID_DATA     *data;
  FLUID_DYNAMIC  *fdyn;

  static INT shnod[6][8] = {
    {0, 1, 2, 3,  8,  9, 10, 11},  /* 1st surf nodes */
    {4, 5, 1, 0, 12, 17,  8, 16},  /* 2nd surf nodes */
    {5, 6, 2, 1, 13, 18,  9, 17},
    {3, 2, 6, 7, 10, 18, 14, 19},
    {7, 4, 0, 3, 15, 16, 11, 19},
    {7, 6, 5, 4, 14, 13, 12, 15}
  };
  INT *shn;


#ifdef DEBUG
  dstrc_enter("f3_liftdrag");
#endif


  /* integration parameters */
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;
  nir     = ele->e.f3->nGP[0];
  nis     = ele->e.f3->nGP[1];
  nit     = ele->e.f3->nGP[2];
  iel     = ele->numnp;
  typ     = ele->distyp;
  ngsurf  = ele->g.gvol->ngsurf;


  /* setup individual element arrays */
  cc=0;
  for (i=0;i<iel;i++)
    for (j=0;j<3;j++)
      xyze[j][i] = ele->node[i]->x[j];
  if(iel==20)
  {
    xyze[0][12] = ele->node[16]->x[0];
    xyze[1][12] = ele->node[16]->x[1];
    xyze[2][12] = ele->node[16]->x[2];
    xyze[0][13] = ele->node[17]->x[0];
    xyze[1][13] = ele->node[17]->x[1];
    xyze[2][13] = ele->node[17]->x[2];
    xyze[0][14] = ele->node[18]->x[0];
    xyze[1][14] = ele->node[18]->x[1];
    xyze[2][14] = ele->node[18]->x[2];
    xyze[0][15] = ele->node[19]->x[0];
    xyze[1][15] = ele->node[19]->x[1];
    xyze[2][15] = ele->node[19]->x[2];
    xyze[0][16] = ele->node[12]->x[0];
    xyze[1][16] = ele->node[12]->x[1];
    xyze[2][16] = ele->node[12]->x[2];
    xyze[0][17] = ele->node[13]->x[0];
    xyze[1][17] = ele->node[13]->x[1];
    xyze[2][17] = ele->node[13]->x[2];
    xyze[0][18] = ele->node[14]->x[0];
    xyze[1][18] = ele->node[14]->x[1];
    xyze[2][18] = ele->node[14]->x[2];
    xyze[0][19] = ele->node[15]->x[0];
    xyze[1][19] = ele->node[15]->x[1];
    xyze[2][19] = ele->node[15]->x[2];
  }


  /* read nodal stresses */
  for (j=0;j<iel;j++)
  {
    for (i=0; i<12; i++)
      nsigma[i][j] = ele->e.f3->stress_ND.a.da[j][i];
  }



  /* loop all gsurfs of the element */
  for (surf=0; surf<ngsurf; surf++)
  {
    actgsurf = ele->g.gvol->gsurf[surf];

    /* confirm that the gsurf is the right one that lies on the body */
    if (actgsurf->dsurf==NULL) continue;
    if (actgsurf->dsurf->liftdrag == NULL) continue;


    /* get center of area exposed to lift & drag */
    center[0] = actgsurf->dsurf->liftdrag->ld_center[0];
    center[1] = actgsurf->dsurf->liftdrag->ld_center[1];
    center[2] = actgsurf->dsurf->liftdrag->ld_center[2];

    ld_id     = actgsurf->dsurf->liftdrag->liftdrag;


    /* check number of nodes on surf */
    ngnode = actgsurf->ngnode;

    /* make coordinates and weights of integration points */
    ngr = nir;
    ngs = nis;
    ngt = nit;
    for (i=0; i<ngr; i++)
    {
      xgp[i] = data->qxg[i][nir-1];
      wgx[i] = data->qwgt[i][nir-1];
    }
    for (i=0; i<ngs; i++)
    {
      ygp[i] = data->qxg[i][nis-1];
      wgy[i] = data->qwgt[i][nis-1];
    }
    for (i=0; i<ngt; i++)
    {
      zgp[i] = data->qxg[i][nit-1];
      wgz[i] = data->qwgt[i][nit-1];
    }
    /*------------ get surface-node topology, det. integration order ---*/
    shn=shnod[surf];
    switch (surf)
    {
      case 0:
        ngt = 1;
        zgp[0] = -1.;
        wgz[0] =  1.;
        break;
      case 1:
        ngr = 1;
        xgp[0] =  1.;
        wgx[0] =  1.;
        break;
      case 2:
        ngs = 1;
        ygp[0] =  1.;
        wgy[0] =  1.;
        break;
      case 3:
        ngr = 1;
        xgp[0] = -1.;
        wgx[0] =  1.;
        break;
      case 4:
        ngs = 1;
        ygp[0] = -1.;
        wgy[0] =  1.;
        break;
      case 5:
        ngt = 1;
        zgp[0] =  1.;
        wgz[0] =  1.;
        break;
    }



    /* start loop over integration points */
    for (gpr=0; gpr<ngr; gpr++) {
      e1   = xgp[gpr];
      facr = wgx[gpr];

      for (gps=0; gps<ngs; gps++) {
        e2   = ygp[gps];
        facs = wgy[gps];

        for (gpt=0; gpt<ngt; gpt++) {
          e3   = zgp[gpt];
          fact = wgz[gpt];


          /* shape functions and derivatives at this point */
          f3_hex2(funct,deriv,e1,e2,e3,typ);
          f3_jaco2(funct,deriv,xjm,&det0,ele,iel);

          fac = facr * facs * fact;


          /* compute normal vector for the surface */
          switch (surf)
          {
            case 0:
            case 5:
              vn[0] = xjm[0][1]*xjm[1][2] - xjm[1][1]*xjm[0][2];
              vn[1] = xjm[0][2]*xjm[1][0] - xjm[1][2]*xjm[0][0];
              vn[2] = xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1];
              break;
            case 1:
            case 3:
              vn[0] = xjm[1][1]*xjm[2][2] - xjm[2][1]*xjm[1][2];
              vn[1] = xjm[1][2]*xjm[2][0] - xjm[2][2]*xjm[1][0];
              vn[2] = xjm[1][0]*xjm[2][1] - xjm[2][0]*xjm[1][1];
              break;
            case 4:
            case 2:
              vn[0] = xjm[0][1]*xjm[2][2] - xjm[2][1]*xjm[0][2];
              vn[1] = xjm[0][2]*xjm[2][0] - xjm[2][2]*xjm[0][0];
              vn[2] = xjm[0][0]*xjm[2][1] - xjm[2][0]*xjm[0][1];
              break;
            default:
              dserror("Unknown type of integration");
              break;
          }


          /* interpolate stress tensor to Gauss point */
          for (i=0;i<12;i++) sigmaint[i] = ZERO;

          for (i=0;i<12;i++)
            for (j=0;j<ngnode;j++)
              sigmaint[i] += funct[shn[j]] * nsigma[i][shn[j]];


          /* compute stress vector at gauss point *
             force = sigma * n  (Cauchy's law) */
          /* viscous */
          forceline[0] = sigmaint[0]*vn[0] + sigmaint[3]*vn[1] + sigmaint[5]*vn[2];
          forceline[1] = sigmaint[3]*vn[0] + sigmaint[1]*vn[1] + sigmaint[4]*vn[2];
          forceline[2] = sigmaint[5]*vn[0] + sigmaint[4]*vn[1] + sigmaint[2]*vn[2];

          /* only pressure */
          forceline[3] =  sigmaint[6]*vn[0] +  sigmaint[9]*vn[1] + sigmaint[11]*vn[2];
          forceline[4] =  sigmaint[9]*vn[0] +  sigmaint[7]*vn[1] + sigmaint[10]*vn[2];
          forceline[5] = sigmaint[11]*vn[0] + sigmaint[10]*vn[1] +  sigmaint[8]*vn[2];


          /* compute vector from center to Gauss point */
          for (j=0; j<3; j++)
          {
            xyz[j] = -center[j];
            for (i=0; i<ngnode; i++)
            {
              xyz[j] += funct[shn[i]] * xyze[j][shn[i]];
            }
          }


          /* viscous */
          /* lift and drag forces in global coordinate system */
          container->liftdrag[ld_id*6+0] -= forceline[0] * fac;
          container->liftdrag[ld_id*6+1] += forceline[1] * fac;
          container->liftdrag[ld_id*6+2] += forceline[2] * fac;

          /* add momentums */
          container->liftdrag[ld_id*6+3] +=
            ( forceline[1]*xyz[2] - forceline[2]*xyz[1] ) * fac;

          container->liftdrag[ld_id*6+4] +=
            - ( forceline[0]*xyz[2] - forceline[2]*xyz[0] ) * fac;

          container->liftdrag[ld_id*6+5] +=
            ( forceline[0]*xyz[1] - forceline[1]*xyz[0] ) * fac;


          /* only pressure */
          /* lift and drag forces in global coordinate system */
          container->liftdrag[(FLUID_NUM_LD+1+ld_id)*6+0] -= forceline[3] * fac;
          container->liftdrag[(FLUID_NUM_LD+1+ld_id)*6+1] += forceline[4] * fac;
          container->liftdrag[(FLUID_NUM_LD+1+ld_id)*6+2] += forceline[5] * fac;

          /* add momentums */
          container->liftdrag[(FLUID_NUM_LD+1+ld_id)*6+3] +=
            ( forceline[4]*xyz[2] - forceline[5]*xyz[1] ) * fac;

          container->liftdrag[(FLUID_NUM_LD+1+ld_id)*6+4] +=
            - ( forceline[3]*xyz[2] - forceline[5]*xyz[0] ) * fac;

          container->liftdrag[(FLUID_NUM_LD+1+ld_id)*6+5] +=
            ( forceline[3]*xyz[1] - forceline[4]*xyz[0] ) * fac;

        }
      }
    }/* end loop gpr over integration points */

  }/* end loop surf over surfs */


  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_liftdrag */



static DOUBLE Q12 = ONE/TWO;
static DOUBLE Q14 = ONE/FOUR;
static DOUBLE Q18 = ONE/EIGHT;
/*!---------------------------------------------------------------------
  \brief shape functions and their natural derivatives for hexaeder

  <pre>                                                         mn 03/04

  In this routine the shape functions and their natural first
  derivatives with respect to r/s/t are evaluated for
  H E X A H E D E R


  </pre>
  \param   funct[]     DOUBLE   (o)    shape functions
  \param   deriv[][]   DOUBLE   (o)    1st natural deriv. of shape funct.
  \param   r 	       DOUBLE   (i)    coordinate
  \param   s 	       DOUBLE   (i)    coordinate
  \param   t 	       DOUBLE   (i)    coordinate
  \param   typ 	       DIS_TYP  (i)    element type

  \return void
  \warning shape functions for hex20/hex27/tet10 not implemented yet!!!

  ------------------------------------------------------------------------*/
void f3_hex2(
    DOUBLE      funct[MAXNOD_F3],
    DOUBLE      deriv[3][MAXNOD_F3],
    DOUBLE      r,
    DOUBLE      s,
    DOUBLE      t,
    DIS_TYP     typ
    )
{
  DOUBLE rp,rm,sp,sm,tp,tm;
  DOUBLE rrm,ssm,ttm;
  DOUBLE drm1,dr00,drp1,dsm1,ds00,dsp1,dtm1,dt00,dtp1;
  DOUBLE rm1,r00,rp1,sm1,s00,sp1,tm1,t00,tp1;

#ifdef DEBUG
  dstrc_enter("f3_hex2");
#endif


  /* selection of polynomial interpolation */
  switch (typ)
  {

    case hex8: /* LINEAR shape functions and their natural derivatives */
      rp=ONE+r;
      rm=ONE-r;
      sp=ONE+s;
      sm=ONE-s;
      tp=ONE+t;
      tm=ONE-t;

      funct[0]=Q18*rp*sm*tm;
      funct[1]=Q18*rp*sp*tm;
      funct[2]=Q18*rm*sp*tm;
      funct[3]=Q18*rm*sm*tm;
      funct[4]=Q18*rp*sm*tp;
      funct[5]=Q18*rp*sp*tp;
      funct[6]=Q18*rm*sp*tp;
      funct[7]=Q18*rm*sm*tp;

      /* first derivative */
      deriv[0][0]= Q18*sm*tm  ;
      deriv[0][1]= Q18*sp*tm  ;
      deriv[0][2]=-deriv[0][1];
      deriv[0][3]=-deriv[0][0];
      deriv[0][4]= Q18*sm*tp  ;
      deriv[0][5]= Q18*sp*tp  ;
      deriv[0][6]=-deriv[0][5];
      deriv[0][7]=-deriv[0][4];

      deriv[1][0]=-Q18*tm*rp  ;
      deriv[1][1]=-deriv[1][0];
      deriv[1][2]= Q18*tm*rm  ;
      deriv[1][3]=-deriv[1][2];
      deriv[1][4]=-Q18*tp*rp  ;
      deriv[1][5]=-deriv[1][4];
      deriv[1][6]= Q18*tp*rm  ;
      deriv[1][7]=-deriv[1][6];

      deriv[2][0]=-Q18*rp*sm  ;
      deriv[2][1]=-Q18*rp*sp  ;
      deriv[2][2]=-Q18*rm*sp  ;
      deriv[2][3]=-Q18*rm*sm  ;
      deriv[2][4]=-deriv[2][0];
      deriv[2][5]=-deriv[2][1];
      deriv[2][6]=-deriv[2][2];
      deriv[2][7]=-deriv[2][3];
      break;


    case hex20: /* QUADRATIC shape functions and their natural derivatives
                   without central nodes                      ----*/

      dserror("shape functions for hex20 not implemented yet - see f3_calfuncderiv!\n");
      /*--------------------------------------------------- form basic values */
      rp=ONE+r;
      rm=ONE-r;
      sp=ONE+s;
      sm=ONE-s;
      tp=ONE+t;
      tm=ONE-t;
      rrm=ONE-r*r;
      ssm=ONE-s*s;
      ttm=ONE-t*t;

      funct[0] =Q18*rp*sm*tm*(rp+sm+tm-FIVE);
      funct[1] =Q18*rp*sp*tm*(rp+sp+tm-FIVE);
      funct[2] =Q18*rm*sp*tm*(rm+sp+tm-FIVE);
      funct[3] =Q18*rm*sm*tm*(rm+sm+tm-FIVE);
      funct[4] =Q18*rp*sm*tp*(rp+sm+tp-FIVE);
      funct[5] =Q18*rp*sp*tp*(rp+sp+tp-FIVE);
      funct[6] =Q18*rm*sp*tp*(rm+sp+tp-FIVE);
      funct[7] =Q18*rm*sm*tp*(rm+sm+tp-FIVE);
      funct[8] =Q14*rp*ssm*tm;
      funct[9] =Q14*rrm*sp*tm;
      funct[10]=Q14*rm*ssm*tm;
      funct[11]=Q14*rrm*sm*tm;
      funct[12]=Q14*rp*ssm*tp;
      funct[13]=Q14*rrm*sp*tp;
      funct[14]=Q14*rm*ssm*tp;
      funct[15]=Q14*rrm*sm*tp;
      funct[16]=Q14*rp*sm*ttm;
      funct[17]=Q14*rp*sp*ttm;
      funct[18]=Q14*rm*sp*ttm;
      funct[19]=Q14*rm*sm*ttm;

      /* first derivative evaluation */
      deriv[0][0] = Q18*sm*tm*(TWO*rp+sm+tm-FIVE);
      deriv[1][0] =-Q18*tm*rp*(TWO*sm+tm+rp-FIVE);
      deriv[2][0] =-Q18*rp*sm*(TWO*tm+rp+sm-FIVE);

      deriv[0][1] = Q18*sp*tm*(TWO*rp+sp+tm-FIVE);
      deriv[1][1] = Q18*tm*rp*(TWO*sp+tm+rp-FIVE);
      deriv[2][1] =-Q18*rp*sp*(TWO*tm+rp+sp-FIVE);

      deriv[0][2] =-Q18*sp*tm*(TWO*rm+sp+tm-FIVE);
      deriv[1][2] = Q18*tm*rm*(TWO*sp+tm+rm-FIVE);
      deriv[2][2] =-Q18*rm*sp*(TWO*tm+rm+sp-FIVE);

      deriv[0][3] =-Q18*sm*tm*(TWO*rm+sm+tm-FIVE);
      deriv[1][3] =-Q18*tm*rm*(TWO*sm+tm+rm-FIVE);
      deriv[2][3] =-Q18*rm*sm*(TWO*tm+rm+sm-FIVE);

      deriv[0][4] = Q18*sm*tp*(TWO*rp+sm+tp-FIVE);
      deriv[1][4] =-Q18*tp*rp*(TWO*sm+tp+rp-FIVE);
      deriv[2][4] = Q18*rp*sm*(TWO*tp+rp+sm-FIVE);

      deriv[0][5] = Q18*sp*tp*(TWO*rp+sp+tp-FIVE);
      deriv[1][5] = Q18*tp*rp*(TWO*sp+tp+rp-FIVE);
      deriv[2][5] = Q18*rp*sp*(TWO*tp+rp+sp-FIVE);

      deriv[0][6] =-Q18*sp*tp*(TWO*rm+sp+tp-FIVE);
      deriv[1][6] = Q18*tp*rm*(TWO*sp+tp+rm-FIVE);
      deriv[2][6] = Q18*rm*sp*(TWO*tp+rm+sp-FIVE);

      deriv[0][7] =-Q18*sm*tp*(TWO*rm+sm+tp-FIVE);
      deriv[1][7] =-Q18*tp*rm*(TWO*sm+tp+rm-FIVE);
      deriv[2][7] = Q18*rm*sm*(TWO*tp+rm+sm-FIVE);

      deriv[0][8] = Q14*ssm*tm;
      deriv[1][8] =-Q12*s*tm*rp;
      deriv[2][8] =-Q14*ssm*rp;

      deriv[0][9] =-Q12*r*sp*tm;
      deriv[1][9] = Q14*rrm*tm;
      deriv[2][9] =-Q14*rrm*sp;

      deriv[0][10]=-deriv[0][8];
      deriv[1][10]=-Q12*s*tm*rm;
      deriv[2][10]=-Q14*ssm*rm;

      deriv[0][11]=-Q12*r*sm*tm;
      deriv[1][11]=-deriv[1][9];
      deriv[2][11]=-Q14*rrm*sm;

      deriv[0][12]= Q14*ssm*tp;
      deriv[1][12]=-Q12*s*tp*rp;
      deriv[2][12]=-deriv[2][8];

      deriv[0][13]=-Q12*r*sp*tp;
      deriv[1][13]= Q14*rrm*tp;
      deriv[2][13]=-deriv[2][8];

      deriv[0][14]=-deriv[0][12];
      deriv[1][14]=-Q12*s*tp*rm;
      deriv[2][14]=-deriv[2][10];

      deriv[0][15]=-Q12*r*sm*tp;
      deriv[1][15]=-deriv[1][13];
      deriv[2][15]=-deriv[2][11];

      deriv[0][16]= Q14*sm*ttm;
      deriv[1][16]=-Q14*ttm*rp;
      deriv[2][16]=-Q12*t*rp*sm;

      deriv[0][17]= Q14*sp*ttm;
      deriv[1][17]=-deriv[1][16];
      deriv[2][17]=-Q12*t*rp*sp;

      deriv[0][18]=-deriv[0][17];
      deriv[1][18]= Q14*ttm*rm;
      deriv[2][18]=-Q12*t*rm*sp;

      deriv[0][19]=-deriv[0][16];
      deriv[1][19]=-deriv[1][18];
      deriv[2][19]=-Q12*t*rm*sm;
      break;


    case hex27: /* QUADRATIC shape functions and their natural derivatives
                   with central nodes                         ----*/

      dserror("shape functions for hex20 not implemented yet - see f3_calfuncderiv!\n");
      /*--------------------------------------------------- form basic values */
      rm1=Q12*r*(r - ONE);
      r00=(ONE - r*r);
      rp1=Q12*r*(r + ONE);
      sm1=Q12*s*(s - ONE);
      s00=(ONE - s*s);
      sp1=Q12*s*(s + ONE);
      tm1=Q12*t*(t - ONE);
      t00=(ONE - t*t);
      tp1=Q12*t*(t + ONE);

      funct[0 ]= rp1 * sm1 * tm1;
      funct[1 ]= rp1 * sp1 * tm1;
      funct[2 ]= rm1 * sp1 * tm1;
      funct[3 ]= rm1 * sm1 * tm1;
      funct[4 ]= rp1 * sm1 * tp1;
      funct[5 ]= rp1 * sp1 * tp1;
      funct[6 ]= rm1 * sp1 * tp1;
      funct[7 ]= rm1 * sm1 * tp1;
      funct[8 ]= rp1 * s00 * tm1;
      funct[9 ]= r00 * sp1 * tm1;
      funct[10]= rm1 * s00 * tm1;
      funct[11]= r00 * sm1 * tm1;
      funct[12]= rp1 * s00 * tp1;
      funct[13]= r00 * sp1 * tp1;
      funct[14]= rm1 * s00 * tp1;
      funct[15]= r00 * sm1 * tp1;
      funct[16]= rp1 * sm1 * t00;
      funct[17]= rp1 * sp1 * t00;
      funct[18]= rm1 * sp1 * t00;
      funct[19]= rm1 * sm1 * t00;
      funct[20]= rp1 * s00 * t00;
      funct[21]= r00 * sp1 * t00;
      funct[22]= rm1 * s00 * t00;
      funct[23]= r00 * sm1 * t00;
      funct[24]= r00 * s00 * tp1;
      funct[25]= r00 * s00 * tm1;
      funct[26]= r00 * s00 * t00;

      /* first derivative evaluation */
      drm1 = r - Q12;
      dr00 = -TWO * r;
      drp1 = r + Q12;
      dsm1 = s - Q12;
      ds00 = -TWO * s;
      dsp1 = s + Q12;
      dtm1 = t - Q12;
      dt00 = -TWO * t;
      dtp1 = t + Q12;

      deriv[0][0 ]= drp1 * sm1 * tm1;
      deriv[0][1 ]= drp1 * sp1 * tm1;
      deriv[0][2 ]= drm1 * sp1 * tm1;
      deriv[0][3 ]= drm1 * sm1 * tm1;
      deriv[0][4 ]= drp1 * sm1 * tp1;
      deriv[0][5 ]= drp1 * sp1 * tp1;
      deriv[0][6 ]= drm1 * sp1 * tp1;
      deriv[0][7 ]= drm1 * sm1 * tp1;
      deriv[0][8 ]= drp1 * s00 * tm1;
      deriv[0][9 ]= dr00 * sp1 * tm1;
      deriv[0][10]= drm1 * s00 * tm1;
      deriv[0][11]= dr00 * sm1 * tm1;
      deriv[0][12]= drp1 * s00 * tp1;
      deriv[0][13]= dr00 * sp1 * tp1;
      deriv[0][14]= drm1 * s00 * tp1;
      deriv[0][15]= dr00 * sm1 * tp1;
      deriv[0][16]= drp1 * sm1 * t00;
      deriv[0][17]= drp1 * sp1 * t00;
      deriv[0][18]= drm1 * sp1 * t00;
      deriv[0][19]= drm1 * sm1 * t00;
      deriv[0][20]= drp1 * s00 * t00;
      deriv[0][21]= dr00 * sp1 * t00;
      deriv[0][22]= drm1 * s00 * t00;
      deriv[0][23]= dr00 * sm1 * t00;
      deriv[0][24]= dr00 * s00 * tp1;
      deriv[0][25]= dr00 * s00 * tm1;
      deriv[0][26]= dr00 * s00 * t00;

      deriv[1][0 ]= rp1 * dsm1 * tm1;
      deriv[1][1 ]= rp1 * dsp1 * tm1;
      deriv[1][2 ]= rm1 * dsp1 * tm1;
      deriv[1][3 ]= rm1 * dsm1 * tm1;
      deriv[1][4 ]= rp1 * dsm1 * tp1;
      deriv[1][5 ]= rp1 * dsp1 * tp1;
      deriv[1][6 ]= rm1 * dsp1 * tp1;
      deriv[1][7 ]= rm1 * dsm1 * tp1;
      deriv[1][8 ]= rp1 * ds00 * tm1;
      deriv[1][9 ]= r00 * dsp1 * tm1;
      deriv[1][10]= rm1 * ds00 * tm1;
      deriv[1][11]= r00 * dsm1 * tm1;
      deriv[1][12]= rp1 * ds00 * tp1;
      deriv[1][13]= r00 * dsp1 * tp1;
      deriv[1][14]= rm1 * ds00 * tp1;
      deriv[1][15]= r00 * dsm1 * tp1;
      deriv[1][16]= rp1 * dsm1 * t00;
      deriv[1][17]= rp1 * dsp1 * t00;
      deriv[1][18]= rm1 * dsp1 * t00;
      deriv[1][19]= rm1 * dsm1 * t00;
      deriv[1][20]= rp1 * ds00 * t00;
      deriv[1][21]= r00 * dsp1 * t00;
      deriv[1][22]= rm1 * ds00 * t00;
      deriv[1][23]= r00 * dsm1 * t00;
      deriv[1][24]= r00 * ds00 * tp1;
      deriv[1][25]= r00 * ds00 * tm1;
      deriv[1][26]= r00 * ds00 * t00;

      deriv[2][0 ]= rp1 * sm1 * dtm1;
      deriv[2][1 ]= rp1 * sp1 * dtm1;
      deriv[2][2 ]= rm1 * sp1 * dtm1;
      deriv[2][3 ]= rm1 * sm1 * dtm1;
      deriv[2][4 ]= rp1 * sm1 * dtp1;
      deriv[2][5 ]= rp1 * sp1 * dtp1;
      deriv[2][6 ]= rm1 * sp1 * dtp1;
      deriv[2][7 ]= rm1 * sm1 * dtp1;
      deriv[2][8 ]= rp1 * s00 * dtm1;
      deriv[2][9 ]= r00 * sp1 * dtm1;
      deriv[2][10]= rm1 * s00 * dtm1;
      deriv[2][11]= r00 * sm1 * dtm1;
      deriv[2][12]= rp1 * s00 * dtp1;
      deriv[2][13]= r00 * sp1 * dtp1;
      deriv[2][14]= rm1 * s00 * dtp1;
      deriv[2][15]= r00 * sm1 * dtp1;
      deriv[2][16]= rp1 * sm1 * dt00;
      deriv[2][17]= rp1 * sp1 * dt00;
      deriv[2][18]= rm1 * sp1 * dt00;
      deriv[2][19]= rm1 * sm1 * dt00;
      deriv[2][20]= rp1 * s00 * dt00;
      deriv[2][21]= r00 * sp1 * dt00;
      deriv[2][22]= rm1 * s00 * dt00;
      deriv[2][23]= r00 * sm1 * dt00;
      deriv[2][24]= r00 * s00 * dtp1;
      deriv[2][25]= r00 * s00 * dtm1;
      deriv[2][26]= r00 * s00 * dt00;
      break;


    default:
      dserror("distyp unknown\n");
  } /* end switch (typ) */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_hex2 */




/*!---------------------------------------------------------------------
  \brief jacobian matrix

  <pre>                                                         mn 03/04

  In this routine the jacobian matrix and its determinant is calculated.

  </pre>
  \param   funct     DOUBLE[]    (i)    natural shape functions
  \param   deriv     DOUBLE[][]  (i)    natural deriv. of shape funcs
  \param   xjm       DOUBLE[][]  (o)    jacobian matrix
  \param  *det       DOUBLE      (o)    determinant of jacobian matrix
  \param  *ele       ELEMENT     (i)    actual element
  \param   iel       INT         (i)    num. of nodes of act. ele

  \return void

  ------------------------------------------------------------------------*/
void f3_jaco2(
    DOUBLE      funct[MAXNOD_F3],
    DOUBLE      deriv[3][MAXNOD_F3],
    DOUBLE      xjm[3][3],
    DOUBLE     *det,
    ELEMENT    *ele,
    INT         iel
    )

{
  INT i,j,l;
  DOUBLE dum;

#ifdef DEBUG
  dstrc_enter("f3_jaco2");
#endif


  /* determine jacobian at point r,s,t */
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
    }
  }


  /* determinant of jacobian */
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


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_jaco2 */



#endif
/*! @} (documentation module close)*/

