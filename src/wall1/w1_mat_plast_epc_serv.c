/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1iwadi' which initializes d-components
       d[0][3], d[1][3], d[2][3], d[3][3] on elastic
       contains the routine 'w1yicsr' which calculates the yield
       criterion for planes strain - Drucker Prager
       contains the routine 'w1yiccap' which calculates the yield
       criterion in the cap-region of the cap-model
       contains the routine 'w1preds' which calculates the gradients of
       the elastic predictor deviatoric/hydrostatic stresses
       contains the routine 'w1pres' which calculates the deviatoric/
       hydrostatic stress components of the elastic predictor

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------|
|      initialize d-components d[0][3], d[1][3], d[2][3], d[3][3]        |
|      on elastic                                                        |
|-----------------------------------------------------------------------*/
void w1iwadi (DOUBLE ym,
              DOUBLE pv,
              DOUBLE *di)
{
/*----------------------------------------------------------------------*/
    static DOUBLE a1, b1, c1;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_enter("w1iwadi");
    #endif
/*----------------------------------------------------------------------*/
    c1 = ym / (pv + 1.);
    b1 = c1 * pv / (1. - pv * 2.);
    a1 = b1 + c1;
/*----- initial components of the constitutive tensor set on elastic ---*/
    di[0] = b1;
    di[1] = b1;
    di[2] = 0.;
    di[3] = a1;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1iwadi */
/*-----------------------------------------------------------------------|
|      topic : wall1 - yield criterion for plane strain                  |
|              yield - criterion:  drucker - prager                      |
|-----------------------------------------------------------------------*/
void w1yicsr (DOUBLE dev,      /* norm der deviatorischen spannungen    */
              DOUBLE hyd,      /* 1. invariante                         */
              DOUBLE sigym,    /* modifizierte fliesspannung            */
              DOUBLE alpha,    /* neigungswinkel im invariantenraum     */
              DOUBLE *ft)      /* fliessbedingung                       */
{
/*----------------------------------------------------------------------*/
    static DOUBLE ro23;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_enter("w1yicsr");
    #endif
/*------------------------------------------------------- initialize ---*/
    *ft = 0.;
    ro23 = sqrt(2./3.);
/*---------------------------------- yield condition: drucker-prager ---*/
    *ft  = dev + alpha * hyd - ro23 * sigym;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1yicsr */
/*-----------------------------------------------------------------------|
|      topic : wall1 - yield criterion (cap - region)                    |
|              yield - criterion:  cap - model                           |
|-----------------------------------------------------------------------*/
void w1yiccap (DOUBLE dev,   /* norm of deviatoric predictor stresses   */
               DOUBLE hyd,   /* 1. invariant  of the predictor stresses */
               DOUBLE *sigym,/* yield stresses                          */
               DOUBLE alpha, /* factor for the first invariants         */
               DOUBLE *ft)   /* yield function                          */
{
/*----------------------------------------------------------------------*/
    static DOUBLE dum;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_enter("w1yiccap");
    #endif
/*------------------------------------------------------- initialize ---*/
    *ft = 0.;
/*---------------------------------- yield condition: drucker-prager ---*/
    dum = hyd-alpha*sigym[2];
    *ft = sqrt(dev*dev + (dum*dum)/9.);
    *ft -= sigym[3];
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return;
} /* end of w1yiccap */
/*-----------------------------------------------------------------------|
|    topic: calculate the gradients of the elastic predictor             |
|           deviatoric/hydrostatic stresses                              |
|-----------------------------------------------------------------------*/
void w1preds(DOUBLE *sigma, /* elastic predictor projected              */
             DOUBLE *alpha, /* neigungswinkel im invariantenraum        */
             DOUBLE *devsig,/* deviatoric stresses                      */
             DOUBLE *dev,   /* norm der deviatorischen spannungen       */
             DOUBLE *dn,    /* gradient in deviatoric direction         */
             DOUBLE *dcom,  /* gradient in hdrostatic direction         */
             DOUBLE *grad)  /* total gradient                           */
{
/*----------------------------------------------------------------------*/
    static INT i;
    static DOUBLE half, ro23, q13, xsi1, xsi2, xsi3, xsi4;
    static DOUBLE dsig[4], fac[4];
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_enter("w1preds");
    #endif
/*----------------------------------------------------------------------*/
      half = 1./2.;
      ro23 = sqrt(2./3.);
      q13  = 1./3.;
/*-------------------------------------------- initialized variables ---*/
    for (i = 0; i < 3; i++) fac[i]  = 1.;
    fac[3]  = 0.;

    xsi1 = sigma[0];
    xsi2 = sigma[1];
    xsi3 = sigma[3];
    xsi4 = sigma[2];

    for (i = 0; i < 4; i++) dsig[i]  = devsig[i];
/*------------------------------------------------ stress components ---*/
      if((sigma[0]==sigma[1])&&(sigma[0]==sigma[3]))
      {
        if(fabs(sigma[3])-1.0E-4<=0.)
        {
          xsi3 = 1.0E-4;
        }
        else
        {
          xsi3 = 1.001 * sigma[3];
        }

       *dev  =    q13 * (xsi1*(2. *xsi1-xsi2-xsi3) +
                         xsi2*(2. *xsi2-xsi1-xsi3) +
                         xsi3*(2. *xsi3-xsi1-xsi2))
                  + 2. *(xsi4*xsi4);
       *dev = sqrt(*dev);

        dsig[0] = ( 2. *xsi1 - xsi2 - xsi3 ) / 3.;
        dsig[1] = (-xsi1 + 2. *xsi2 - xsi3 ) / 3.;
        dsig[2] = (-xsi1  -xsi2 + 2. *xsi3 ) / 3.;
        dsig[3] = xsi4;
      }
/*--------------------------------------- components ot the gradient ---*/
    for (i = 0; i < 4; i++)
    {
      dn[i]   = dsig[i] / *dev  ;
      dcom[i] = fac[i] * *alpha ;
      grad[i] = dn[i] + dcom[i];
    }
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1preds */
/*-----------------------------------------------------------------------|
|    topic: calculate the deviatoric/hydrostatic stress components       |
|           of the elastic predictor                                     |
|-----------------------------------------------------------------------*/
void w1pres (DOUBLE *sigma ,    /*  elastic predictor projected         */
             DOUBLE *devsig,    /*  deviatoric stresses                 */
             DOUBLE *sm    ,    /*  hydrrostatic stresses               */
             DOUBLE *dev   ,    /*  norm der deviatorischen spannungen  */
             DOUBLE *hyd)       /*  1. invariante                       */
{
/*----------------------------------------------------------------------*/
    static DOUBLE half, dsig[4], dum;
    static INT i;
    static DOUBLE q13, ro23, xsi1, xsi2, xsi3, xsi4;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_enter("w1pres");
    #endif
/*----------------------------------------------------------------------*/
    half = .5;
    ro23 = sqrt(.66666666666666663);
    q13 = .33333333333333331;
/*-------------------------------------------- initialized variables ---*/
    *dev = 1e-4;

    xsi1 = sigma[0];
    xsi2 = sigma[1];
    xsi3 = sigma[3];
    xsi4 = sigma[2];
/*---------------------------------- deviatoric invariant / stresses ---*/
    dum  = q13 * (xsi1 * (xsi1 * 2. - xsi2 - xsi3)
               +  xsi2 * (xsi2 * 2. - xsi1 - xsi3)
               +  xsi3 * (xsi3 * 2. - xsi1 - xsi2))
               +  xsi4 *  xsi4 * 2.;
    *dev = pow(dum, half);

    devsig[0] = ( xsi1 * 2. - xsi2 - xsi3) / 3.;
    devsig[1] = (-xsi1 + xsi2 * 2. - xsi3) / 3.;
    devsig[2] = (-xsi1 - xsi2 + xsi3 * 2.) / 3.;
    devsig[3] =   xsi4;

    for (i = 0; i < 3; i++) dsig[i] = devsig[i];

    if (*dev < 1e-4) *dev = 1e-4;
/*--------------------------------- hydrostatic invariant / stresses ---*/
    *hyd = xsi1 + xsi2 + xsi3;

    sm[0] = *hyd / 3.;
    sm[1] = *hyd / 3.;
    sm[2] = *hyd / 3.;
    sm[3] = 0.;
/*----------------------------------------------------------------------*/
    #ifdef DEBUG
    dstrc_exit();
    #endif
/*----------------------------------------------------------------------*/
    return ;
} /* end of w1pres */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
