/*!----------------------------------------------------------------------
\file
\brief contains the routine 'saxi_eleload' which computes the right hand
       side of the symmetric shell element
       
<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*-----------------------------------------------------------------------*/
#ifdef D_AXISHELL
#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"

/*! 
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  Computation of the rhs of the axialsymmetric shell element 

<pre>                                                              mn 05/03 
This routine computes the right hand side of the axialsymmetric shell element
for loads acting on the actual element. The element loads used here are:

 /-----------------------------------------------------------------------*
 |  From Input-file:                                                     |
 |    px...load parallel to xsi direction (meridian direction)           |
 |    pw...load normal to xsi direction (normal direction)               |
 |    pv...load in vertical direction                                    |
 |    ph...load in horizontal direction                                  |
 |  Computed:                                                            |
 |    qx...resultant load in xsi direction (meridian direction)          |
 |    pw...resultant load normal to xsi direction (normal direction)     |
 *-----------------------------------------------------------------------/

 /----------------- SIGN CONVENTION FOR ELEMENT FORCES ------------------*
 |                                                                       |
 |                                   node i                              |
 |                +pv    +pw           o          |                      |
 |                 ^        \         /           |                      |
 |                 |          \      /            |                      |
 |                 |            \   /             |                      |
 |                 |              \/              |                      |
 |                 |------>+ph   //               |                      |
 |                              //                |                      |
 |                             //                 |                      |
 |                            //                  |                      |
 |                         +px/                   |                      |
 |                           /                    |                      |
 |                          o                     |                      |
 |                       node j                   rotational axis        |     
 |    +y                                                                 |
 |    | global                                                           |
 |    | coordinate                                                       |
 |    | system                                                           |
 |    |----->+x                                                          |
 |                                                                       |
 *-----------------------------------------------------------------------/
                                                                         
</pre>                                                                   
\param *ele           ELEMENT   (i)   my element
\param *data          SAXI_DATA (i)   structure containing gaussian point and weight
\param *loadvec       DOUBLE    (i)   global element load vector
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void                                               
\sa calling:   saxiintg;
    called by: axishell();

*----------------------------------------------------------------------*/
void saxi_eleload(ELEMENT    *ele,    
                  SAXI_DATA  *data,   
                  DOUBLE     *loadvec,
                  INT	      init)   
{
INT     i,lxsi,nixsi;
DOUBLE  xsi,fac;
DOUBLE  ri,rj,dr,dz,dl,cosa,sina,r,px,pw,pv,ph,qx,qw;
DOUBLE  px_ij[2],pw_ij[2],pv_ij[2],ph_ij[2];    /* 0..node i, 1..node j */
DOUBLE  loadvec_7[7];     /*local  loadvector before static condensation*/
DOUBLE  loadvec_local[6];   /*local loadvector after static condensation*/

static ARRAY eload_a; 
static DOUBLE **eload;                    /* static element load vector */
DOUBLE *statcond;

/*----------------------------------------------------------------------*/
/* init phase        (init=1)                                           */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("saxi_eleload");
#endif

if (init==1)
{
  eload   = amdef("eload"  ,&eload_a,MAXDOFPERNODE,MAXNOD_AXISHELL,"DA");
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/

else if (init==-1)
{
amdel(&eload_a);
goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase        (init=0)                                    */
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- initialize eload */
amzero(&eload_a);
/*------------------------------------------- integration parameters ---*/
saxiintg(data);
/*--------------------------------------------------- some element infos*/
ri = ele->node[0]->x[0];
rj = ele->node[1]->x[0];
dr = rj - ri;                                 /* delta r */
dz = ele->node[1]->x[1] - ele->node[0]->x[1]; /* delta z */
dl = sqrt(dr*dr+dz*dz);                       /* delta l */
cosa = dr/dl;
sina = dz/dl;

/* general case for a linear distribution of the loads in the element
 ------------------------------------------------------------------------*/
/* Value of the loads at node i and j */
pv_ij[0] = ele->e.saxi->pv[0];
pv_ij[1] = ele->e.saxi->pv[1];
ph_ij[0] = ele->e.saxi->ph[0];
ph_ij[1] = ele->e.saxi->ph[1];
px_ij[0] = ele->e.saxi->px[0];
px_ij[1] = ele->e.saxi->px[1];
pw_ij[0] = ele->e.saxi->pw[0];
pw_ij[1] = ele->e.saxi->pw[1];

for (i=0; i<7; i++)
{
   loadvec_7[i] = 0.0;
}
nixsi   = 5; /* number of Gauss-points */
/*================================================= integration loop ===*/
for (lxsi=0; lxsi<nixsi; lxsi++)
{
   /*==================== gaussian point, weight and thickness at it ===*/
   xsi = data->xgr[lxsi];
   fac = data->wgt[lxsi];
   r = ele->node[1]->x[0] + cosa*xsi*dl;
   
   pv = pv_ij[0]+xsi/1.0*(pv_ij[1]-pv_ij[0]);
   ph = ph_ij[0]+xsi/1.0*(ph_ij[1]-ph_ij[0]);
   px = px_ij[0]+xsi/1.0*(px_ij[1]-px_ij[0]);
   pw = pw_ij[0]+xsi/1.0*(pw_ij[1]-pw_ij[0]);

   qx = px - ph*cosa + pv*sina;  
   qw = pw + ph*sina - pv*cosa;  

   loadvec_7[0] += dl*r*qx*(-1.0+3.0*xsi-2.0*xsi*xsi)        *fac;
   /*zusaetzlicher term bei loadvec_7[0] ?*/
   loadvec_7[1] += dl*r*qw*(-1.0+3.0*xsi*xsi-2.0*xsi*xsi*xsi)*fac;
   loadvec_7[2] += dl*r*qw*dl*xsi*(-1.0+2.0*xsi-xsi*xsi)     *fac;
   loadvec_7[3] += dl*r*qx*xsi*(1.0-2.0*xsi)                 *fac;
   /*zusaetzlicher term bei loadvec_7[3] ?*/
   loadvec_7[4] += dl*r*qw*xsi*xsi*(-3.0+2.0*xsi)            *fac;
   loadvec_7[5] += dl*r*qw*dl*xsi*xsi*(1.0-xsi)              *fac;
   loadvec_7[6] += dl*r*qx*4.0*xsi*(-1.0+xsi)                *fac;
}
/*========================================== end of integration loop ===*/
ele->e.saxi->Sk *= loadvec_7[6]; /*Sk = S^0_7/k77 (siehe Skript S.11.11)*/

/*--------------------- reduction of the rhs-vector from 7 to 6 entries */
statcond = ele->e.saxi->statcond;
for (i=0; i<6; i++)
{
  /*auf statcond soll schon -k_17/k77 stehen!*/
  loadvec_local[i] = loadvec_7[i] + statcond[i] * loadvec_7[6];/*S_ausFlaechenlasten)*/
  loadvec_local[i] = -loadvec_local[i];/* K*u=(S_ausKnotenkraeften-S_ausFlaechenlasten)*/
}

/*------------------------------ equivalent loads in global coordinates */
if (ele->node[0]->gnode->ondesigntyp==ondnode && ele->node[0]->gnode->d.dnode->cos_type == 1)
{
  /* local cos */
  loadvec[0] = loadvec_local[0];
  loadvec[1] = loadvec_local[1];
  loadvec[2] = loadvec_local[2];
}
else
{
  /* global cos */
  loadvec[0] =  loadvec_local[0]*cosa + loadvec_local[1]*sina;/*load in horizontal direction*/
  loadvec[1] =  loadvec_local[0]*sina - loadvec_local[1]*cosa;/*load in vertical direction*/
  loadvec[2] = -loadvec_local[2];                             /*moment */     
}
if (ele->node[1]->gnode->ondesigntyp==ondnode && ele->node[1]->gnode->d.dnode->cos_type == 1)
{
  /* local cos */
  loadvec[3] = loadvec_local[3];
  loadvec[4] = loadvec_local[4];
  loadvec[5] = loadvec_local[5];
}
else
{
  /* global cos */
  loadvec[3] =  loadvec_local[3]*cosa + loadvec_local[4]*sina;/*load in horizontal direction*/
  loadvec[4] =  loadvec_local[3]*sina - loadvec_local[4]*cosa;/*load in vertical direction*/
  loadvec[5] = -loadvec_local[5];                             /*moment */     
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of saxi_eleload */

/*----------------------------------------------------------------------*/
#endif /*D_AXISHELL*/
/*! @} (documentation module close)*/
