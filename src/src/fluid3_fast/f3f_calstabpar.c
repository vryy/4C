/*-----------------------------------------------------------------------*/
/*!
\file
\brief Brief description.

  Very detailed description.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

#ifndef CCADISCRET
/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "f3f_prototypes.h"
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

static DOUBLE Q13  = ONE/THREE;
static DOUBLE Q112 = ONE/TWELVE;
static FLUID_DYNAMIC *fdyn;




/*-----------------------------------------------------------------------*/
/*!
  \brief routine to calculate stabilisation parameter


  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param velint         *DOUBLE   (i) vel at int point
  \param tau            *DOUBLE   (o) stabilisation parameter
  \param visc            DOUBLE   (i) viscosity
  \param inttyp          INT      (i) typ of integration
  \param iflag           INT      (i) flag for evaluation
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalstabpar(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *velint,
    DOUBLE          *tau,
    DOUBLE           visc,
    INT              inttyp,
    INT              iflag,
    INT              sizevec[6]
    )
{
  INT    j,isp;
  DOUBLE hdiv=ONE;
  DOUBLE velno;
  DOUBLE c_mk=0.0;
  DOUBLE dt;
  DOUBLE re;
  DOUBLE hk;
  DOUBLE aux1;
  STAB_PAR_GLS *gls;    /* pointer to GLS stabilisation parameters      */

#ifdef DEBUG
  dstrc_enter("f3_calstabpar");
#endif

  /* initialise */
  gls    = ele[0]->e.f3->stabi.gls;
  fdyn   = alldyn[genprob.numff].fdyn;

  if (ele[0]->e.f3->stab_type != stab_gls)
    dserror("routine with no or wrong stabilisation called");


  /* higher order element diameter modifications ? */
  for(j=0;j<sizevec[4];j++)
  {

    switch(gls->mk)
    {
      case -1:
        c_mk = Q13;
        if (inttyp==20 || inttyp==27)
        {
          if (sizevec[1]<32)
            hdiv = TWO;
          else
            hdiv = THREE;
        }
        else if (inttyp==10)
        {
          if (sizevec[1]==10)
            hdiv = TWO;
          else
            hdiv = THREE;
        }
        break;
      case 0:
        if (sizevec[1]>=8)
          c_mk=Q112;
        else
          c_mk=Q13;
        break;
      default:
        dserror("mk > 0 not implemented yet!");
    } /* end swtich (ele->e.f3->mk) */


    /* choose stability-parameter version */
    switch(gls->istapa)
    {
      case 35: /* version diss. Wall - instationary */
        velno = sqrt(velint[j]*velint[j] \
            +velint[LOOPL+j]*velint[LOOPL+j] \
            +velint[2*LOOPL+j]*velint[2*LOOPL+j]); /*norm of vel */
        dt = fdyn->dta;     /* check if dta or dt has to be chosen!!!!!!!! */
        for (isp=0;isp<3;isp++)
        {
          if (gls->itau[isp]!=iflag)
            continue;
          hk = ele[j]->e.f3->hk[isp]/hdiv;
          switch(isp)
          {
            case 2:/* continiuty stabilisation */
              re = c_mk*hk*velno/TWO/visc;  /* element reynolds number */
              /*         fdyn->tau[isp] = (ele[j]->e.f3->clamb)*velno*hk/TWO*DMIN(ONE,re);  */
              tau[isp*LOOPL+j] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
              break;
            default: /* velocity / pressure stabilisation */
              if (velno>EPS15)
              {
                aux1 = DMIN(hk/TWO/velno , c_mk*hk*hk/FOUR/visc);
                tau[isp*LOOPL+j] = DMIN(dt , aux1);
              }
              else
                tau[isp*LOOPL+j] = DMIN(dt , c_mk*hk*hk/FOUR/visc);
              break;
          } /* end switch (isp) */
        } /* end of loop over isp */
        break;

      case 36: /* version diss. Wall - stationary */
        dserror("stationary stabilisation not checked yet!!!");
        velno = sqrt(velint[0]*velint[0] \
            +velint[1]*velint[1] \
            +velint[2]*velint[2]); /*norm of vel */
        aux1= velno*c_mk/FOUR/visc;
        for (isp=0;isp<3;isp++)
        {
          if (gls->itau[isp]!=iflag)
            continue;
          hk = ele[j]->e.f3->hk[isp]/hdiv;
          re = aux1*hk;
          switch(isp)
          {
            case 2: /* continiuty stabilisation ### TWO VERSIONS ??? ###*/
              /*fdyn->tau[isp] = (ele->e.f3->clamb)*velno*hk/TWO*DMIN(ONE,re);*/
              tau[isp*LOOPL+j] = (gls->clamb)*velno*hk/TWO*DMIN(ONE,re);
              /*         fdyn->tau[isp] = velno*hk/TWO*DMIN(ONE,re);*/
              break;
            default: /* velocity / pressure stabilisation */
              if (re<ONE)
                /*fdyn->tau[isp] = c_mk*hk*hk/FOUR/visc;*/
                tau[isp*LOOPL+j] = c_mk*hk*hk/FOUR/visc;
              else
                /*fdyn->tau[isp] = hk/TWO/velno;*/
                tau[isp*LOOPL+j] = hk/TWO/velno;
              break;
          }  /* end switch (isp) */
        } /* end loop over isp */
        break;

      default:
        dserror("stability parameter version ISTAP unknown!");
    } /* end switch (ele->e.f3->istapa) */

  }/*loop*/

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3_calstabpar*/


#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/


#endif
