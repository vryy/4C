/*======================================================================*/
/*!
\file
\brief Integration parameters for THERM3 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 09/06
*/
/*======================================================================*/
/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/* limit to D_THERM3 */
#ifdef D_THERM3


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "therm3.h"


/*----------------------------------------------------------------------*/
/*!
\brief general problem data

global variable GENPROB genprob is defined in global_control.c

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
extern GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief Pointer to allocate dynamic variables if needed

defined in global_control.c

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
extern ALLDYNA *alldyn;


/*----------------------------------------------------------------------*/
/*!
\brief Set type, total amount of Gauss points, etc of element

\param ELEMENT    *ele      (i)   pointer to element
\parem INT         ierr     (o)   error flag: 1 ==> success 

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
void th3_intg_eleinp(ELEMENT *actele, 
                     INT *ierr)
{

  INT i, j;
  INT *actGP;
  INT *intc;
  THERM3 *actth3;

#ifdef DEBUG
  dstrc_enter("th3_intgeleinp");
#endif

  /*--------------------------------------------------------------------*/
  /* init error flag */
  *ierr = 1;  /* start without an error */

  /*--------------------------------------------------------------------*/
  /* get local copies / pointers */
  actth3 = actele->e.th3;
  actGP = &(actth3->gpnum[0]);
  intc = &(actth3->gpintc[0]);

  /*--------------------------------------------------------------------*/
  /* set Gauss number type and total Gauss points */
  switch (actele->distyp)
  {
    /*------------------------------------------------------------------*/
    /* hexahedra */
    case hex8: case hex20: case hex27:
      /* check if first Gauss point number is set larger than 0 */
      if ( (actGP[0] < 1) || (actGP[0] > GLINTC_THERM3) )
      {
        dserror("First GP must be 1, 2, 3, 4, 5 or 6 (hexahedra)"); 
      }
      /* set Gauss point numbers in s- and t-direction if not set */
      for (i=1; i<NDIM_THERM3; i++)
      {
        if ( (actGP[i] < 1) || (actGP[i] > GLINTC_THERM3) )
        {
          actGP[i] = actGP[0];
        }
      }
      /* reduce numbers by 1 to achieve 0-based indices */
      for (i=0; i<NDIM_THERM3; i++)
      {
        intc[i] = actGP[i] - 1;
      }
      break;
    /*------------------------------------------------------------------*/
    /* tetrahedra */
    case tet4: case tet10:
      /*----------------------------------------------------------------*/
      /* 1st stage : verify/set Gauss point numbers */
      /* check if first Gauss point number (for domain) is permissible */
      if ( !( (actGP[0] == 1) || (actGP[0] == 4) || (actGP[0] == 5) ) )
      {
        dserror("First GP must be 1, 4, or 5 (total GP of tetrahedra domain)"); 
      }
      /* set 2nd Gauss point number (for sides) */
      if ( !( (actGP[1] == 1) || (actGP[1] == 3) 
              || (actGP[1] == 4) ||  (actGP[1] == 6) ) )
      {
        if (actGP[0] == 1)
        {
          actGP[1] = 1;
        }
        else if (actGP[0] == 4)
        {
          actGP[1] = 3;
        }
        else if (actGP[0] == 5)
        {
          actGP[1] = 4;
        }
      }
      /* set 2nd Gauss point number (for edges) */
      if ( (actGP[2] < 1) || (actGP[2] > GLINTC_THERM3) )
      {
        if (actGP[0] == 1)
        {
          actGP[2] = 1;
        }
        else if (actGP[0] == 4)
        {
          actGP[2] = 2;
        }
        else if (actGP[0] == 5)
        {
          actGP[2] = 2;
        }
      }
      /*----------------------------------------------------------------*/
      /* 2nd stage :  associate indices to Gauss point numbers */
      /* Domain Gauss point numbers to index of corresponding GP set */
      switch (actGP[0])
      {
          case 1: 
            gintc[0] = 0;
            break;
          case 4:
            intc[0] = 1;
            break;
          case 5:
            intc[0] = 2;
            break;
          default:
            dserror("Impermissible Gauss point number");
            break;
      }
      /* Side Gauss point numbers to index of corresponding GP set */
      switch (actGP[1])
      {
          case 1:
            intc[1] = 0;
            break;
          case 3:
            intc[1] = 1;
            break;
          case 4:
            intc[1] = 3;
            break;
          case 6:
            intc[1] = 4;
            break;
          default:
            dserror("Impermissible Gauss point number");
            break;
      }
      /* Edge Gauss point numbers to index of corresponding GP set */
      intc[2] = actGP[2] - 1;
      break;
    /*------------------------------------------------------------------*/
    /* catch erroneous disc. types */
    default:
      dserror("Discretisation type is impossible!");
      break;
  }  /* end of switch (actele->distyp) */

  /*--------------------------------------------------------------------*/
  /* finish off */
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}  /* end of th3_intg_eleinp(...)  */


/*----------------------------------------------------------------------*/
/*!
\brief integration parameters of THERM3 element

this routine is a try to organise the integration parameters
different. ALL paramters are stored in FLUID_DATA, so that this
routine has to be (hopefully) called only once!!!

\param  *data   TH3_DATA  (o)  prepare quadrature data 
\return void

\author bborn
\date 09/06
*/
/*----------------------------------------------------------------------*/
void th3_intg_init(TH3_DATA *data)

{
  INT i, j, k;  /* counters */

#ifdef DEBUG
  dstrc_enter("th3_intg");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize arrays */
  /* hexahedra and tetrahedra lines */
  for (k=0; k<GLINTC_THERM3; k++)
  {
    for (i=0; i<GLMAXP_THERM3; i++)
    {
      data.ghlc[k][i] = 0.0;
      data.ghlw[k][i] = 0.0;
      data.gtlc[k][i] = 0.0;
      data.gtlw[k][i] = 0.0;
    }
  }
  /* tetrahedra domains */
  for (k=0; k<GTINTC_THERM3; k++)
  {
    for (i=0; i<GTMAXP_THERM3; i++)
    {
      data.gtdcr[k][i] = 0.0;
      data.gtdcs[k][i] = 0.0;
      data.gtdcr[k][i] = 0.0;
      data.gtdw[k][i] = 0.0;
    }
  }
  /* tetrahedra sides */
  for (k=0; k<GSINTC_THERM3; k++)
  {
    for (i=0; i<GSMAXP_THERM3; i++)
    {
      data.gsdcr[k][i] = 0.0;
      gdata.sdcs[k][i] = 0.0;
      data.gsdw[k][i] = 0.0;
    }
  }

  /*--------------------------------------------------------------------*/
  /* hexahedra lines [-1,+1] Gauss points and weights */
  /* 1 Gauss point (1st order) */
  if ( (GLINTC_THERM3 >= 1) && (GLMAXP_THERM3 >= 1) )
  {
    data.ghlc[0][0] = 0.0;  /* r-, s- or t-coordinate respectively */
    data.ghlw[0][0] = 2.0;  /* weight */
  }
  /* 2 Gauss points (3rd order)  */
  if ( (GLINTC_THERM3 >= 2) && (GLMAXP_THERM3 >= 2) )
  {
    data.ghlc[1][0] = -1.0/sqrt(3.0);
    data.ghlc[1][1] = 1.0/sqrt(3.0);
    data.ghlw[1][0] = 1.0;
    data.ghlw[1][1] = 1.0;
  }
  /* 3 Gauss points (5th order) */
  if ( (GLINTC_THERM3 >= 3) && (GLMAXP_THERM3 >= 3) )
  {
    data.ghlc[2][0] = -sqrt(3.0/5.0);
    data.ghlc[2][1] = 0.0;
    data.ghlc[2][2] = sqrt(3.0/5.0);
    data.ghlw[2][0] = 5.0/9.0;
    data.ghlw[2][1] = 8.0/9.0;
    data.ghlw[2][2] = 5.0/9.0;
  }
  /* 4 Gauss points (7th order) */
  if ( (GLINTC_THERM3 >= 4) && (GLMAXP_THERM3 >= 4) )
  {
    data.ghlc[3][0] = -sqrt((15.0+sqrt(120.0))/35.0);
    data.ghlc[3][1] = -sqrt((15.0-sqrt(120.0))/35.0);
    data.ghlc[3][2] = sqrt((15.0-sqrt(120.0))/35.0);
    data.ghlc[3][3] = sqrt((15.0+sqrt(120.0))/35.0);
    data.ghlw[3][0] = (18.0-sqrt(30.0))/36.0;
    data.ghlw[3][1] = (18.0+sqrt(30.0))/36.0;
    data.ghlw[3][2] = (18.0+sqrt(30.0))/36.0;
    data.ghlw[3][3] = (18.0-sqrt(30.0))/36.0;
  }
  /* 5 Gauss points (9th order) */
  if ( (GLINTC_THERM3 >= 5) && (GLMAXP_THERM3 >= 5) )
  {
    data.ghlc[4][0] = -0.9061798459387;
    data.ghlc[4][1] = -0.5384693101057;
    data.ghlc[4][2] = 0.0;
    data.ghlc[4][3] = 0.5384693101057;
    data.ghlc[4][4] = 0.9061798459387;
    data.ghlw[4][0] = 0.2369268850562;
    data.ghlw[4][1] = 0.4786286704994;
    data.ghlw[4][2] = 0.5688888888889;
    data.ghlw[4][3] = 0.4786286704994;
    data.ghlw[4][4] = 0.2369268850562;
  }
  /* 6 Gauss points (11th order) */
  if ( (GLINTC_THERM3 >= 6) && (GLMAXP_THERM3 >= 6) )
  {
    data.ghlc[5][0] = -0.9324695142032;
    data.ghlc[5][1] = -0.6612093864663;
    data.ghlc[5][2] = -0.2386191860832;
    data.ghlc[5][3] = 0.2386191860832;
    data.ghlc[5][4] = 0.6612093864663;
    data.ghlc[5][5] = 0.9324695142032;
    data.ghlw[5][0] = 0.1713244923792;
    data.ghlw[5][1] = 0.3607615730481;
    data.ghlw[5][2] = 0.4679139345727;
    data.ghlw[5][3] = 0.4679139345727;
    data.ghlw[5][3] = 0.3607615730481;
    data.ghlw[5][5] = 0.1713244923792;
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron domains */
  /* 1 Gauss point (1st order) */
  if ( (GTINTC_THERM3 >= 1) && (GTMAXP_THERM3 >= 1) )
  {
    data.gtdcr[0][0] = 0.25;  /* r-coordinate */
    data.gtdcs[0][0] = 0.25;  /* s-coordinate */
    data.gtdct[0][0] = 0.25;  /* t-coordinate */
    data.gtdw[0][0]  = 1.0/6.0;  /* weight */
  }
  /* 4 Gauss points (2nd order) */
  if ( (GTINTC_THERM3 >= 2) && (GTMAXP_THERM3 >= 4) )
  {
    data.gtdcr[1][0] = 0.13819660;
    data.gtdcs[1][0] = 0.13819660;
    data.gtdct[1][0] = 0.13819660;
    data.gtdw[1][0] = 1.0/24.0;
    data.gtdcr[1][1] = 0.58541020;
    data.gtdcs[1][1] = 0.13819660;
    data.gtdct[1][1] = 0.13819660;
    data.gtdw[1][1] = 1.0/24.0;
    data.gtdcr[1][2] = 0.13819660;
    data.gtdcs[1][2] = 0.58541020;
    data.gtdct[1][2] = 0.13819660;
    data.gtdw[1][2] = 1.0/24.0;
    data.gtdcr[1][3] = 0.13819660;
    data.gtdcs[1][3] = 0.13819660;
    data.gtdct[1][3] = 0.58541020;
    data.gtdw[1][3] = 1.0/24.0;
  }
  /* 5 Gauss points (3rd order) */
  if ( (GTINTC_THERM3 >= 3) && (GTMAXP_THERM3 >= 5) )
  {
    data.gtdcr[2][0] = 1.0/6.0;
    data.gtdcs[2][0] = 1.0/6.0;
    data.gtdct[2][0] = 1.0/6.0;
    data.gtdw[2][0] = 1.0/36.0;
    data.gtdcr[2][1] = 1.0/2.0;
    data.gtdcs[2][1] = 1.0/6.0;
    data.gtdct[2][1] = 1.0/6.0;
    data.gtdw[2][1] = 1.0/36.0;
    data.gtdcr[2][2] = 1.0/6.0;
    data.gtdcs[2][2] = 1.0/2.0;
    data.gtdct[2][2] = 1.0/6.0;
    data.gtdw[2][2] = 1.0/36.0;
    data.gtdcr[2][3] = 1.0/6.0;
    data.gtdcs[2][3] = 1.0/6.0;
    data.gtdct[2][3] = 1.0/2.0;
    data.gtdw[2][3] = 1.0/36.0;
    data.gtdcr[2][4] = 1.0/4.0;
    data.gtdcs[2][4] = 1.0/4.0;
    data.gtdct[2][4] = 1.0/4.0;
    data.gtdw[2][4] = -2.0/15.0;
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron sides */
  /* 1 Gauss point (1st order) */
  if ( (GSINTC_THERM3 >= 1) && (GSMAXP_THERM3 >= 1) )
  {
    data.gtscr[0][0] = 1.0/3.0;  /* r-coordinate (or other) */
    data.gtscs[0][0] = 1.0/3.0;  /* s-coordinate (or other) */
    data.gtsw[0][0] = 1.0/2.0;  /* weight */
  }
  /* 3 Gauss points (2nd order) -- kind 1 */
  if ( (GTINTC_THERM3 >= 2) && (GTMAXP_THERM3 >= 3) )
  {
    data.gtscr[1][0] = 1.0/6.0;
    data.gtscs[1][0] = 1.0/6.0;
    data.gtsw[1][0] = 1.0/6.0;
    data.gtscr[1][1] = 2.0/3.0;
    data.gtscs[1][1] = 1.0/6.0; 
    data.gtsw[1][1] = 1.0/6.0;
    data.gtscr[1][2] = 1.0/6.0;
    data.gtscs[1][2] = 2.0/3.0;
    data.gtsw[1][2] = 1.0/6.0;
  }
  /* 3 Gauss points (2nd order) -- kind 2 -- unused */
  if ( (GTINTC_THERM3 >= 3) && (GTMAXP_THERM3 >= 3) )
  {
    data.gtscr[2][0] = 0.5;
    data.gtscs[2][0] = 0.0;
    data.gtsw[2][0] = 1.0/6.0;
    data.gtscr[2][1] = 0.5;
    data.gtscs[2][1] = 0.5;
    data.gtsw[2][1] = 1.0/6.0;
    data.gtscr[2][2] = 0.0;
    data.gtscs[2][2] = 0.5;
    data.gtsw[2][2] = 1.0/6.0;
  }
  /* 4 Gauss points (3rd order) */
  if ( (GTINTC_THERM3 >= 4) && (GTMAXP_THERM3 >= 4) )
  {
    data.gtscr[3][0] = 0.2;
    data.gtscs[3][0] = 0.2;
    data.gtsw[3][0] = 25.0/96.0;
    data.gtscr[3][1] = 0.6;
    data.gtscs[3][1] = 0.2;
    data.gtsw[3][1] = 25.0/96.0;
    data.gtscr[3][2] = 0.2;
    data.gtscs[3][2] = 0.6;
    data.gtsw[3][2] = 25.0/96.0;
    data.gtscr[3][3] = 1.0/3.0;
    data.gtscs[3][3] = 1.0/3.0;
    data.gtsw[3][3] = 9.0/32.0;
  }
  /* 6 Gauss points (4th order) */
  if ( (GTINTC_THERM3 >= 5) && (GTMAXP_THERM3 >= 6) )
  {
    data.gtscr[4][0] = 0.091576213509771;
    data.gtscs[4][0] = 0.091576213509771;
    data.gtsw[4][0] = 0.05497587182766;
    data.gtscr[4][1] = 0.81684757298045;
    data.gtscs[4][1] = 0.091576213509771;
    data.gtsw[4][1] = 0.054975871827661;
    data.gtscr[4][2] = 0.09157621350977;
    data.gtscs[4][2] = 0.816847572980459;
    data.gtsw[4][2] = 0.054975871827661;
    data.gtscr[4][3] = 0.445948490915965;
    data.gtscs[4][3] = 0.108103018168070;
    data.gtsw[4][3] = 0.111690794839006;
    data.gtscr[4][4] = 0.445948490915965;
    data.gtscs[4][4] = 0.445948490915965;
    data.gtsw[4][4] = 0.111690794839006;
    data.gtscr[4][5] = 0.108103018168070;
    data.gtscs[4][5] = 0.445948490915965;
    data.gtsw[4][5] = 0.111690794839006;
  }
  /*--------------------------------------------------------------------*/
  /* tetrahedron edges [0,1] */
  /* similar to hexahedra but on different parameter space */
  for (j=0; j<GLMAXP_THERM3; j++)
  {
    for (i=0; i<=j; i++)
    {
      /* r-coordinate (or other) */
      data.gtlc[j][i] = (1.0 + data.ghlc[j][i])/2.0;
      /* weight */
      data.gtlw[j][i] = data.ghlw[j][i]/2.0;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of th3_intg_init() */




/*----------------------------------------------------------------------*/
#endif  /* end of #ifdef D_THERM3 */

/*! @} (documentation module close)*/
