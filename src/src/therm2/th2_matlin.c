/*======================================================================*/
/*!
\file
\brief Linear heat conduction laws dubbed Fourier's law

\author bborn
\date 03/06
*/
#ifndef CCADISCRET
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"


/*======================================================================*/
/*!
\brief Construct constitutive matrix for isotropic linear Fourier law

\param con         DOUBLE           (i)   isotropic thermal conductivity
\param *ele        ELEMENT          (i)   pointer to current element
\param **bop       DOUBLE           (i)   B-operator at GP
\param *heatflux   DOUBLE           (o)   heat flux
\param **cmat      DOUBLE           (o)   constitutive matrix
                                          conductivity matrix
\return void

\author bborn
\date 03/06
*/
void th2_matlin_iso(DOUBLE con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat)
{
  TH2_PLANESTATES plst;  /* plane state */
  DOUBLE tmgr[NUMTMGR_THERM2];  /* temperature gradient vector */
  INT itmgr; /* temperature gradient counter */
  DOUBLE tmgrsum;  /* temp. grad. dummy sum */
  INT inode;  /* element node counter */
  INT ihflux;  /* heat flux component counter */
  DOUBLE hfluxsum;   /* dummy heat flux sum */


  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_matlin_iso");
#endif

  /*--------------------------------------------------------------------*/
  plst = ele->e.th2->planestate;
  switch (plst)
  {
    case th2_plane_tempgrad: case th2_plane_heatflux:
      /* temperature gradient at GP */
      for (itmgr=0; itmgr<NUMTMGR_THERM2; itmgr++)
      {
        tmgrsum = 0.0;
        for (inode=0; inode<ele->numnp; inode++)
        {
          tmgrsum += bop[itmgr][inode]*ele->node[inode]->sol.a.da[0][0];
        }
        tmgr[itmgr] = tmgrsum;
      }
      /* constitutive matrix */
      cmat[0][0] = -con;
      cmat[0][1] = 0.0;
      cmat[1][0] = 0.0;
      cmat[1][1] = -con;
      /* compute heat fluxes in xy-plane */
      for (ihflux=0; ihflux<NUMTMGR_THERM2; ihflux++)
      {
        hfluxsum = 0.0;
        for (itmgr=0; itmgr<NUMTMGR_THERM2; itmgr++)
        {
          hfluxsum += cmat[ihflux][itmgr] * tmgr[itmgr];
        }
        heatflux[ihflux] = hfluxsum;
      }
      /* compute heat flux perpendicular to xy-plane, which is zero for
       * both plane states
       * plane temper. grad. : zero due to diagonal constitutive matrix
       * planar heat flux : zero per definition
       */
      heatflux[NUMHFLX_THERM2-1] = 0.0;
      break;
  }  /* end of switch (plst) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief  Construct constitutive matrix for general linear Fourier law

\param **con       DOUBLE           (i)   thermal conductivity matrix
                                          general 3dim setup

\param *heatflux   DOUBLE           (o)   heatflux
\param **cmat      DOUBLE           (o)   constitutive matrix
                                          conductivity matrix

\author bborn
\date 03/06
*/
void th2_matlin_gen(DOUBLE **con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat)
{
  const INT heatminus = -1.0;  /* minus sign occuring in heat conduction */
  INT plst;  /* plane state */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_matlin_gen");
#endif


  /*--------------------------------------------------------------------*/
  /* distinguish plane state */
  plst = ele->e.th2->planestate;
  switch (plst)
  {
    case th2_plane_tempgrad:
      cmat[0][0] = heatminus * con[0][0];
      cmat[0][1] = heatminus * con[0][1];
      cmat[1][0] = heatminus * con[1][0];
      cmat[1][1] = heatminus * con[1][1];
      break;
    case th2_plane_heatflux:
      cmat[0][0] = heatminus * (con[0][0] + con[0][2]*con[2][0]/con[2][2]);
      cmat[0][1] = heatminus * (con[0][1] + con[0][2]*con[2][1]/con[2][2]);
      cmat[1][0] = heatminus * (con[1][0] + con[1][2]*con[2][0]/con[2][2]);
      cmat[1][1] = heatminus * (con[1][1] + con[1][2]*con[2][1]/con[2][2]);
      break;
  }  /* end of switch (plst) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
#endif  /* end of #ifdef D_THERM2 */

#endif
