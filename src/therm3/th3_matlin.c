/*======================================================================*/
/*!
\file
\brief Linear heat conduction laws dubbed Fourier's law

\author bborn
\date 03/06
*/
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"


/*======================================================================*/
/*!
\brief Construct constitutive matrix for isotropic linear Fourier law

\param con         DOUBLE           (i)   isotropic thermal conductivity
\param *ele        ELEMENT          (i)   pointer to current element
\param **bop       DOUBLE           (i)   B-operator at GP
\param *heatflux   DOUBLE           (io)  heat flux
\param **cmat      DOUBLE           (o)   constitutive matrix
                                          conductivity matrix
\return void

\author bborn
\date 09/06
*/
void th2_matlin_iso(DOUBLE con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat)
{
  DOUBLE tmgr[NUMTMGR_THERM3];  /* temperature gradient vector */
  INT itmgr; /* temperature gradient counter */
  DOUBLE tmgrsum;  /* temp. grad. dummy sum */
  INT inode;  /* element node counter */
  INT ihflux;  /* heat flux component counter */
  DOUBLE hfluxsum;   /* dummy heat flux sum */
  

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_matlin_iso");
#endif

  /*--------------------------------------------------------------------*/
  /* temperature gradient at GP */
  for (itmgr=0; itmgr<NUMTMGR_THERM3; itmgr++)
  {
    tmgrsum = 0.0;
    for (inode=0; inode<ele->numnp; inode++)
    {
      tmgrsum += bop[itmgr][inode] * ele->node[inode]->sol.a.da[0][0];
    }
    tmgr[itmgr] = tmgrsum;
  }
  /* constitutive matrix */
  cmat[0][0] = -con;
  cmat[0][1] = 0.0;
  cmat[0][2] = 0.0;
  cmat[1][0] = 0.0;
  cmat[1][1] = -con;
  cmat[1][2] = 0.0;
  cmat[2][0] = 0.0;
  cmat[2][1] = 0.0;
  cmat[2][2] = -con;
  /* compute heat fluxes in xy-plane */
  for (ihflux=0; ihflux<NUMHFLX_THERM3; ihflux++)
  {
    hfluxsum = 0.0;
    for (itmgr=0; itmgr<NUMTMGR_THERM3; itmgr++)
    {
      hfluxsum += cmat[ihflux][itmgr] * tmgr[itmgr];
    }
    heatflux[ihflux] = hfluxsum;
  }

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

\param *heatflux   DOUBLE           (io)  heatflux
\param **cmat      DOUBLE           (o)   constitutive matrix
                                          conductivity matrix

\author bborn
\date 03/06
*/
void th3_matlin_gen(DOUBLE **con,
                    ELEMENT *ele,
                    DOUBLE **bop,
                    DOUBLE *heatflux,
                    DOUBLE **cmat)
{
  const INT heatminus = -1.0;  /* minus sign occuring in heat conduction */
  INT itmgr; /* temperature gradient counter */
  DOUBLE tmgrsum;  /* temp. grad. dummy sum */
  INT ihflx;  /* heat flux component counter */
  DOUBLE hfluxsum;   /* dummy heat flux sum */
  INT inode;  /* element node counter */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_matlin_gen");
#endif

  /*--------------------------------------------------------------------*/
  /* temperature gradient at GP */
  for (itmgr=0; itmgr<NUMTMGR_THERM3; itmgr++)
  {
    tmgrsum = 0.0;
    for (inode=0; inode<ele->numnp; inode++)
    {
      tmgrsum += bop[itmgr][inode] * ele->node[inode]->sol.a.da[0][0];
    }
    tmgr[itmgr] = tmgrsum;
  }
  /*--------------------------------------------------------------------*/
  /* conductivity matrix */
  for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
  {
    for (itmgr=0; itmgr<NUMTMGR_THERM3; itmgr++)
    {
      cmat[ihflx][itmgr] = heatminus * con[ihflx][itmgr];
    }
  }
  /*--------------------------------------------------------------------*/
  /* compute heat fluxes */
  for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
  {
    hfluxsum = 0.0;
    for (itmgr=0; itmgr<NUMTMGR_THERM3; itmgr++)
    {
      hfluxsum += cmat[ihflx][itmgr] * tmgr[itmgr];
    }
    heatflux[ihflx] = hfluxsum;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close) */
