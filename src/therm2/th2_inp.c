/*----------------------------------------------------------------------*/
/*!
\file
\brief Contains the routine 'th2_inp' which reads the planar thermal field
element (THERM2)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifdef D_THERM2


/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "therm2.h"


/*!
\addtogroup THERM2
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Input of therm2 element

This routines reads the therm2 elements. It is called by inpfield() 
(in input_mesh.c)

\param   ele     *ELEMENT   (o)   pointer to element

\author bborn
\date 03/06
*/
void th2_inp(ELEMENT *ele)
{
  INT i;
  INT ierr=0;
  char buffer[50];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_inp");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate the element */
  ele->e.th2 = (THERM2*) CCACALLOC(1, sizeof(THERM2));
  /* check allocation */
  if (ele->e.th2 == NULL)
  {
    dserror("Allocation of element failed");
  }

  /*--------------------------------------------------------------------*/
  /* read element topology */
  frchk("QUAD4", &ierr);
  if (ierr == 1)
  {
    ele->distyp = quad4;
    ele->numnp = 4;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("QUAD4", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of ELEMENT Topology failed");
    }
  }

  frchk("QUAD8", &ierr);
  if (ierr == 1)
  {
    ele->distyp = quad8;
    ele->numnp = 8;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("QUAD8", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of ELEMENT Topology failed");
    }
  }

  frchk("QUAD9", &ierr);
  if (ierr == 1)
  {
    ele->distyp = quad9;
    ele->numnp = 9;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("QUAD9", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1) 
    {
      dserror("Reading of ELEMENT Topology failed");
    }
  }

  frchk("TRI3", &ierr);
  if (ierr == 1)
  {
    ele->distyp = tri3;
    ele->numnp = 3;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed\n");
    }
    frint_n("TRI3", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of ELEMENT Topology failed\n");
    }
  }

  frchk("TRI6", &ierr);
  if (ierr == 1)
  {
    ele->distyp = tri6;
    ele->numnp = 6;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed\n");
    }
    frint_n("TRI6", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of ELEMENT Topology failed\n");
    }
  }

  /*--------------------------------------------------------------------*/
  /* reduce node numbers by one */
  for (i=0; i<ele->numnp; i++) 
  {
    (ele->lm[i])--;
  }

  /*--------------------------------------------------------------------*/
  /* read the material number */
  frint("MAT", &(ele->mat), &ierr);
  if (ierr != 1)
  {
    dserror("Reading of THERM2 element failed: MAT was not found!");
  }

  /*--------------------------------------------------------------------*/
  /* read element thickness */
  frdouble("THICK", &(ele->e.th2->thick), &ierr);
  if (ierr != 1)
  {
    dserror("Reading of THERM2 element failed: THICK was not found!");
  }

  /*--------------------------------------------------------------------*/
  /* read Gaussian points for THERM2 elements */
  frint_n("GP", &(ele->e.th2->nGP[0]), 2, &ierr);
  if (ierr != 1)
  {
    dserror("Reading of THERM2 element failed: GP was not found!");
  }
  th2_intg_eleinp(ele, &ierr);
  if (ierr != 1)
  {
    dserror("Gauss point combination is not possible!");
  }
  /*--------------------------------------------------------------------*/
  /* plane state */
  ele->e.th2->planestate = th2_plane_heatflux;
  frchk("PLANE_HEATFLUX", &ierr);
  if (ierr == 1)
  { 
    ele->e.th2->planestate = th2_plane_heatflux;
  }
  frchk("PLANE_TEMPGRAD", &ierr);
  if (ierr == 1)
  { 
    ele->e.th2->planestate = th2_plane_tempgrad;
  }
  /*--------------------------------------------------------------------*/
  /* read kinematic type */
  ele->e.th2->kintype = th2_geo_lin;   /* default */
  frchk("KIN_GEOLIN", &ierr);
  if (ierr == 1)
  {
    ele->e.th2->kintype = th2_geo_lin;
  }
  frchk("KIN_TOTLAG", &ierr);
  if (ierr == 1)
  {
    ele->e.th2->kintype = th2_total_lagr;
  }
  frchk("KIN_UPDLAG", &ierr);
  if (ierr == 1)
  {
    ele->e.th2->kintype = th2_updated_lagr;
    dserror("Updated Lagrange for THERM2 is not implemented!");
  }
  /*--------------------------------------------------------------------*/
  /* heat flux output type */
  /* set default */
  ele->e.th2->hfluxtype = th2_hflux_none;
  /* select output type */
  frchar("HFLUX", buffer, &ierr);
  if (ierr == 1)
  {
    /* at Gauss point */
    if (strncmp(buffer, "Gpxy", 4) == 0)
    {
      ele->e.th2->hfluxtype = th2_hflux_gpxy;
    }
    if (strncmp(buffer, "Gp12", 4) == 0)
    {
      ele->e.th2->hfluxtype = th2_hflux_gp12;
    }
    if (strncmp(buffer, "Gpxy12", 6) == 0)
    {
      ele->e.th2->hfluxtype = th2_hflux_gpxy12;
    }
    /* at element nodes */
    if (strncmp(buffer, "Ndxy", 4) == 0)
    {
      ele->e.th2->hfluxtype = th2_hflux_ndxy;
    }
    if (strncmp(buffer, "Nd12", 4) == 0)
    {
      ele->e.th2->hfluxtype = th2_hflux_nd12;
    }
    if (strncmp(buffer, "Ndxy12", 6) == 0)
    {
      ele->e.th2->hfluxtype = th2_hflux_ndxy12;
    }
  }

  /*--------------------------------------------------------------------*/
  /* finalise */

#ifdef DEBUG
  dstrc_exit();
#endif

return;}
/* end of th2_inp() */

/*======================================================================*/
#endif /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
