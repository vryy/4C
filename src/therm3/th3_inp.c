/*----------------------------------------------------------------------*/
/*!
\file
\brief Input of 3dim thermal element (THERM3)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Input of therm3 element

This routines reads the therm3 element description in input file. 
It is called by inpfield() (in input_mesh.c)

\param   ele     *ELEMENT   (o)   pointer to element

\author bborn
\date 09/06
*/
void th3_inp(ELEMENT *ele)
{
  INT i;
  INT ierr = 0;
  char buffer[50];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_inp");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate the element */
  ele->e.th3 = (THERM3*) CCACALLOC(1, sizeof(THERM3));
  /* check allocation */
  if (ele->e.th3 == NULL)
  {
    dserror("Allocation of element failed");
  }

  /*--------------------------------------------------------------------*/
  /* read element topology */
  frchk("HEX8", &ierr);  /* 8-node hexahedron */
  if (ierr == 1)
  {
    ele->distyp = hex8;
    ele->numnp = 8;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("HEX8", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of THERM3 HEX8 element topology failed");
    }
  }

  frchk("HEX20", &ierr);  /* 20-node hexahedron */
  if (ierr == 1)
  {
    ele->distyp = hex20;
    ele->numnp = 20;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("HEX20", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of THERM3 HEX20 element topology failed");
    }
  }
  
  frchk("HEX27", &ierr);  /* 27-node hexahedron */
  if (ierr==1)
  {
    ele->distyp = hex27;
    ele->numnp = 27;
    ele->lm = (INT*) CCACALLOC(ele->numnp,sizeof(INT));
    if (ele->lm==NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("HEX27", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1) 
    {
      dserror("Reading of THERM3 HEX27 element topology failed");
    }
  }

  frchk("TET4", &ierr);  /* 4-node tetrahedron */
  if (ierr == 1)
  {
    ele->distyp = tet4;
    ele->numnp = 4;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("TET4", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of THERM3 TET10 element topology failed");
    }
  }

  frchk("TET10", &ierr);  /* 10-node tetrahedron */
  if (ierr == 1)
  {
    ele->distyp = tet10;
    ele->numnp = 10;
    ele->lm = (INT*) CCACALLOC(ele->numnp, sizeof(INT));
    if (ele->lm == NULL)
    {
      dserror("Allocation of lm in ELEMENT failed");
    }
    frint_n("TET10", &(ele->lm[0]), ele->numnp, &ierr);
    if (ierr != 1)
    {
      dserror("Reading of THERM3 TET10 element topology failed");
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
    dserror("Reading of THERM3 element failed");
  }
  
  /*--------------------------------------------------------------------*/
  /* read the gaussian points */
  frint_n("GP", &(ele->e.th3->gpnum[0]), NDIM_THERM3, &ierr);
  if (ierr != 1)
  {
    dserror("Reading of THERM3 element failed: GP not found!");
  }
  th3_intg_eleinp(ele, &ierr);
  if (ierr != 1)
  {
    dserror("Gauss point combination is not possible!");
  }

  /*--------------------------------------------------------------------*/
  /* read kinematic type */
  ele->e.th3->kintype = th3_geo_lin;   /* default */
  frchar("KINEM", buffer, &ierr);
  if (ierr == 1)
  {
    /* geometrically linear */
    if (strncmp(buffer,"Geolin",6) == 0)
    {
      ele->e.th3->kintype = th3_geo_lin;
    }
    /* geometrically non-linear with Total Lagrangean approach */
    if (strncmp(buffer,"Totlag",6) == 0)
    {
      ele->e.th3->kintype = th3_total_lagr;
    }
    /* geometrically non-linear with Updated Lagrangean approach */
    if (strncmp(buffer,"Updlag",6) == 0)
    {
      ele->e.th3->kintype = th3_updated_lagr;
      dserror("Updated Lagrange for THERM3 is not implemented!");
    }
  }

  /*--------------------------------------------------------------------*/
  /* read local or global stresses */
  ele->e.th3->hfluxtype = th3_hflux_none; /* set default */
  frchar("HFLUX", buffer, &ierr);
  if (ierr == 1)
  {
    if (strncmp(buffer,"None",4) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_none;
    }
    if (strncmp(buffer,"Gpxyz",5) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_gpxyz;
    }
    if (strncmp(buffer,"Gprst",5) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_gprst;
    }
    if (strncmp(buffer,"Gp123",5) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_gp123;
    }
    if (strncmp(buffer,"Ndxyz",5) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_ndxyz;
    }
    if (strncmp(buffer,"Ndrst",5) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_ndrst;
    }
    if (strncmp(buffer,"Nd123",5) == 0)
    {
      ele->e.th3->hfluxtype = th3_hflux_nd123;
    }
  }


  /*--------------------------------------------------------------------*/
  /* finalise */
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of th3_inp */

#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/

