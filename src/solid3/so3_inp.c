/*----------------------------------------------------------------------*/
/*!
\file
\brief Input of 3dim thermal element (SOLID3)

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>
*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Input of solid3 element

This routines reads the solid3 element description in input file. 
It is called by inpfield() (in input_mesh.c)

\param   ele     *ELEMENT   (o)   pointer to element

\author mf
\date 10/06
*/
void so3_inp(ELEMENT *ele)
{
  INT i;
  INT ierr = 0;
  char buffer[50];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_inp");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate the element */
  ele->e.so3 = (SOLID3*) CCACALLOC(1, sizeof(SOLID3));
  /* check allocation */
  if (ele->e.so3 == NULL)
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
      dserror("Reading of ELEMENT Topology failed");
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
      dserror("Reading of ELEMENT Topology failed");
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
      dserror("Reading of ELEMENT Topology failed");
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
      dserror("Reading of ELEMENT Topology failed");
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
      dserror("Reading of ELEMENT Topology failed");
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
    dserror("Reading of MAT for SOLID3 element failed");
  }
  
  /*--------------------------------------------------------------------*/
  /* read the gaussian points */
  frint_n("GP", &(ele->e.so3->gpnum[0]), NDIM_SOLID3, &ierr);
  if (ierr != 1)
  {
    dserror("Reading of SOLID3 element failed: GP not found!");
  }
  so3_intg_eleinp(ele, &ierr);
  if (ierr != 1)
  {
    dserror("Gauss point combination is not possible!");
  }

  /*--------------------------------------------------------------------*/
  /* read kinematic type */
  frchar("KINEM", buffer, &ierr);
  if (ierr == 1)
  {
    /* geometrically linear */
    if (strncmp(buffer,"Geolin",6) == 0)
    {
      ele->e.so3->kintype = so3_geo_lin;
    }
    /* geometrically non-linear with Total Lagrangean approach */
    if (strncmp(buffer,"Totlag",6) == 0)
    {
      ele->e.so3->kintype = so3_total_lagr;
    }
    /* geometrically non-linear with Updated Lagrangean approach */
    if (strncmp(buffer,"Updlag",6) == 0)
    {
      ele->e.so3->kintype = so3_updated_lagr;
      dserror("Updated Lagrange for SOLID3 is not implemented!");
    }
  }
  else
  {
    /* default */
    ele->e.so3->kintype = so3_total_lagr; 
  }

  /*--------------------------------------------------------------------*/
  /* read local or global stresses */
  frchar("STRESS", buffer, &ierr);
  if (ierr == 1)
  {
    if (strncmp(buffer,"None",4) == 0)
    {
      ele->e.so3->stresstype = so3_stress_none;
    }
    if (strncmp(buffer,"Gpxyz",5) == 0)
    {
      ele->e.so3->stresstype = so3_stress_gpxyz;
    }
    if (strncmp(buffer,"Gprst",5) == 0)
    {
      ele->e.so3->stresstype = so3_stress_gprst;
    }
    if (strncmp(buffer,"Gp123",5) == 0)
    {
      ele->e.so3->stresstype = so3_stress_gp123;
    }
    if (strncmp(buffer,"Ndxyz",5) == 0)
    {
      ele->e.so3->stresstype = so3_stress_ndxyz;
    }
    if (strncmp(buffer,"Ndrst",5) == 0)
    {
      ele->e.so3->stresstype = so3_stress_ndrst;
    }
    if (strncmp(buffer,"Nd123",5) == 0)
    {
      ele->e.so3->stresstype = so3_stress_nd123;
    }
  }
  else
  {
    /* set default */
    ele->e.so3->stresstype = so3_stress_none;
  }

  /*--------------------------------------------------------------------*/
  /* read TSI coupling */
#ifdef D_TSI
  ele->e.so3->tsi_couptyp = tsi_coup_none;  /* default */
  frchar("TSI_COUPTYP",buffer,&ierr);
  if (ierr)
  {
    if (strncmp(buffer,"None",4)==0)
    {
      ele->e.so3->tsi_couptyp = tsi_coup_none;
    }
    if (strncmp(buffer,"Thermconf",9)==0)
    {
      ele->e.so3->tsi_couptyp = tsi_coup_thermconf;
    }
    if (strncmp(buffer,"Thermcreate",11)==0)
    {
      ele->e.so3->tsi_couptyp = tsi_coup_thermcreate;
    }
  }
#endif

  /*--------------------------------------------------------------------*/
  /* finalise */
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of so3_inp */

#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/

