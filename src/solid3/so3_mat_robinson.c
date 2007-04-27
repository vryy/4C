/*======================================================================*/
/*!
\file
\brief Robinson's visco-plastic material

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de/Members/bornemann
	    089-289-15237
</pre>

\author bborn
\date 03/07
*/
#ifdef D_SOLID3
#ifdef D_TSI

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/07
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Alldyn dynamic control

pointer to allocate dynamic variables if needed
dedfined in global_control.c

\auther bborn
\date 03/06
*/
extern ALLDYNA* alldyn;

/*======================================================================*/
/*!
\brief Allocate material internal variables (MIVs)
\author bborn
\date 03/07
*/
void so3_mat_robinson_init(ELEMENT* ele)  /*!< current element */
{
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_init");
#endif

  /*--------------------------------------------------------------------*/
  /* check if we are in a sane context */
  if (genprob.probtyp != prb_tsi)
  {
    dserror("Trying to use Robinson's material in non-TSI calculation!");
  }

  /*--------------------------------------------------------------------*/
  /* allocate MIV */
  if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
  {
    so3_mat_robinson_fb4_init(actso3);
  }
  else if (tsidyn->kind == tsi_therm_stat_struct_genalp)
  {
    so3_mat_robinson_be_init(actso3);
  }
  else
  {
    dserror("Trying to use Robinson's material with unimplemented time integration");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end so3_mat_robinson_init */


/*======================================================================*/
/*!
\brief Deallocate material internal variables (MIVs)
\author bborn
\date 03/07
*/
void so3_mat_robinson_final(ELEMENT* ele)  /*!< current element */
{
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */ 
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocate MIV */
  if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
  {
    so3_mat_robinson_fb4_final(actso3);
  }
  else if (tsidyn->kind == tsi_therm_stat_struct_genalp)
  {
    so3_mat_robinson_be_final(actso3);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end so3_mat_robinson_init */

/*======================================================================*/
/*!
\brief Get temperature-dependent material parameter at current temperature
\param   ipl    MAT_PARAM_INTPOL  (i)  interpolation type
\param   prm_n      INT           (i)  magnitude param data
\param   prm        DOUBLE*       (i)  parameter data
\param   tmpr       DOUBLE        (i)  curr. temperature
\param   prmbytempr DOUBLE*       (o)   param. at current temperature
\return  void
\author bborn
\date 04/07
*/
void so3_mat_robinson_prmbytmpr(const MAT_PARAM_INTPOL ipl,
                                const INT prm_n,
                                const DOUBLE* prm,
                                const DOUBLE tmpr,
                                DOUBLE* prmbytempr)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_prmbytmpr");
#endif

  /*--------------------------------------------------------------------*/
  /* constant */
  if (ipl == mat_param_ipl_const)
  {
    *prmbytempr = prm[0];
  }
  /* polynomial */
  else if (ipl == mat_param_ipl_poly)
  {
    *prmbytempr = 0.0;  /* initialise */
    INT tmpr_pow = 1.0;
    INT i;
    for (i=0; i<prm_n; ++i)
    {
      *prmbytempr += prm[i] * tmpr_pow;
      tmpr_pow *= tmpr;
    }
  }
  /* piece-wise linear */
  else if (ipl == mat_param_ipl_pcwslnr)
  {
    /* constant if lower than smallest provided temperature */
    if (tmpr <= prm[0])
    {
      *prmbytempr = prm[1];
    }
    /* constant if greater than largest provided temperature */
    else if (tmpr > prm[prm_n-2])
    {
      *prmbytempr = prm[prm_n-1];
    }
    /* linear interpolation */
    else
    {
      INT i;
      for (i=0; i<prm_n-2; i+=2)
      {
        if (tmpr <= prm[i+2])
        {
          /* we got the correct interval */
          DOUBLE x1 = prm[i];
          DOUBLE y1 = prm[i+1];
          DOUBLE x2 = prm[i+2];
          DOUBLE y2 = prm[i+3];
          DOUBLE x = tmpr;
          *prmbytempr = (y2*(x2-x) - y1*(x1-x))/(x2-x1);
          break;
        }
      }
    }
  }
  else if (ipl == mat_param_ipl_none)
  {
    dserror("Interpolation type is unknown!");
  }
  else
  {
    dserror("Interpolation type is unknown!");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Select Robinson's material and integrate internal variables
\param   container   CONTAINER*     (i)   container.h data
\param   ele         ELEMENT*       (i)   curr. element
\param   mat_robin   VP_ROBINSON*   (i)   curr. Robinson material data
\param   ip          INT            (i)   total GP index
\param   gds         SO3_GEODEFSTR* (i)   tensors at curr. GP
\param   stress      DOUBLE[]       (o)   stress
\param   cmat        DOUBLE[][]     (o)   elasticity tensor
\return  void
\author bborn
\date 04/07
*/
void so3_mat_robinson_sel(const CONTAINER* container,
                          const ELEMENT* ele,
                          const VP_ROBINSON* mat_robin,
                          const INT ip,
                          const SO3_GEODEFSTR* gds,
                          DOUBLE stress[NUMSTR_SOLID3],
                          DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_sel");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish time integration scheme */
  if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
  {
    so3_mat_robinson_fb4_sel(container,
                             ele,
                             mat_robin,
                             ip,
                             gds,
                             stress,
                             cmat);
  }
  else if (tsidyn->kind == tsi_therm_stat_struct_genalp)
  {
    so3_mat_robinson_be_sel(container,
                            ele,
                            mat_robin,
                            ip,
                            gds,
                            stress,
                            cmat);
  }
  else
  {
    dserror("TSI time integration is not available for Robinson's material");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update Robinson's internal material variables 
\param   ele          ELEMENT*        (io)   curr. elem.
\parm    mat_robin    VP_ROBINSON*    ()     elem. mater.
\author bborn
\date 04/07
*/
void so3_mat_robinson_mivupd(ELEMENT* ele,
                             VP_ROBINSON* mat_robin)
{
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_mivupd");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish time integration scheme */
  if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
  {
    so3_mat_robinson_fb4_mivupd(ele, mat_robin);
  }
  else if (tsidyn->kind == tsi_therm_stat_struct_genalp)
  {
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Select Robinson's material and return stress
\param   container   CONTAINER*     (i)   container.h data
\param   ele         ELEMENT*       (i)   curr. element
\param   mat_robin   VP_ROBINSON*   (i)   curr. Robinson material data
\param   ip          INT            (i)   total GP index
\param   gds         SO3_GEODEFSTR* (i)   tensors at curr. GP
\param   stress      DOUBLE[]       (o)   stress
\param   cmat        DOUBLE[][]     (o)   elasticity tensor
\return  void
\author bborn
\date 04/07
*/
void so3_mat_robinson_stress(const CONTAINER* container,
                             const ELEMENT* ele,
                             const VP_ROBINSON* mat_robin,
                             const INT ip,
                             const SO3_GEODEFSTR* gds,
                             DOUBLE stress[NUMSTR_SOLID3],
                             DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* distinguish time integration scheme */
  if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
  {
    so3_mat_robinson_fb4_stress(container,
                                ele,
                                mat_robin,
                                ip,
                                gds,
                                stress,
                                cmat);
  }
  else if (tsidyn->kind == tsi_therm_stat_struct_genalp)
  {
    so3_mat_robinson_be_stress(container,
                               ele,
                               mat_robin,
                               ip,
                               gds,
                               stress,
                               cmat);
  }
  else
  {
    dserror("Cannot use Robinson's material"
            " with selected TSI time integration");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif  /* end D_TSI */
#endif  /* end D_SOLID3 */
