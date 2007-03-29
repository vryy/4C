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
#include "../tsi_full/tsi_fehlbg4.h"

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

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_init");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate MIV */
  /* check if we are in a thermo-structure interaction problem */
  if (genprob.probtyp == prb_tsi)
  {
    INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
    TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
    if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
    {
      SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
      /* allocation of MIVar */
      if (actso3->miv_rob == NULL)
      {
        /* allocation of MIVar */
        actso3->miv_rob 
          = (SO3_MIV_ROBINSON*) CCACALLOC(1, sizeof(SO3_MIV_ROBINSON));
        /* allocation of arrays in MIV */
        amdef("miv_rob_vicstn", &(actso3->miv_rob->vicstn),
              actso3->gptot, NUMSTR_SOLID3, "DA");
        amdef("miv_rob_vicstnn", &(actso3->miv_rob->vicstnn),
              actso3->gptot, NUMSTR_SOLID3, "DA");
        amdef("miv_rob_bacsts", &(actso3->miv_rob->bacsts),
              actso3->gptot, NUMSTR_SOLID3, "DA");
        amdef("miv_rob_bacstsn", &(actso3->miv_rob->bacstsn),
              actso3->gptot, NUMSTR_SOLID3, "DA");
        INT fb4_stg = tsi_fehlbg4_stages();
        am4def("miv_rob_dvicstn", &(actso3->miv_rob->dvicstn),
               actso3->gptot, fb4_stg, NUMSTR_SOLID3, 0, "D3");
        am4def("miv_rob_dbacsts", &(actso3->miv_rob->dbacsts),
               actso3->gptot, fb4_stg, NUMSTR_SOLID3, 0, "D3");
      }
      else
      {
        dserror("Material internal variables already allocated!");
      }
    }
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

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocate MIV */
  /* check if we are in a thermo-structure interaction problem */
  if (genprob.probtyp == prb_tsi)
  {
    INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
    TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
    if (tsidyn->kind == tsi_therm_stat_struct_fehlbg)
    {
      SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
      /* deallocation of MIVar */
      if (actso3->miv_rob != NULL)
      {
        amdel(&(actso3->miv_rob->vicstn));
        amdel(&(actso3->miv_rob->vicstnn));
        amdel(&(actso3->miv_rob->bacsts));
        amdel(&(actso3->miv_rob->bacstsn));
        am4del(&(actso3->miv_rob->dvicstn));
        am4del(&(actso3->miv_rob->dbacsts));
        CCAFREE(actso3->miv_rob);
      }
      else
      {
        dserror("Material internal variables are not allocated!");
      }
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end so3_mat_robinson_init */

#endif  /* end D_TSI */
#endif  /* end D_SOLID3 */
