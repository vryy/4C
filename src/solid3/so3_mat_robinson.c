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

/*---------------------------------------------------------------------*/
/*!
\brief Fehlberg4 parameters

\author bborn
\date 03/07
*/
extern TSI_FEHLBG4 tsi_fehlbg4;

/*======================================================================*/
/*!
\brief Verify context of Robinson's material, ie must be
       - thermo-structure interaction
       - Fehlberg4 time integration
\author bborn
\date 04/07
*/
void so3_mat_robinson_ctxtvrfy()
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_ctxtvrfy");
#endif

  /*--------------------------------------------------------------------*/
  if (genprob.probtyp != prb_tsi)
  {
    dserror("Trying to use Robinson's material in non-TSI calculation!");
  }

  /*--------------------------------------------------------------------*/
  {
    INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
    TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
    if (tsidyn->kind != tsi_therm_stat_struct_fehlbg)
    {
      dserror("Robinson's material only applicable in Fehlberg4 TSI");
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

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
  /* check if we are in a sane context */
  so3_mat_robinson_ctxtvrfy();

  /*--------------------------------------------------------------------*/
  /* allocate MIV */
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
      /* Fehlberg4 staged data --- OK to assume 'coz checked above */
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
  /* check if we are in a sane context */
  so3_mat_robinson_ctxtvrfy();

  /*--------------------------------------------------------------------*/
  /* deallocate MIV */
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

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end so3_mat_robinson_init */

/*======================================================================*/
/*!
\brief Get temperature-dependent material parameter at current temperature
\author bborn
\date 04/07
*/
void so3_mat_robinson_prmbytmpr(VP_ROBINSON_INTPOL ipl,  /*!< interpolation 
                                                           type */
                                INT prm_n,  /*!< magnitude param data */
                                DOUBLE* prm,  /*!< parameter data */
                                DOUBLE tmpr,  /*!< curr. temperature */
                                DOUBLE* prmbytempr)  /*!< param. at current
                                                      temperature */
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_prmbytmpr");
#endif

  /*--------------------------------------------------------------------*/
  /* constant */
  if (ipl == vp_robinson_ipl_const)
  {
    *prmbytempr = prm[0];
  }
  /* polynomial */
  else if (ipl == vp_robinson_ipl_poly)
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
  else if (ipl == vp_robinson_ipl_pcwslnr)
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
  else if (ipl == vp_robinson_ipl_none)
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
\brief Select Robinson's material
\author bborn
\date 04/07
*/
void so3_mat_robinson_sel(CONTAINER* container,  /*!< container */
                          ELEMENT* ele,  /*!< curr. elem. */
                          VP_ROBINSON* mat_robin,  /*!< elem. mater. */
                          INT ip,  /*!< total Gauss pnt. index */
                          SO3_GEODEFSTR* gds,  /*!< elem. data at Gauss 
                                                 point */
                          DOUBLE stress[NUMSTR_SOLID3],  /*!< stress */
                          DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])  /*!< elasticity tensor */
{
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  DOUBLE stntot[NUMSTR_SOLID3];  /* total strain */
  DOUBLE stnela[NUMSTR_SOLID3];  /* elastic strain */
  DOUBLE stnthr[NUMSTR_SOLID3];  /* thermal strain */
  DOUBLE stnvsc[NUMSTR_SOLID3];  /* viscous strain */
  
  DOUBLE stsbac[NUMSTR_SOLID3];  /* back stress */
  DOUBLE stsdev[NUMSTR_SOLID3];  /* stress deviator */
  DOUBLE stsovr[NUMSTR_SOLID3];  /* overstress */

  INT fb4_stg = tsi_fehlbg4_stages();
  INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
  INT actstg = tsidyn->actstg;  /* curr. RK stage */
  DOUBLE dt = tsidyn->dt;  /* time step size */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_sel");
#endif

  /*--------------------------------------------------------------------*/
  /* check if we are in a sane context */
  so3_mat_robinson_ctxtvrfy();

  /*--------------------------------------------------------------------*/
  /* total strain at t_{n+c_i} */
  if (actso3->kintype == so3_geo_lin)
  {
    /* linear (engineering) strain vector */
    INT istrn;
    for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
    {
      stntot[istrn] = gds->stnengv[istrn];
    }
  } 
  else if (actso3->kintype == so3_total_lagr)
  {
    /* Green-Lagrange strain vector */
    INT istrn;
    for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
    {
      stntot[istrn] = gds->stnglv[istrn];
    }
  }
  else
  {
    dserror("Cannot digest chosen type of spatial kinematic\n");
  }

  /*--------------------------------------------------------------------*/
  /* thermal strain */
#ifdef D_THERM3
  {
    DOUBLE tem;  /* temperature */
    /* temperature at Gauss point */
    so3_tsi_temper(container, ele,
                   gds->gpc[0], gds->gpc[1], gds->gpc[2], 
                   &tem);
    /* coefficient of linear thermal expansion */
    DOUBLE thermexpans = mat_robin->thermexpans;
    /* thermal strain vector */
    stnthr[0] = thermexpans * tem;  /* E_xx */
    stnthr[1] = thermexpans * tem;  /* E_yy */
    stnthr[2] = thermexpans * tem;  /* E_zz */
    stnthr[3] = 0.0;  /* E_xy */
    stnthr[4] = 0.0;  /* E_yz */
    stnthr[5] = 0.0;  /* E_zx */
  }
#endif

  /*--------------------------------------------------------------------*/
  /* viscous strain at t_{n+c_i} */
  {
    /* initialise viscous strain at t_{n+c+i} */
    INT istn;
    for (istn=0; istn<NUMSTR_SOLID3; istn++)
    {
      stnvsc[istn] = actso3->miv_rob->vicstn.a.da[ip][istn];
    }
    /* create value (this is Fehlberg4) */
    INT istg;
    for (istg=0; istg<actstg; istg++)
    {
      for (istn=0; istn<NUMSTR_SOLID3; istn++)
      {
        stnvsc[istn] += dt * tsi_fehlbg4.a[actstg][istn]
          * actso3->miv_rob->dvicstn.a.d3[ip][istg][istn];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* elastic strain at t_{n+c_i} */
  {
    INT istr;  /* strain index */
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      stnela[istr] = stntot[istr] - stnvsc[istr] - stnthr[istr];
    }
  }

  /*--------------------------------------------------------------------*/
  /* elasticity tensor */
  {
    DOUBLE emod = mat_robin->youngmodul;
    DOUBLE prat = mat_robin->possionratio;
    DOUBLE mfac = emod/((1.0+prat)*(1.0-2.0*prat));  /* factor */
    /* zero content */
    INT istr, jstr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      for (jstr=0; jstr<NUMSTR_SOLID3; jstr++)
      {
        cmat[istr][jstr] = 0.0;
      }
    }
    /* non-zero content --- axial */
    cmat[0][0] = mfac*(1.0-prat);
    cmat[0][1] = mfac*prat;
    cmat[0][2] = mfac*prat;
    cmat[1][0] = mfac*prat;
    cmat[1][1] = mfac*(1.0-prat);
    cmat[1][2] = mfac*prat;
    cmat[2][0] = mfac*prat;
    cmat[2][1] = mfac*prat;
    cmat[2][2] = mfac*(1.0-prat);
    /* non-zero content --- shear */
    cmat[3][3] = mfac*0.5*(1.0-2.0*prat);
    cmat[4][4] = mfac*0.5*(1.0-2.0*prat);
    cmat[5][5] = mfac*0.5*(1.0-2.0*prat);
  }

  /*--------------------------------------------------------------------*/
  /* stress at t_{n+c_i} */
  {
    INT istss;  /* stress index */
    for (istss=0; istss<NUMSTR_SOLID3; istss++)
    {
      INT istrn;  /* strain index */
      DOUBLE stresssum = 0.0;  /* intermediate row * column sum */
      for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
      {
        stresssum += cmat[istss][istrn] * stnela[istrn];
      }
      stress[istss] = stresssum;
    }
  }

  /*--------------------------------------------------------------------*/
  /* back stress at t_{n+c_i} */
  {
    /* initialise viscous strain at t_{n+c+i} */
    INT istn;
    for (istn=0; istn<NUMSTR_SOLID3; istn++)
    {
      stsbac[istn] = actso3->miv_rob->bacsts.a.da[ip][istn];
    }
    /* create value (this is Fehlberg4) */
    INT istg;
    for (istg=0; istg<actstg; istg++)
    {
      for (istn=0; istn<NUMSTR_SOLID3; istn++)
      {
        stsbac[istn] += dt * tsi_fehlbg4.a[actstg][istn]
          * actso3->miv_rob->dbacsts.a.d3[ip][istg][istn];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* deviatoric stress at t_{n+c_i} */
  so3_vct6_dev(stress, stsdev);

  /*--------------------------------------------------------------------*/
  /* overstress */
  so3_vct6_sub(stsdev, stsbac, stsovr);

  /*--------------------------------------------------------------------*/
  /* viscous strain rate at t_{n+c_i} */


  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Select Robinson's material
\author bborn
\date 04/07
*/
void so3_mat_robinson_mivupd(CONTAINER* container,  /*!< container */
                             ELEMENT* ele,  /*!< curr. elem. */
                             MATERIAL* mat,  /*!< elem. mater. */
                             INT ip,  /*!< total Gauss pnt. index */
                             SO3_GEODEFSTR* gds,  /*!< elem. data at Gauss 
                                                 point */
                          DOUBLE stress[NUMSTR_SOLID3],  /*!< stress */
                          DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])  /*!< elasticity tensor */
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_mivupd");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif  /* end D_TSI */
#endif  /* end D_SOLID3 */
