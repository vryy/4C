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
\brief Allocate material internal variables (MIVs)
\author bborn
\date 04/07
 */
void so3_mat_robinson_fb4_init(SOLID3* actso3)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_init");
#endif

  /*--------------------------------------------------------------------*/
  if (actso3->miv_rob == NULL)
  {
    /* allocation of MIVar */
    actso3->miv_rob 
      = (SO3_MIV_ROBINSON*) CCACALLOC(1, sizeof(SO3_MIV_ROBINSON));
    /* allocation of arrays in MIV */
    amdef("miv_rob_vicstn", &(actso3->miv_rob->vicstn),
          actso3->gptot, NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->vicstn));
    amdef("miv_rob_bacsts", &(actso3->miv_rob->bacsts),
          actso3->gptot, NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->bacsts));
    /* Fehlberg4 staged data --- OK to assume 'coz checked above */
    INT fb4_stg = tsi_fehlbg4_stages();
    am4def("miv_rob_dvicstn", &(actso3->miv_rob->dvicstn),
           actso3->gptot, fb4_stg, NUMSTR_SOLID3, 0, "D3");
    am4zero(&(actso3->miv_rob->dvicstn));
    am4def("miv_rob_dbacsts", &(actso3->miv_rob->dbacsts),
           actso3->gptot, fb4_stg, NUMSTR_SOLID3, 0, "D3");
    am4zero(&(actso3->miv_rob->dbacsts));
  }
  else
  {
    dserror("Material internal variables already allocated!");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Deallocate material internal variables (MIVs)
\author bborn
\date 04/07
 */
void so3_mat_robinson_fb4_final(SOLID3* actso3)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_final");
#endif

  /*--------------------------------------------------------------------*/
  if (actso3->miv_rob != NULL)
  {
    amdel(&(actso3->miv_rob->vicstn));
    amdel(&(actso3->miv_rob->bacsts));
    am4del(&(actso3->miv_rob->dvicstn));
    am4del(&(actso3->miv_rob->dbacsts));
    CCAFREE(actso3->miv_rob);
  }
  else
  {
    dserror("Material internal variables are not allocated!");
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
void so3_mat_robinson_fb4_sel(const CONTAINER* container,
                              const ELEMENT* ele,
                              const VP_ROBINSON* mat_robin,
                              const INT ip,
                              const SO3_GEODEFSTR* gds,
                              DOUBLE stress[NUMSTR_SOLID3],
                              DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  const SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  DOUBLE stntot[NUMSTR_SOLID3];  /* total strain */
  DOUBLE stnela[NUMSTR_SOLID3];  /* elastic strain */
  DOUBLE stnthr[NUMSTR_SOLID3];  /* thermal strain */
  DOUBLE stnvsc[NUMSTR_SOLID3];  /* viscous strain */
  
  DOUBLE stsbac[NUMSTR_SOLID3];  /* back stress */
  DOUBLE stsdev[NUMSTR_SOLID3];  /* stress deviator */
  DOUBLE stsovr[NUMSTR_SOLID3];  /* overstress */

  DOUBLE tem;  /* temperature */

  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
  const INT numstg = tsi_fehlbg4.stg;  /* number of RK stages */
  const INT actstg = tsidyn->actstg;  /* curr. RK stage */
  const DOUBLE dt = tsidyn->dt;  /* time step size */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_sel");
#endif

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
  /* temperature and thermal strain */
#ifdef D_THERM3
  {
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
    /* initialise viscous strain at t_{n+c_i} */
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
    DOUBLE emod;
    so3_mat_robinson_prmbytmpr(mat_robin->youngmodul_ipl,
                               mat_robin->youngmodul_n,
                               mat_robin->youngmodul,
                               tem, &emod);
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
  /* deviatoric stress 's' at t_{n+c_i} */
  so3_vct6_dev(stress, stsdev);

  /*--------------------------------------------------------------------*/
  /* overstress 'Sig' = 's' - 'alpha' */
  so3_vct6_sub(stsdev, stsbac, stsovr);

  /*--------------------------------------------------------------------*/
  /* viscous strain rate at t_{n+c_i} */
  so3_mat_robinson_fb4_stnvscrat(mat_robin, tem, stsdev, stsovr,
                                 actso3->miv_rob->dvicstn.a.d3[ip][actstg]);

  /*--------------------------------------------------------------------*/
  /* back stress rate */
  so3_mat_robinson_fb4_stsbckrat(mat_robin, tem, stsdev, stsbac,
                                 actso3->miv_rob->dvicstn.a.d3[ip][actstg],
                                 actso3->miv_rob->dbacsts.a.d3[ip][actstg]);

  /*--------------------------------------------------------------------*/
  /* update internal variables if we are in last RK stage */
  if (actstg == (numstg - 1) )
  {
    so3_mat_robinson_fb4_mivupd(ele, mat_robin, ip);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Viscous strain rate at t_{n+c_i} at Gauss point
\param   mat_robin    VP_ROBINSON*  (i)   element Robinson material
\param   tmpr         DOUBLE        (i)   temperature
\param   stsdev       DOUBLE[]      (i)   stress deviator
\param   stsovr       DOUBLE[]      (i)   over stress
\param   dstnvsc      DOUBLE[]      (o)   viscous strain rate
\author bborn
\date 04/07
*/
void so3_mat_robinson_fb4_stnvscrat(const VP_ROBINSON* mat_robin,
                                    const DOUBLE tmpr,
                                    const DOUBLE stsdev[NUMSTR_SOLID3],
                                    const DOUBLE stsovr[NUMSTR_SOLID3],
                                    DOUBLE dstnvsc[NUMSTR_SOLID3])
{
  DOUBLE j2;  /* 'J_2' */
  DOUBLE shrthrshld;  /* shear threshold 'K^2' */
  DOUBLE ff;  /* 'F' */
  DOUBLE ss;  /* */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_stnvscrat");
#endif

  /*--------------------------------------------------------------------*/
  /* preliminaries */
  /* J_2 = 1/2 * Sig : Sig  with Sig...overstress */
  so3_vct6_dblctr(stsovr, stsovr, &j2);
  j2 *= 0.5;
  /* Bingham-Prager shear stress threshold at current temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld_ipl,
                             mat_robin->shrthrshld_n,
                             mat_robin->shrthrshld,
                             tmpr, &shrthrshld);
  /* F = (J_2 - K^2)/K_2 */
  if (fabs(shrthrshld) <= EPS10)
  {
    ff = -1.0;
    dserror("Division by zero: Shear threshold very close to zero");
  }
  else
  {
    /* debug: */
    /* printf("J_2 %f;  K^2 %f\n", j2, shrthrshld); */
    ff = (j2 - shrthrshld)/shrthrshld;
  }
  /* ss = 1/2 * s : Sig  with  Sig...overstress, s...deviat.stress */
  so3_vct6_dblctr(stsovr, stsdev, &ss);
  ss *= 0.5;

  /*--------------------------------------------------------------------*/
  /* viscous strain rate at t_{n+c_i} */
  if ( (ff > 0.0) && (ss > 0.0) )
  {
    DOUBLE fct;
    if (mat_robin->kind == 1)  /* Butler */
    {
      fct = mat_robin->hrdn_fact * pow(ff, mat_robin->hrdn_expo) / sqrt(j2);
    }
    else if (mat_robin->kind == 2)  /* Arya */
    {
      fct = mat_robin->hrdn_fact * pow(ff, mat_robin->hrdn_expo) / sqrt(j2);
    }
    else
    {
      fct = -1.0;
      dserror("Kind of Robinson material is unknown");
    }
    so3_vct6_2_assscl(fct, stsovr, dstnvsc);
  }
  else
  {
    so3_vct6_zero(dstnvsc);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Back stress rate at t_{n+c_i} at Gauss point
\param   mat_robin    VP_ROBINSON*  (i)   element Robinson material
\param   tmpr         DOUBLE        (i)   temperature
\param   stsdev       DOUBLE[]      (i)   deviatoric stress
\param   stsbck       DOUBLE[]      (i)   back stress
\param   dstnvsc      DOUBLE[]      (i)   viscous strain rate
\param   dstsbck      DOUBLE[]      (o)   back stress rate
\author bborn
\date 04/07
*/
void so3_mat_robinson_fb4_stsbckrat(const VP_ROBINSON* mat_robin,
                                    const DOUBLE tmpr,
                                    const DOUBLE stsdev[NUMSTR_SOLID3],
                                    const DOUBLE stsbck[NUMSTR_SOLID3],
                                    const DOUBLE dstnvsc[NUMSTR_SOLID3],
                                    DOUBLE dstsbck[NUMSTR_SOLID3])
{
  DOUBLE i2;  /* 'I_2' */
  DOUBLE tem0;  /* activation temperature */
  DOUBLE shrthrshld0;  /* shear threshold 'K_0^2' */
  DOUBLE hh;  /* 'H' */
  DOUBLE beta;  /* 'beta' */
  DOUBLE mm;  /* 'm' */
  DOUBLE rr0;  /* recovery factor 'R_0' */
  DOUBLE rr;  /* recovery term 'R' */
  DOUBLE gg0;  /* 'G_0' */
  DOUBLE gg;  /* 'G' */
  DOUBLE sa;  /* */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_stsbckrat");
#endif

  /*--------------------------------------------------------------------*/
  /* preliminaries */
  /* I_2 = 1/2 * Alpha : Alpha  with Alpha...back stress */
  so3_vct6_dblctr(stsbck, stsbck, &i2);
  i2 *= 0.5;
  /* debug: */
/*   printf("I_2 %f\n", sqrt(i2)); */
/*   printf("back stress %f %f %f %f %f %f\n", stsbck[0], stsbck[1], stsbck[2], stsbck[3], stsbck[4], stsbck[5]); */
  /* activation temperature */
  tem0 = mat_robin->actv_tmpr;
  /* Bingham-Prager shear stress threshold at activation temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld_ipl,
                             mat_robin->shrthrshld_n,
                             mat_robin->shrthrshld,
                             tem0,
                             &shrthrshld0);
  /* 'H' at current temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->h_ipl,
                             mat_robin->h_n,
                             mat_robin->h,
                             tmpr, 
                             &hh);
  /* 'beta' at current temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->beta_ipl,
                             mat_robin->beta_n,
                             mat_robin->beta,
                             tmpr, 
                             &beta);
  /* 'm' */
  mm = mat_robin->m;
  /* recovery factor 'R_0' */
  so3_mat_robinson_prmbytmpr(mat_robin->rcvry_ipl,
                             mat_robin->rcvry_n,
                             mat_robin->rcvry,
                             tmpr, 
                             &rr0);
  /* 'R' */
  rr = rr0 * exp(mat_robin->actv_ergy*(tmpr-tem0)/(tmpr*tem0));
  /* 'G_0' */
  gg0 = mat_robin->g0;
  /* G = I_2/K_0^2 */
  if (fabs(shrthrshld0) <= EPS10)
  {
    dserror("Division by zero: Shear threshold very close to zero");
  }
  else
  {
    if (mat_robin->kind == 1)
    {
      gg = i2/shrthrshld0;
    }
    else if (mat_robin->kind == 2)
    {
      gg = sqrt(i2/shrthrshld0);
    }
    else
    {
      gg = -1.0;
    }
  }
  /* ss = 1/2 * s : Alpha  with  Alpha...backstress, s...deviat.stress */
  so3_vct6_dblctr(stsbck, stsdev, &sa);
  sa *= 0.5;

  /*--------------------------------------------------------------------*/
  /* back stress rate at t_{n+c_i} */
  if ( (gg > gg0) && (sa > 0.0) )
  {
    DOUBLE fcte;
    if (mat_robin->kind == 1)
    {
      fcte = hh / pow(gg, beta);
    }
    else if (mat_robin->kind == 2)
    {
      fcte = hh / pow(gg, beta);
    }
    else
    {
      fcte = -1.0;
    }
    so3_vct6_05_assscl(fcte, dstnvsc, dstsbck);
    DOUBLE fcta;
    if (mat_robin->kind == 1)
    {
      fcta = -rr * pow(gg, (mm-beta)) / sqrt(i2);
    }
    else if (mat_robin->kind == 2)
    {
      fcta = -rr * pow(gg, (mm-beta)) / sqrt(i2);
    }
    else
    {
      fcta = -1.0;
    }
    so3_vct6_updscl(fcta, stsbck, dstsbck);
  }
  else
  {
    DOUBLE fcte;
    if (mat_robin->kind == 1)
    {
      fcte = hh / pow(gg, beta);
    }
    else if (mat_robin->kind == 2)
    {
      fcte = hh / pow(gg0, beta);
    }
    else
    {
      fcte = -1.0;
    }
    so3_vct6_05_assscl(fcte, dstnvsc, dstsbck);
    DOUBLE fcta;
    if (mat_robin->kind == 1)
    {
      fcta = -rr * pow(gg0, (mm-beta)) / sqrt(i2);
    }
    else if (mat_robin->kind == 2)
    {
      if (sqrt(i2) < EPS10)
      {
        DOUBLE ii2 = 1.e6;
        fcta = -rr * pow(gg0, (mm-beta)) / ii2;
      }
      else
      {
        fcta = -rr * pow(gg0, (mm-beta)) / sqrt(i2);
      }
    }
    else
    {
      fcta = -1.0;
    }
    so3_vct6_updscl(fcta, stsbck, dstsbck);
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
\parm    mat_robin    VP_ROBINSON*    (i)    elem. mater.
\param   ip           INT             (i)    curr. Gauss point index
\author bborn
\date 04/07
*/
void so3_mat_robinson_fb4_mivupd(ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const INT ip)
{
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  INT gptot = actso3->gptot;  /* total number of GPs in domain */
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
  const DOUBLE dt = tsidyn->dt;  /* time step size */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_mivupd");
#endif

  /*--------------------------------------------------------------------*/
  /* HINT: The precondition if we are in a Fehlberg4 time integration is
   *       not explictly tested here. It has to be done before. */

  /*--------------------------------------------------------------------*/
  /* update viscous strain */
  {
    DOUBLE stnvscn[NUMSTR_SOLID3];
    so3_vct6_ass(actso3->miv_rob->vicstn.a.da[ip], stnvscn);
    INT istg;
    for (istg=0; istg<tsi_fehlbg4.stg; istg++)
    {
      DOUBLE fct = dt * tsi_fehlbg4.b[istg];
      so3_vct6_updscl(fct, actso3->miv_rob->dvicstn.a.d3[ip][istg], 
                      stnvscn);
    }
    so3_vct6_ass(stnvscn, actso3->miv_rob->vicstn.a.da[ip]);
  }

  /*--------------------------------------------------------------------*/
  /* update back stress 'Alpha' */
  {
    DOUBLE stsbckn[NUMSTR_SOLID3];
    so3_vct6_ass(actso3->miv_rob->bacsts.a.da[ip], stsbckn);
    INT istg;
    for (istg=0; istg<tsi_fehlbg4.stg; istg++)
    {
      DOUBLE fct = dt * tsi_fehlbg4.b[istg];
      so3_vct6_updscl(fct, actso3->miv_rob->dbacsts.a.d3[ip][istg], 
                      stsbckn);
    }
    so3_vct6_ass(stsbckn, actso3->miv_rob->bacsts.a.da[ip]);
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
void so3_mat_robinson_fb4_stress(const CONTAINER* container,
                                 const ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const INT ip,
                                 const SO3_GEODEFSTR* gds,
                                 DOUBLE stress[NUMSTR_SOLID3],
                                 DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  DOUBLE stntot[NUMSTR_SOLID3];  /* total strain */
  DOUBLE stnela[NUMSTR_SOLID3];  /* elastic strain */
  DOUBLE stnthr[NUMSTR_SOLID3];  /* thermal strain */
  DOUBLE stnvsc[NUMSTR_SOLID3];  /* viscous strain */
  
  DOUBLE tem;  /* temperature */

  INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_fb4_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* total strain at t_{n} */
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
  /* temperature and thermal strain */
#ifdef D_THERM3
  {
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
  /* viscous strain at t_{n} */
  {
    /* set viscous strain at t_{n} */
    INT istn;
    for (istn=0; istn<NUMSTR_SOLID3; istn++)
    {
      stnvsc[istn] = actso3->miv_rob->vicstn.a.da[ip][istn];
    }
  }

  /*--------------------------------------------------------------------*/
  /* elastic strain at t_{n} */
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
    DOUBLE emod;  /* modulus of elasticity */
    so3_mat_robinson_prmbytmpr(mat_robin->youngmodul_ipl,
                               mat_robin->youngmodul_n,
                               mat_robin->youngmodul,
                               tem, &emod);
    DOUBLE prat = mat_robin->possionratio;  /* Poisson's ratio */
    DOUBLE mfac = emod/((1.0+prat)*(1.0-2.0*prat));  /* factor */
    /* zero out content */
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
  /* stress at t_{n} */
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
  /* debug: */
  if (ip == 0)
  {
    DOUBLE stsdev[NUMSTR_SOLID3];
    /* stress deviator */
    so3_vct6_dev(stress, stsdev);
    /* overstress 'Sig' = 's' - 'alpha' */
    DOUBLE stsovr[NUMSTR_SOLID3];
    so3_vct6_sub(stsdev, actso3->miv_rob->bacsts.a.da[ip], stsovr);
    FILE* oo = fopen("so3_mat_robin_ip0.dat", "a+");
    fprintf(oo, "# Time step %d Gauss point %d\n", container->kstep, ip);
    fprintf(oo, "# step strain-tot strain-ela strain-vis strain-thr"
            "  stress-tot stress-dev stress-back stress-over\n");
    INT i;
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      fprintf(oo, "%d %f %f %f %f  %f %f %f %f\n", 
              container->kstep,
              stntot[i], stnela[i], stnvsc[i], stnthr[i],
              stress[i], 
              stsdev[i], 
              actso3->miv_rob->bacsts.a.da[ip][i], 
              stsovr[i]);
    }
    /* fprintf(oo, "\n"); */
    fclose(oo);
  }
  if (ip == 7)
  {
    DOUBLE stsdev[NUMSTR_SOLID3];
    /* stress deviator */
    so3_vct6_dev(stress, stsdev);
    /* overstress 'Sig' = 's' - 'alpha' */
    DOUBLE stsovr[NUMSTR_SOLID3];
    so3_vct6_sub(stsdev, actso3->miv_rob->bacsts.a.da[ip], stsovr);
    FILE* oo = fopen("so3_mat_robin_ip7.dat", "a+");
    fprintf(oo, "# Time step %d Gauss point %d\n", container->kstep, ip);
    fprintf(oo, "# step strain-tot strain-ela strain-vis strain-thr"
            "  stress-tot stress-dev stress-back stress-over\n");
    INT i;
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      fprintf(oo, "%d %f %f %f %f  %f %f %f %f\n", 
              container->kstep,
              stntot[i], stnela[i], stnvsc[i], stnthr[i],
              stress[i], stsdev[i], 
              actso3->miv_rob->bacsts.a.da[ip][i], stsovr[i]);
    }
    /* fprintf(oo, "\n"); */
    fclose(oo);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif  /* end D_TSI */
#endif  /* end D_SOLID3 */
