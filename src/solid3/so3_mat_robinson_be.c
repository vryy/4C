/*======================================================================*/
/*!
\file
\brief Robinson's visco-plastic material
       with Backward Euler time integration

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
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 05/07
*/
extern FIELD *field;

/*----------------------------------------------------------------------*/
/*!
\brief Alldyn dynamic control

pointer to allocate dynamic variables if needed
dedfined in global_control.c

\auther bborn
\date 03/06
*/
extern ALLDYNA* alldyn;

/*----------------------------------------------------------------------*/
/*!
\brief IO files
\author bborn
\date 05/07
*/
#ifdef TESTROBIN_SOLID3
extern FILES allfiles;
#endif

/*======================================================================*/
/*!
\brief Allocate material internal variables (MIVs)
\author bborn
\date 04/07
 */
void so3_mat_robinson_be_init(SOLID3* actso3)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_init");
#endif

  /*--------------------------------------------------------------------*/
  /* allocation of MIVar */
  if (actso3->miv_rob == NULL)
  {
    /* allocation of MIVar */
    actso3->miv_rob 
      = (SO3_MIV_ROBINSON*) CCACALLOC(1, sizeof(SO3_MIV_ROBINSON));
    /* allocation of arrays in MIV */
    amdef("miv_rob_vscstns", &(actso3->miv_rob->vscstns),
          actso3->gptot, 1, "IV");
    amzero(&(actso3->miv_rob->vscstns));
    amdef("miv_rob_vicstn", &(actso3->miv_rob->vicstn),
          actso3->gptot, NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->vicstn));
    amdef("miv_rob_vicstnn", &(actso3->miv_rob->vicstnn),
          actso3->gptot, NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->vicstn));
    amdef("miv_rob_bckstss", &(actso3->miv_rob->bckstss),
          actso3->gptot, 1, "IV");
    amzero(&(actso3->miv_rob->bckstss));
    amdef("miv_rob_bacsts", &(actso3->miv_rob->bacsts),
          actso3->gptot, NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->bacsts));
    amdef("miv_rob_bacstsn", &(actso3->miv_rob->bacstsn),
          actso3->gptot, NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->bacsts));
    amdef("miv_rob_kvarva", &(actso3->miv_rob->kvarva),
          actso3->gptot, 2*NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->kvarva));
    amdef("miv_rob_kvakvae", &(actso3->miv_rob->kvakvae),
          actso3->gptot, 2*NUMSTR_SOLID3 * NUMSTR_SOLID3, "DA");
    amzero(&(actso3->miv_rob->kvakvae));
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
void so3_mat_robinson_be_final(SOLID3* actso3)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocation of MIVar */
  if (actso3->miv_rob != NULL)
  {
    amdel(&(actso3->miv_rob->vscstns));
    amdel(&(actso3->miv_rob->vicstn));
    amdel(&(actso3->miv_rob->vicstnn));
    amdel(&(actso3->miv_rob->bckstss));
    amdel(&(actso3->miv_rob->bacsts));
    amdel(&(actso3->miv_rob->bacstsn));
    amdel(&(actso3->miv_rob->kvarva));
    amdel(&(actso3->miv_rob->kvakvae));
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
\param   ele         ELEMENT*       (io)  curr. element
\param   mat_robin   VP_ROBINSON*   (i)   curr. Robinson material data
\param   ip          INT            (i)   total GP index
\param   gds         SO3_GEODEFSTR* (i)   tensors at curr. GP
\param   stress      DOUBLE[]       (o)   stress
\param   cmat        DOUBLE[][]     (o)   elasticity tensor
\return  void
\author bborn
\date 04/07
*/
void so3_mat_robinson_be_sel(const CONTAINER* container,
                             ELEMENT* ele,
                             const VP_ROBINSON* mat_robin,
                             SO3_GPSHAPEDERIV* gpshade,
                             const INT ip,
                             SO3_GEODEFSTR* gds,
                             DOUBLE stress[NUMSTR_SOLID3],
                             DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  const INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  const TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */
  const DOUBLE dt = tsidyn->dt;  /* time step size */

  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  DOUBLE stntotn[NUMSTR_SOLID3];  /* total strain */
  DOUBLE stnela[NUMSTR_SOLID3];  /* elastic strain */
  DOUBLE stnthr[NUMSTR_SOLID3];  /* thermal strain */
  DOUBLE stnvscn[NUMSTR_SOLID3];  /* viscous strain */
  
  DOUBLE devstsn[NUMSTR_SOLID3];  /* stress deviator */
  DOUBLE stsovrn[NUMSTR_SOLID3];  /* over stress */

  DOUBLE vscstnr[NUMSTR_SOLID3];  /* visc. strain residual */
  DOUBLE bckstsr[NUMSTR_SOLID3];  /* back stress residual */

  DOUBLE kev[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd sig}{\pd eps^v} */
  DOUBLE kea[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd sig}{\pd al} */
  DOUBLE kve[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd res^v}{\pd eps} */
  DOUBLE kvv[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd res^v}{\pd eps^v} */
  DOUBLE kva[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd res^v}{\pd al} */
  DOUBLE kae[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd res^al}{\pd eps} */
  DOUBLE kav[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd res^al}{\pd eps^v} */
  DOUBLE kaa[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* \frac{\pd res^al}{\pd al} */

  DOUBLE tem;  /* temperature */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_sel");
#endif

  /*--------------------------------------------------------------------*/
  /* Material internal variables (visc. strain and back stress) are
   * updated by their iterative increments.
   * Their iterative increments are expressed in terms of the iterative
   * increment of the total strain.
   * Here the reduction matrices --- stored in the previous iteration ---
   * are used. */
  if ( (actso3->miv_rob->vscstns.a.iv[ip] != so3_mat_robinson_state_vague)
       && (actso3->miv_rob->bckstss.a.iv[ip] != so3_mat_robinson_state_vague) )
  {
    so3_mat_robinson_be_mivupditr(container, ele, 
                                  gpshade->gpderiv[ip],
                                  gds,
                                  actso3->miv_rob->kvarva.a.da[ip],
                                  actso3->miv_rob->kvakvae.a.da[ip],
                                  actso3->miv_rob->vicstnn.a.da[ip],
                                  actso3->miv_rob->bacstsn.a.da[ip]);
  }

  /*--------------------------------------------------------------------*/
  /* total strain eps_{n+1} at t_{n+1} */
  if (actso3->kintype == so3_geo_lin)
  {
    /* linear (engineering) strain vector */
    so3_mv6_v_ass(gds->stnengv, stntotn);
  } 
  else if (actso3->kintype == so3_total_lagr)
  {
    /* Green-Lagrange strain vector */
    so3_mv6_v_ass(gds->stnglv, stntotn);
  }
  else
  {
    dserror("Cannot digest chosen type of spatial kinematics\n");
  }

  /*--------------------------------------------------------------------*/
  /* temperature and thermal strain */
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

  /*--------------------------------------------------------------------*/
  /* viscous strain eps_{n+1}^{i} at t_{n+1} */
  so3_mv6_v_ass(actso3->miv_rob->vicstnn.a.da[ip], stnvscn);

  /*--------------------------------------------------------------------*/
  /* elastic strain at t_{n+c_i} */
  {
    INT istr;  /* strain index */
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      stnela[istr] = stntotn[istr] - stnvscn[istr] - stnthr[istr];
    }
  }

  /*--------------------------------------------------------------------*/
  /* elasticity tensor */
  so3_mat_robinson_elmat(mat_robin, tem, cmat);  /* cmat == kee */

  /*--------------------------------------------------------------------*/
  /* tangents of stress equation */
  /* so3_mv6_ass(cmat, kee); */
  so3_mv6_m_assscl(-1.0, cmat, kev);
  so3_mv6_m_id(kea);

  /*--------------------------------------------------------------------*/
  /* stress sig_{n+1}^i at t_{n+1} */
  so3_mv6_v_assmvp(cmat, stnela, stress);

  /*--------------------------------------------------------------------*/
  /* deviatoric stress s_{n+1}^i at t_{n+1} */
  so3_mv6_v_dev(stress, devstsn);

  /*--------------------------------------------------------------------*/
  /* overstress Sig_{n+1}^i = s_{n+1}^i - al_{n+1}^i */
  so3_mv6_v_sub(devstsn, actso3->miv_rob->bacstsn.a.da[ip], stsovrn);

  /*--------------------------------------------------------------------*/
  /* residual of visc. strain eps_{n+1}
   * and its conistent tangent */
  so3_mat_robinson_be_rvscstn(ele, mat_robin, dt, tem,
                              actso3->miv_rob->vicstn.a.da[ip],
                              actso3->miv_rob->vicstnn.a.da[ip],
                              devstsn, stsovrn,
                              &(actso3->miv_rob->vscstns.a.iv[ip]),
                              vscstnr, kve, kvv, kva);

  /*--------------------------------------------------------------------*/
  /* residual of back stress al_{n+1}
   * and its consistent tangent */
  so3_mat_robinson_be_rbcksts(ele, mat_robin, dt, tem,
                              actso3->miv_rob->vicstn.a.da[ip],
                              actso3->miv_rob->vicstnn.a.da[ip],
                              devstsn,
                              actso3->miv_rob->bacsts.a.da[ip],
                              actso3->miv_rob->bacstsn.a.da[ip],
                              &(actso3->miv_rob->bckstss.a.iv[ip]),
                              bckstsr, kae, kav, kaa);

  /*--------------------------------------------------------------------*/
  /* build reduced stress and tangent
   * ==> static condensation */
  so3_mat_robinson_be_red(stress, cmat, kev, kea,
                          vscstnr, kve, kvv, kva,
                          bckstsr, kae, kav, kaa,
                          actso3->miv_rob->kvarva.a.da[ip],
                          actso3->miv_rob->kvakvae.a.da[ip]);



  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Residual of BE-discretised viscous strain rate at Gauss point
\param   ele          ELEMENT*      (i)   element
\param   mat_robin    VP_ROBINSON*  (i)   element Robinson material
\param   dt           DOUBLE        (i)   time step size
\param   tmpr         DOUBLE        (i)   temperature
\param   vscstn       DOUBLE[]      (i)   viscous strain at t_n
\param   vscstnn      DOUBLE[]      (i)   viscous strain at t_{n+1}
\param   devstsn      DOUBLE[]      (i)   stress deviator at t_{n+1}
\param   ovrstsn      DOUBLE[]      (i)   over stress at t_{n+1}
\param   vscstns      INT*          (io)  viscous strain rate mode
\param   vscstnr      DOUBLE[]      (o)   viscous strain rate
\param   kve          DOUBLE[][]    (o)
\param   kvv          DOUBLE[][]    (o)
\param   kva          DOUBLE[][]    (o)
\author bborn
\date 04/07
*/
void so3_mat_robinson_be_rvscstn(ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const DOUBLE dt,
                                 const DOUBLE tmpr,
                                 const DOUBLE vscstn[NUMSTR_SOLID3],
                                 const DOUBLE vscstnn[NUMSTR_SOLID3],
                                 const DOUBLE devstsn[NUMSTR_SOLID3],
                                 const DOUBLE ovrstsn[NUMSTR_SOLID3],
                                 INT* vscstns,
                                 DOUBLE vscstnr[NUMSTR_SOLID3],
                                 DOUBLE kve[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kvv[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kva[NUMSTR_SOLID3][NUMSTR_SOLID3])

{
  DOUBLE j2;  /* 'J_2' */
  DOUBLE shrthrshld;  /* shear threshold 'K^2' */
  DOUBLE aa;  /* hardening factor 'A' */
  DOUBLE nn;  /* hardening exponent 'N' */
  DOUBLE ff;  /* 'F' */
  DOUBLE ss;  /* ss = 1/2 * s : Sig */

  DOUBLE dvscstn[NUMSTR_SOLID3];  /* viscous strain rate */
  DOUBLE kvs[NUMSTR_SOLID3][NUMSTR_SOLID3];  /* d eps^v/d Sig */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_rvscstn");
#endif

  /*--------------------------------------------------------------------*/
  /* preliminaries */
  /* J_2 = 1/2 * Sig : Sig  with Sig...overstress */
  so3_mv6_v_dblctr(ovrstsn, ovrstsn, &j2);
  j2 *= 0.5;
  /* Bingham-Prager shear stress threshold at current temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld, tmpr, &shrthrshld);
  /* F = (J_2 - K^2)/K_2 */
  if (fabs(shrthrshld) <= EPS10)
  {
    ff = -1.0;
    dserror("Division by zero: Shear threshold very close to zero");
  }
  else
  {
    ff = (j2 - shrthrshld)/shrthrshld;
  }
  /* ss = 1/2 * s : Sig  with  Sig...overstress, s...deviat.stress */
  so3_mv6_v_dblctr(ovrstsn, devstsn, &ss);
  ss *= 0.5;
  /* hardening factor aa */
  aa = mat_robin->hrdn_fact;
  /* hardening exponent nn */
  nn = mat_robin->hrdn_expo;

  /*--------------------------------------------------------------------*/
  /* determine mode
   * The mode is determined for every time step once. The iteration
   * sticks to the initially determined mode. */
  if (*vscstns == so3_mat_robinson_state_vague)
  {
    /* viscous mode */
    if ( (ff > 0.0) && (ss > 0.0) )
    {
      *vscstns = so3_mat_robinson_state_inelastic;
    }
    /* elastic mode */
    else
    {
      *vscstns = so3_mat_robinson_state_elastic;
    }
  }

  /*--------------------------------------------------------------------*/
  /* residual of viscous strain rate at t_{n+1} */
  if (*vscstns == so3_mat_robinson_state_elastic)
  {
    so3_mv6_v_zero(dvscstn);
  }
  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    DOUBLE fct = aa * pow(ff, nn) / sqrt(j2);
    so3_mv6_v2_assscl(fct, ovrstsn, dvscstn);
  }
  else
  {
    dserror("Oh, no");
  }

  /*--------------------------------------------------------------------*/
  /* viscous residual
   *    res_{n+1}^v = (eps_{n+1}^v - eps_n^v)/dt
   *                - deps_{n+1}^v (viscous strain rate) */
  {
    INT istr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      vscstnr[istr] = (vscstnn[istr] - vscstn[istr])/dt - dvscstn[istr];
    }
  }

  /*--------------------------------------------------------------------*/
  /* derivative of viscous residual with respect to over stress \Sig */
  if (*vscstns == so3_mat_robinson_state_elastic)
  {
    so3_mv6_m_zero(kvs);
  }
  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    const DOUBLE rhssign = -1.0;
    /* kvs = \frac{\pd res^v}{\pd Sig} */
    DOUBLE fctu = rhssign * aa * pow(ff, nn) / sqrt(j2);
    so3_mv6_m_idscl(fctu, kvs);
    /* kvs = \frac{\pd res^v}{\pd Sig} */
    DOUBLE fcto
      = rhssign * nn * aa * pow(ff, nn-1.0) / (sqrt(j2) * shrthrshld)
      - rhssign * aa * pow(ff, nn) / (2.0 * pow(j2, 1.5));
    so3_mv6_m_upddydscl(fcto, ovrstsn, ovrstsn, kvs);
    /* multiply by 2 last 3 rows to conform with definition
     * of strain vectors */
    so3_mv6_m2_updmtom2(kvs);
  }

  /*--------------------------------------------------------------------*/
  /* derivative with respect to total strain eps
   *    kve = \frac{\pd res_{n+1}^v}{\pd \eps_{n+1}}|^i */
  if (*vscstns == so3_mat_robinson_state_elastic)
  {
    so3_mv6_m_zero(kve);
  }
  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    DOUBLE kse[NUMSTR_SOLID3][NUMSTR_SOLID3];
    so3_mat_robinson_elmat(mat_robin, tmpr, kse);
    DOUBLE iv[NUMSTR_SOLID3];
    so3_mv6_v_id(iv);
    DOUBLE civ[NUMSTR_SOLID3];
    so3_mv6_v_assmvp(kse, iv, civ);
    so3_mv6_m_upddydscl(-1.0/3.0, iv, civ, kse);
    /* kve = kvs . kse */
    so3_mv6_m_mprd(kvs, kse, kve);
  }

  /*--------------------------------------------------------------------*/
  /* derivative with respect to visc. strain eps^v */
  if (*vscstns == so3_mat_robinson_state_elastic)
  {
    so3_mv6_m_idscl(1.0/dt, kvv);
  }
  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    /* matrix  Im - 1/3 Iv . Iv^T */
    DOUBLE fm[NUMSTR_SOLID3][NUMSTR_SOLID3];
    so3_mv6_m_id(fm);
    so3_mv6_m_updonescl(-1.0/3.0, fm);
    /* elasticity matrix */
    DOUBLE cm[NUMSTR_SOLID3][NUMSTR_SOLID3];
    so3_mat_robinson_elmat(mat_robin, tmpr, cm);
    /* matrix  ksv = factm . (-cmat) */
    DOUBLE ksv[NUMSTR_SOLID3][NUMSTR_SOLID3];
    so3_mv6_m_mprdscl(-1.0, fm, cm, ksv);
    /* matrix  kvv part I */
    so3_mv6_m_idscl(1.0/dt, kvv);
    /* matrix  kvv part II */
    so3_mv6_m_updmprd(kvs, ksv, kvv);
  }

  /*--------------------------------------------------------------------*/
  /* derivative with respect to back stress al */
  if (*vscstns == so3_mat_robinson_state_elastic)
  {
    so3_mv6_m_zero(kva);
  }
  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    so3_mv6_m_assscl(-1.0, kvs, kva);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Residual of BE-discretised back stress rate flow rule at Gauss point
\param   mat_robin    VP_ROBINSON*  (i)   element Robinson material
\param   tmpr         DOUBLE        (i)   temperature
\param   stsdev       DOUBLE[]      (i)   deviatoric stress
\param   stsbck       DOUBLE[]      (i)   back stress
\param   dstnvsc      DOUBLE[]      (i)   viscous strain rate
\param   dstsbck      DOUBLE[]      (o)   back stress rate
\author bborn
\date 04/07
*/
void so3_mat_robinson_be_rbcksts(ELEMENT* ele,
                                 const VP_ROBINSON* mat_robin,
                                 const DOUBLE dt,
                                 const DOUBLE tmpr,
                                 const DOUBLE vscstn[NUMSTR_SOLID3],
                                 const DOUBLE vscstnn[NUMSTR_SOLID3],
                                 const DOUBLE devstsn[NUMSTR_SOLID3],
                                 const DOUBLE bacsts[NUMSTR_SOLID3],
                                 const DOUBLE bacstsn[NUMSTR_SOLID3],
                                 INT* bckstss,
                                 DOUBLE bckstsr[NUMSTR_SOLID3],
                                 DOUBLE kae[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kav[NUMSTR_SOLID3][NUMSTR_SOLID3],
                                 DOUBLE kaa[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
/*   const SOLID3* actso3 = ele->e.so3; */

  DOUBLE i2;  /* 'I_2' */
  DOUBLE tem0;  /* activation temperature */
  DOUBLE kk0;  /* shear threshold 'K_0^2' */
  DOUBLE hh;  /* 'H' */
  DOUBLE beta;  /* 'beta' */
  DOUBLE mm;  /* 'm' */
  DOUBLE rr0;  /* recovery factor 'R_0' */
  DOUBLE rr;  /* recovery term 'R' */
  DOUBLE gg0;  /* 'G_0' */
  DOUBLE gg;  /* 'G' */
  DOUBLE sa;  /* */

  DOUBLE vscstnd05[NUMSTR_SOLID3];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_rbcksts");
#endif

  /*--------------------------------------------------------------------*/
  /* preliminaries */
  /* I_2 = 1/2 * Alpha : Alpha  with Alpha...back stress */
  so3_mv6_v_dblctr(bacstsn, bacstsn, &i2);
  i2 *= 0.5;
  /* activation temperature */
  tem0 = mat_robin->actv_tmpr;
  /* Bingham-Prager shear stress threshold at activation temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld, tem0, &kk0);
  /* 'H' at current temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->h, tmpr, &hh);
  /* 'beta' at current temperature */
  so3_mat_robinson_prmbytmpr(mat_robin->beta, tmpr, &beta);
  /* 'm' */
  mm = mat_robin->m;
  /* recovery factor 'R_0' */
  so3_mat_robinson_prmbytmpr(mat_robin->rcvry, tmpr, &rr0);
  /* 'R' */
  rr = rr0 * exp(mat_robin->actv_ergy*(tmpr-tem0)/(tmpr*tem0));
  /* 'G_0' */
  gg0 = mat_robin->g0;
  /* G = I_2/K_0^2 */
  if (fabs(kk0) <= EPS10)
  {
    dserror("Division by zero: Shear threshold very close to zero");
  }
  else
  {
    gg = sqrt(i2/kk0);
  }
  /* ss = 1/2 * s : Alpha  with  Alpha...backstress, s...deviat.stress */
  so3_mv6_v_dblctr(bacstsn, devstsn, &sa);
  sa *= 0.5;

  /*--------------------------------------------------------------------*/
  /* determine mode */
  if (*bckstss == so3_mat_robinson_state_vague)
  {
    /* viscous mode */
    if ( (gg > gg0) && (sa > 0.0) )
    {
      *bckstss = so3_mat_robinson_state_inelastic;
    }
    /* elastic mode */
    else
    {
      *bckstss = so3_mat_robinson_state_elastic;
    }
  }

  /*--------------------------------------------------------------------*/
  /* difference of current and last viscous strains
   *    \incr \eps^v = \eps_{n+1}^v - \eps_{n}^v
   * with halved entries to conform with stress vectors */
  so3_mv6_v_sub(vscstnn, vscstn, vscstnd05);
  so3_mv6_v05_updvtov05(vscstnd05);

  /*--------------------------------------------------------------------*/
  /* residual of back stress rate */
  if (*bckstss == so3_mat_robinson_state_elastic)
  {
    DOUBLE fctv = hh / pow(gg0, beta);
    DOUBLE fcta;
    if (sqrt(i2) < EPS10)
    {
      /* sqrt(i2) := 1.0e6 assures units are OK */
      fcta = rr * pow(gg0, (mm-beta)) / 1.0e6;
    }
    else
    {
      fcta = rr * pow(gg0, (mm-beta)) / sqrt(i2);
    }
    INT istr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      bckstsr[istr] = ( bacstsn[istr] - bacsts[istr]
                        - fctv * vscstnd05[istr]
                        + dt * fcta * bacstsn[istr] ) / dt;
    }
  }
  else if (*bckstss == so3_mat_robinson_state_inelastic)
  {
    DOUBLE fctv = hh / pow(gg, beta);
    DOUBLE fcta = rr * pow(gg, (mm-beta)) / sqrt(i2);
    INT istr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      bckstsr[istr] = ( bacstsn[istr] - bacsts[istr]
                        - fctv * vscstnd05[istr]
                        + dt * fcta * bacstsn[istr] ) / dt;
    }
  }

  /*--------------------------------------------------------------------*/
  /* derivative of back stress residual with respect to total strains */
  if ( (*bckstss == so3_mat_robinson_state_elastic) 
       || (*bckstss == so3_mat_robinson_state_inelastic) ) 
  {
    so3_mv6_m_zero(kae);
  } 

  /*--------------------------------------------------------------------*/
  /* derivative of back stress residual with respect to viscous strains
   *    kav = \frac{\pd res_{n+1}^al}{\pd eps_{n+1}^v} */
  if (*bckstss == so3_mat_robinson_state_elastic)
  {
    DOUBLE fctv = -hh / (pow(gg0, beta) * dt);
    so3_mv6_m05_idscl(fctv, kav);
  }
  else if (*bckstss == so3_mat_robinson_state_inelastic)
  {
    DOUBLE fctv = -hh / (pow(gg, beta) * dt);
    so3_mv6_m05_idscl(fctv, kav);
  }
  

  /*--------------------------------------------------------------------*/
  /* derivative of back stress residual with respect to back stress
   *    kaa = \frac{\pd res_{n+1}^al}{\pd al_{n+1}} */
  if (*bckstss == so3_mat_robinson_state_elastic)
  {
    DOUBLE fctu = 1.0/dt;
    DOUBLE fcta;
    if (sqrt(i2) < EPS10)
    {
      DOUBLE ii2 = 1.0e12;  /* sqrt(i2) := 1.0e6 assures units are OK */
      fctu += rr * pow(gg0,(mm-beta)) / sqrt(ii2);
      fcta = -rr * pow(gg0,(mm-beta)) / (2.0 * pow(ii2, 1.5));
    }
    else
    {
      fctu += rr * pow(gg0,(mm-beta)) / sqrt(i2);
      fcta = -rr * pow(gg0,(mm-beta)) / (2.0 * pow(i2, 1.5));
    }
    so3_mv6_m_idscl(fctu, kaa);
    so3_mv6_m_upddydscl(fcta, bacstsn, bacstsn, kaa);
  }
  else if (*bckstss == so3_mat_robinson_state_inelastic)
  {
    DOUBLE fctu = 1.0/dt  +  rr * pow(gg,(mm-beta)) / sqrt(i2);
    DOUBLE fctv = beta * hh / ( pow(gg,(beta+1.0)) * dt * kk0 );
    DOUBLE fcta = rr * (mm-beta) * pow(gg,(mm-beta-1.0)) / (sqrt(i2) * kk0)
                - rr * pow(gg,(mm-beta)) / (2.0 * pow(i2, 1.5));
    so3_mv6_m_idscl(fctu, kaa);
    so3_mv6_m_upddydscl(fctv, vscstnd05, bacstsn, kaa);
    so3_mv6_m_upddydscl(fcta, bacstsn, bacstsn, kaa);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}



/*======================================================================*/
/*!
\brief Reduce (statically condense) system in eps,eps^v,al to purely eps

The linearised stress and internal residuals are

      [ sig   ]         [ sig    ]^i
  Lin [ res^v ]       = [ res^v  ]
      [ res^al]_{n+1}   [ res^al ]_{n+1}

                           [ kee  kev  kea ]^i  [ iinc eps   ]^i
                        +  [ kve  kvv  kva ]    [ iinc eps^v ]
                           [ kae  kav  kaa ]    [ iinc al    ]_{n+1}

                        [ sig ]
                      = [  0  ]  on every element (e)
                        [  0  ]  and at each Gauss point <g>

Due to the fact that the internal residuals (the BE-discretised evolution
laws of the viscous strain and the back stress) are C^{-1}-continuous
across element boundaries. We can statically condense this system.
The iterative increments inc eps^v and inc al are expressed in inc eps.
We achieve

  [ iinc eps^v ]   [ kvv  kva ]^{-1} (   [ res^v  ]   [ kve ]                )
  [ iinc al    ] = [ kav  kaa ]      ( - [ res^al ] - [ kae ] * [ iinc eps ] )

thus

                                     [ kvv  kva ]^{-1} [ res^v  ]^i
  sig_red^i = sig^i - [ kev  kea ]^i [ kav  kaa ]      [ res^al ]

and
                                     [ kvv  kva ]^{-1} [ kve ]^i
  kee_red^i = kee^i - [ kev  kea ]^i [ kav  kaa ]      [ kae ]

\param  stress      DOUBLE[]    (io)   stress vector
\param  cmat        DOUBLE[][]  (io)   constitutive tensor
\param  kev         DOUBLE[][]  (i)    \frac{\pd \sig}{\pd \eps^v}
\param  kea         DOUBLE[][]  (i)    \frac{\pd \sig}{\pd \al}
\param  vscstnr     DOUBLE[]    (i)    visc. strain residual
\param  kve         DOUBLE[][]  (i)    \frac{\pd res^v}{\pd \eps}
\param  kvv         DOUBLE[][]  (i)    \frac{\pd res^v}{\pd \eps^v}
\param  kva         DOUBLE[][]  (i)    \frac{\pd res^v}{\pd \al}
\param  kae         DOUBLE[][]  (i)    \frac{\pd res^al}{\pd \eps}
\param  kav         DOUBLE[][]  (i)    \frac{\pd res^al}{\pd \eps^v}
\param  kaa         DOUBLE[][]  (i)    \frac{\pd res^al}{\pd \al}
\param  kvarva      DOUBLE*     (o)    cf. solid3.h
\param  kvakvae     DOUBLE*     (o)    cf. solid3.h

\author bborn
\date 04/07
*/
void so3_mat_robinson_be_red(DOUBLE stress[NUMSTR_SOLID3], 
                             DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kev[NUMSTR_SOLID3][NUMSTR_SOLID3], 
                             DOUBLE kea[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             const DOUBLE vscstnr[NUMSTR_SOLID3],
                             DOUBLE kve[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kvv[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kva[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             const DOUBLE bckstsr[NUMSTR_SOLID3],
                             DOUBLE kae[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE kav[NUMSTR_SOLID3][NUMSTR_SOLID3], 
                             DOUBLE kaa[NUMSTR_SOLID3][NUMSTR_SOLID3],
                             DOUBLE* kvarva,
                             DOUBLE* kvakvae)
{
  INT one = 1;
  INT numstr = NUMSTR_SOLID3;
  INT numstr_2 = 2*NUMSTR_SOLID3;
  INT info;
  CHAR trans;

  /*             [ kvv  kva ]
   * kvvvaavaa = [ kav  kaa ]  and its inverse after factorisation */
  DOUBLE* kvvvaavaa = (DOUBLE*) CCACALLOC(numstr_2*numstr_2, sizeof(DOUBLE));
  /* pivot indices */
  INT* pivot_indices = (INT*) CCACALLOC(numstr_2, sizeof(INT));
  /* kevea = [ kev  kea ] */
  DOUBLE* kevea = (DOUBLE*) CCACALLOC(numstr*numstr_2, sizeof(DOUBLE));

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_red");
#endif

  /*--------------------------------------------------------------------*/
  /* build tangent and right hand side to reduce */
  {
    INT i, j, k;
    /* first NUMSTR_SOLID3 rows */
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      /* first NUMSTR_SOLID3 columns */
      for(j=0; j<NUMSTR_SOLID3; j++)
      {
        /* tangent */
        kvvvaavaa[i+j*numstr_2] = kvv[i][j];
        /* tangent */
        kvakvae[i+j*numstr_2] = kve[i][j];
      }
      /* second NUMSTR_SOLID3 columns */
      for(k=0; k<NUMSTR_SOLID3; k++)
      {
        /* tangent */
        kvvvaavaa[i+(NUMSTR_SOLID3+k)*numstr_2] = kva[i][k];
      }
      /* residual vector */
      kvarva[i] = vscstnr[i];
    }
    /* second NUMSTR_SOLID3 rows */
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      /* first NUMSTR_SOLID3 columns */
      for(j=0; j<NUMSTR_SOLID3; j++)
      {
        /* tangent */
        kvvvaavaa[NUMSTR_SOLID3+i+j*numstr_2] = kav[i][j];
        /* tangent */
        kvakvae[NUMSTR_SOLID3+i+j*numstr_2] = kae[i][j];
      }
      /* second NUMSTR_SOLID3 columns */
      for(k=0; k<NUMSTR_SOLID3; k++)
      {
        /* tangent */
        kvvvaavaa[NUMSTR_SOLID3+i+(NUMSTR_SOLID3+k)*numstr_2] = kaa[i][k];
      }
      /* residual vector */
      kvarva[NUMSTR_SOLID3+i] = bckstsr[i];
    }
  }

  /*--------------------------------------------------------------------*/
  /* factorise kvvvaavaa */
#ifndef AZTEC_PACKAGE
  /* factorise kvvvaavaa */
  dgetrf(&numstr_2,  /* number of rows of matrix */
         &numstr_2,  /* number of columns of matrix */
         kvvvaavaa,  /* matrix stored in Fortranesque vector (column-major) */
         &numstr_2,  /* leading dimension of matrix */
         pivot_indices,
         &info);  /* info: 0=OK, -1=failed */
  if (info != 0)
  {
    dserror("Lapack factorisation failed");
  }
  /* back substitution of residuals */
  trans = 'N';
  dgetrs(&trans,  /* type: 'N'=normal */
         &numstr_2,  /* number of columns of matrix */
         &one,  /* number of RHSs */
         kvvvaavaa,  /* factorised matrix */
         &numstr_2,  /* leading dimension of matrix */
         pivot_indices,  /* pivot indices */
         kvarva,  /* RHS vector (input/output) solution */
         &numstr_2,  /* leading dimension of RHS vector */
         &info);
  if (info != 0)
  {
    dserror("Lapack back substitution failed (vector)");
  }
  /* back substitution of tangent */
  trans = 'N';
  dgetrs(&trans,  /* type: 'N'=normal */
         &numstr_2,  /* number of columns of matrix */
         &numstr,  /* number of RHSs */
         kvvvaavaa,  /* factorised matrix */
         &numstr_2,  /* leading dimension of matrix */
         pivot_indices,  /* pivot indices */
         kvakvae,  /* RHS Fortranesque matrix (input/output) solution */
         &numstr_2,  /* leading dimension of RHS matrix */
         &info);
  if (info != 0)
  {
    dserror("Lapack back substitution failed (matrix)");
  }
#else
  dserror("solver Lapack conflicts with compilation with -DAZTEC_PACKAGE");
#endif

  /*--------------------------------------------------------------------*/
  /* reduce stress vector */
  {
    INT i;
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      DOUBLE rcsum = 0.0;
      INT j;
      for (j=0; j<numstr_2; j++)
      {
        rcsum += kevea[i*numstr_2+j] * kvarva[j];
      }
      stress[i] -= rcsum;
    }
  }
  
  /*--------------------------------------------------------------------*/
  /* reduce tangent */
  {
    INT i;
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      INT j;
      for (j=0; j<NUMSTR_SOLID3; j++)
      {
        DOUBLE rcsum = 0.0;
        INT k;
        for (k=0; k<numstr_2; k++)
        {
          rcsum += kevea[i*numstr_2+k] * kvakvae[j*numstr_2+k];
        }
        cmat[i][j] -= rcsum;
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* clean up */
  CCAFREE(kvvvaavaa);
  CCAFREE(pivot_indices);
  CCAFREE(kevea);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Iterative update Robinson's internal material variables
\author bborn
\date 05/07
*/
void so3_mat_robinson_be_mivupditr(const CONTAINER* container,
                                   ELEMENT* ele,
                                   DOUBLE gpderiv[MAXNOD_SOLID3][NDIM_SOLID3],
                                   SO3_GEODEFSTR* gds,
                                   const DOUBLE* kvarva,
                                   const DOUBLE* kvakvae,
                                   DOUBLE vscstnn[NUMSTR_SOLID3],
                                   DOUBLE bckstsn[NUMSTR_SOLID3])
{
  const INT numstr_2 = 2*NUMSTR_SOLID3;
  const ARRAY_POSITION_SOL *isol
    = &(field[genprob.numsf].dis[container->disnum_s].ipos.isol);
  const ARRAY_POSITION_SOLRES *isolres
    = &(field[genprob.numsf].dis[container->disnum_s].ipos.isolres);
  const INT nelenod = ele->numnp;  /* number of element nodes */
  const INT neledof = NUMDOF_SOLID3 * nelenod;  /* total elem. DOFs */
  DOUBLE ediso[MAXNOD_SOLID3][NDIM_SOLID3];  /* current elem. displ. at
                                              * last iteration */
  DOUBLE disii[MAXDOF_SOLID3];  /* iterative displacements */
  DOUBLE disgrdv[NUMDFGR_SOLID3];  /* displacement gradient vector */
  DOUBLE defgrd[NDIM_SOLID3][NDIM_SOLID3];  /* deformation gradient tensor */
  DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3];  /* B-operator */
  DOUBLE epsii[NUMSTR_SOLID3];  /* iterative/residual total strain */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_mivupditr");
#endif

  /*--------------------------------------------------------------------*/
  /* element nodal displacements are last iteration dis_{n+1}^<i-1> */
  {
    INT inod;
    for (inod=0; inod<nelenod; inod++)
    {
      NODE* actnode = ele->node[inod];
      INT jdim;
      for (jdim=0; jdim<NDIM_SOLID3; jdim++)
      {
        /* reconstruct displacement increment \inc\D_{n+1}^<i-1> 
         * of previous iteration <i-1> */
        ediso[inod][jdim] = actnode->sol.a.da[isol->disn][jdim]
          - actnode->sol_residual.a.da[isolres->disres][jdim];
        /* build iterative/residual displacements \iinc\D_{n+1}^<i> */
        disii[inod*NDIM_SOLID3+jdim] 
          = actnode->sol_residual.a.da[isolres->disres][jdim];
      }
    }
  }
    
  /*--------------------------------------------------------------------*/
  /* deformation tensor and displacement gradient */
  so3_def_grad(nelenod, ediso, gpderiv, gds->xji, 
               disgrdv, defgrd);
  
  /*--------------------------------------------------------------------*/
  /* calculate B-operator
   * bop differs depending on geometrically linearity/non-linearity */
  if (ele->e.so3->kintype == so3_geo_lin)
  {
    /* loop over element nodes */
    INT inod;
    for (inod=0; inod<nelenod; inod++)
    {
      /* nodally stored B-operator */
      /*           [ ... | N_{,1}^k | ... ]
       *    Bn^T = [ ... | N_{,2}^k | ... ]
       *           [ ... | N_{,3}^k | ... ]
       */
      DOUBLE N_X = gds->bopn[inod][0];
      DOUBLE N_Y = gds->bopn[inod][1];
      DOUBLE N_Z = gds->bopn[inod][2];
      /* address node-column in operator matrix */
      INT idof_X = inod * NUMDOF_SOLID3;
      INT idof_Y = idof_X + 1;
      INT idof_Z = idof_X + 2;
      /* linear B-operator */
      /* 
       *     [ ... | N_{,X}^k                       | ... ]
       *     [ ... |            N_{,Y}^k            | ... ]
       *     [ ... |                       N_{,Z}^k | ... ]
       * B = [ ~~~   ~~~~~~~    ~~~~~~~~   ~~~~~~~~   ~~~ ]
       *     [ ... | N_{,Y}^k   N_{,X}^k            | ... ]
       *     [ ... |            N_{,Z}^k   N_{,Y}^k | ... ]
       *     [ ... | N_{,Z}^k              N_{,X}^k | ....]
       */
      bop[0][idof_X] = N_X;  bop[0][idof_Y] = 0.0;  bop[0][idof_Z] = 0.0;
      bop[1][idof_X] = 0.0;  bop[1][idof_Y] = N_Y;  bop[1][idof_Z] = 0.0;
      bop[2][idof_X] = 0.0;  bop[2][idof_Y] = 0.0;  bop[2][idof_Z] = N_Z;
      /* ~~~ */
      bop[3][idof_X] = N_Y;  bop[3][idof_Y] = N_X;  bop[3][idof_Z] = 0.0;
      bop[4][idof_X] = 0.0;  bop[4][idof_Y] = N_Z;  bop[4][idof_Z] = N_Y;
      bop[5][idof_X] = N_Z;  bop[5][idof_Y] = 0.0;  bop[5][idof_Z] = N_X;   
    } 
  }
  else if (ele->e.so3->kintype == so3_total_lagr)
  {
    /* loop over element nodes */
    INT inod;
    for (inod=0; inod<nelenod; inod++)
    {
      /* nodally stored B-operator */
      /*           [ ... | N_{,1}^k | ... ]
       *    Bn^T = [ ... | N_{,2}^k | ... ]
       *           [ ... | N_{,3}^k | ... ]
       */
      DOUBLE N_X = gds->bopn[inod][0];
      DOUBLE N_Y = gds->bopn[inod][1];
      DOUBLE N_Z = gds->bopn[inod][2];
      /* address node-column in operator matrix */
      INT idof_X = inod * NUMDOF_SOLID3;
      INT idof_Y = idof_X + 1;
      INT idof_Z = idof_X + 2;
      /* non-linear B-operator (may so be called, meaning
       * of B-operator is not so sharp in the non-linear realm) */
      /*
       * B = F . Bl
       *
       *      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
       *      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
       *      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
       * B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
       *      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
       *      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
       *      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
       *      [                                                         ]
       *      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
       *      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
       *      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
       *      [                                                         ]
       *      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
       *      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
       *      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
       */
      bop[0][idof_X] = defgrd[0][0]*N_X;
      bop[0][idof_Y] = defgrd[1][0]*N_X;
      bop[0][idof_Z] = defgrd[2][0]*N_X;
      bop[1][idof_X] = defgrd[0][1]*N_Y;
      bop[1][idof_Y] = defgrd[1][1]*N_Y;
      bop[1][idof_Z] = defgrd[2][1]*N_Y;
      bop[2][idof_X] = defgrd[0][2]*N_Z; 
      bop[2][idof_Y] = defgrd[1][2]*N_Z; 
      bop[2][idof_Z] = defgrd[2][2]*N_Z; 
      /* ~~~ */
      bop[3][idof_X] = defgrd[0][0]*N_Y + defgrd[0][1]*N_X;
      bop[3][idof_Y] = defgrd[1][0]*N_Y + defgrd[1][1]*N_X;
      bop[3][idof_Z] = defgrd[2][0]*N_Y + defgrd[2][1]*N_X;
      bop[4][idof_X] = defgrd[0][1]*N_Z + defgrd[0][2]*N_Y;
      bop[4][idof_Y] = defgrd[1][1]*N_Z + defgrd[1][2]*N_Y;
      bop[4][idof_Z] = defgrd[2][1]*N_Z + defgrd[2][2]*N_Y;
      bop[5][idof_X] = defgrd[0][2]*N_X + defgrd[0][0]*N_Z;
      bop[5][idof_Y] = defgrd[1][2]*N_X + defgrd[1][0]*N_Z;
      bop[5][idof_Z] = defgrd[2][2]*N_X + defgrd[2][0]*N_Z;
    }
  }
  else
  {
    dserror("Cannot compute B-operator for chosen kinematics.");
  }

  /*--------------------------------------------------------------------*/
  /* iterative total strain increment \iinc\eps
   *    \iinc\eps_{n+1}^<i-1> = \B . \iinc\D_{n+1}^<i-1> */
  {
    INT istr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      DOUBLE rcsum = 0.0;
      INT idof;
      for (idof=0; idof<neledof; idof++)
      {
        rcsum += bop[istr][idof] * disii[idof];
      }
      epsii[istr] = rcsum;
    }
  }

  /*--------------------------------------------------------------------*/
  /* update visc. strain */
  {
    INT istr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      DOUBLE rcsum = kvarva[istr];  /* visc. resid. contribution */
      INT jstr;
      for (jstr=0; jstr<NUMSTR_SOLID3; jstr++)
      {
        rcsum += kvakvae[istr+jstr*numstr_2] * epsii[jstr];
      }
      vscstnn[istr] -= rcsum;
    }
  }

  /*--------------------------------------------------------------------*/
  /* update back stress */
  {
    INT istr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      DOUBLE rcsum = kvarva[NUMSTR_SOLID3+istr];  /* back str. resid. contr */
      INT jstr;
      for (jstr=0; jstr<NUMSTR_SOLID3; jstr++)
      {
        rcsum += kvakvae[NUMSTR_SOLID3+istr+jstr*numstr_2] * epsii[jstr];
      }
      bckstsn[istr] -= rcsum;
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
\brief Incremental update Robinson's internal material variables 
\param   ele          ELEMENT*        (io)   curr. elem.
\parm    mat_robin    VP_ROBINSON*    (i)    elem. mater.
\author bborn
\date 05/07
*/
void so3_mat_robinson_be_mivupd(ELEMENT* ele,
                                const VP_ROBINSON* mat_robin)
{
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  INT gptot = actso3->gptot;  /* total number of GPs in domain */
  INT ip;  /* total GP index */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_mivupd");
#endif

  /*--------------------------------------------------------------------*/
  /* update/reset MIVs */
  for (ip=0; ip<gptot; ip++)
  {
    /* reset visous strain mode */
    actso3->miv_rob->vscstns.a.iv[ip] = so3_mat_robinson_state_vague;
    /* update viscous strain
     *    eps_{n}^v := eps_{n+1}^v at every Gauss point <g> */
    so3_mv6_v_ass(actso3->miv_rob->vicstnn.a.da[ip],
                  actso3->miv_rob->vicstn.a.da[ip]);
    /* reset back stress mode */
    actso3->miv_rob->bckstss.a.iv[ip] = so3_mat_robinson_state_vague;
    /* update back stress 'Alpha'
     *    al_{n} := al_{n+1} at every Gauss point <g> */
    so3_mv6_v_ass(actso3->miv_rob->bacstsn.a.da[ip],
                  actso3->miv_rob->bacsts.a.da[ip]);
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
void so3_mat_robinson_be_stress(const CONTAINER* container,
                                const ELEMENT* ele,
                                const VP_ROBINSON* mat_robin,
                                const INT ip,
                                const SO3_GEODEFSTR* gds,
                                DOUBLE stress[NUMSTR_SOLID3],
                                DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  SOLID3* actso3 = ele->e.so3;  /* point to SOLID3 element bits */
  DOUBLE totstn[NUMSTR_SOLID3];  /* total strain */
  DOUBLE elastn[NUMSTR_SOLID3];  /* elastic strain */
  DOUBLE thrstn[NUMSTR_SOLID3];  /* thermal strain */
  DOUBLE vscstn[NUMSTR_SOLID3];  /* viscous strain */
  
  DOUBLE tem;  /* temperature */

  INT itsidyn = genprob.numfld;  /* index of TSI dynamics data */
  TSI_DYNAMIC* tsidyn = alldyn[itsidyn].tsidyn;  /* TSI dynamics data */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_robinson_be_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* total strain at t_{n} */
  if (actso3->kintype == so3_geo_lin)
  {
    /* linear (engineering) strain vector */
    so3_mv6_v_ass(gds->stnengv, totstn);
  } 
  else if (actso3->kintype == so3_total_lagr)
  {
    /* Green-Lagrange strain vector */
    so3_mv6_v_ass(gds->stnglv, totstn);
  }
  else
  {
    dserror("Cannot digest chosen type of spatial kinematic\n");
  }

  /*--------------------------------------------------------------------*/
  /* temperature and thermal strain */
  {
    /* temperature at Gauss point */
    so3_tsi_temper(container, ele,
                   gds->gpc[0], gds->gpc[1], gds->gpc[2], 
                   &tem);
    /* coefficient of linear thermal expansion */
    DOUBLE thermexpans = mat_robin->thermexpans;
    /* thermal strain vector */
    thrstn[0] = thermexpans * tem;  /* E_xx */
    thrstn[1] = thermexpans * tem;  /* E_yy */
    thrstn[2] = thermexpans * tem;  /* E_zz */
    thrstn[3] = 0.0;  /* E_xy */
    thrstn[4] = 0.0;  /* E_yz */
    thrstn[5] = 0.0;  /* E_zx */
  }

  /*--------------------------------------------------------------------*/
  /* viscous strain at t_{n} */
  so3_mv6_v_ass(actso3->miv_rob->vicstn.a.da[ip], vscstn);

  /*--------------------------------------------------------------------*/
  /* elastic strain at t_{n} */
  {
    INT istr;  /* strain index */
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      elastn[istr] = totstn[istr] - vscstn[istr] - thrstn[istr];
    }
  }

  /*--------------------------------------------------------------------*/
  /* elasticity tensor */
  so3_mat_robinson_elmat(mat_robin, tem, cmat);

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
        stresssum += cmat[istss][istrn] * elastn[istrn];
      }
      stress[istss] = stresssum;
    }
  }

  /*--------------------------------------------------------------------*/
  /* debug: */
#ifdef TESTROBIN_SOLID3
  if (ip == 0)
  {
    INT kstep = tsidyn->step;  /* current time step */
    DOUBLE stsdev[NUMSTR_SOLID3];
    /* stress deviator */
    so3_mv6_v_dev(stress, stsdev);
    /* overstress 'Sig' = 's' - 'alpha' */
    DOUBLE stsovr[NUMSTR_SOLID3];
    so3_mv6_v_sub(stsdev, actso3->miv_rob->bacsts.a.da[ip], stsovr);
    FILE* oo = allfiles.gnu;
    fprintf(oo, "# Time step %d Gauss point %d\n", kstep, ip);
    fprintf(oo, "# step strain-tot strain-ela strain-vis strain-thr"
            "  stress-tot stress-dev stress-back stress-over\n");
    INT i;
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      fprintf(oo, "%d %12.4e %12.4e %12.4e %12.4e  %12.4e %12.4e %12.4e %12.4e\n", 
              kstep,
              totstn[i], elastn[i], vscstn[i], thrstn[i],
              stress[i], 
              stsdev[i], 
              actso3->miv_rob->bacsts.a.da[ip][i], 
              stsovr[i]);
    }
    /* fprintf(oo, "\n"); */
  }
  if (ip == 7)
  {
    INT kstep = tsidyn->step;  /* current time step */
    DOUBLE stsdev[NUMSTR_SOLID3];
    /* stress deviator */
    so3_mv6_v_dev(stress, stsdev);
    /* overstress 'Sig' = 's' - 'alpha' */
    DOUBLE stsovr[NUMSTR_SOLID3];
    so3_mv6_v_sub(stsdev, actso3->miv_rob->bacsts.a.da[ip], stsovr);
    FILE* oo = allfiles.gnu;
    fprintf(oo, "# Time step %d Gauss point %d\n", kstep, ip);
    fprintf(oo, "# step strain-tot strain-ela strain-vis strain-thr"
            "  stress-tot stress-dev stress-back stress-over\n");
    INT i;
    for (i=0; i<NUMSTR_SOLID3; i++)
    {
      fprintf(oo, "%d %12.4e %12.4e %12.4e %12.4e  %12.4e %12.4e %12.4e %12.4e\n", 
              kstep,
              totstn[i], elastn[i], vscstn[i], thrstn[i],
              stress[i], stsdev[i], 
              actso3->miv_rob->bacsts.a.da[ip][i], stsovr[i]);
    }
    /* fprintf(oo, "\n"); */
  }
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif  /* end D_TSI */
#endif  /* end D_SOLID3 */
